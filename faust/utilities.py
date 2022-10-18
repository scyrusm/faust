import pandas as pd
import numpy as np
import os


def read_gpp_output(
        folders,
        chipfile=None,
        barcode2gene_dict=None,
        indices=['Construct Barcode', 'Construct IDs', 'UMI', 'Target Gene'],
        dropcols=[]):
    """

    Parameters
    ----------
    folders :
        
    chipfile :
         (Default value = None)
    barcode2gene_dict :
         (Default value = None)
    indices :
         (Default value = ['Construct Barcode')
    'Construct IDs' :
        
    'UMI' :
        
    'Target Gene'] :
        
    dropcols :
         (Default value = [])

    Returns
    -------

    """
    if chipfile is None and barcode2gene_dict is None:
        raise Exception(
            "One of barcode2gene_dict or chipfile must not be None")
    if barcode2gene_dict is None:
        chip = pd.read_table(chipfile)
        barcode2gene_dict = chip.set_index(
            'Barcode Sequence')['Gene Symbol'].to_dict()
    dfss = []
    for folder in folders:
        dfs = []
        for txt in [x for x in os.listdir(folder) if 'MATCH' not in x]:
            df = pd.read_table(folder + '/' + txt)
            umi = txt.split('-')[-1].replace('.txt', '')
            df['UMI'] = umi
            dfs.append(df)
        dfs = pd.concat(dfs)
        dfs['Target Gene'] = [
            barcode2gene_dict[x] for x in dfs['Construct Barcode']
        ]
        #for col in [x for x in dfs.columns if 'Tumor' in x or 'LN' in x]:
        #    dfs[col+' IsPos'] = dfs[col]>0
        #dfs['Input Mean'] = dfs[[x for x in dfs.columns if 'Input' in x]].mean(axis=1)
        #dfs = dfs[dfs['Input Mean']>0] # optional
        dfss.append(dfs.copy())
    dfss = [
        x[[y for y in x.columns if y not in dropcols]].set_index(indices)
        for x in dfss
    ]
    return sum(dfss).reset_index()


def get_summary_df(df,
                   controls,
                   inputs,
                   outputs,
                   input_type='single',
                   verbose=True,
                   count_threshold=1,
                   estimate_cells=True):
    """

    Parameters
    ----------
    df :
        
    controls :
        
    inputs :
        
    outputs :
        

    Returns
    -------

    """
    from tqdm import tqdm
    import pandas as pd
    from faust.utilities import nan_fdrcorrection_q
    from scipy.stats import mannwhitneyu

    if input_type not in ['single', 'matched']:
        raise Exception("input_type must be one of 'single', 'matched'")
    df = df.copy()
    if input_type == 'single':
        df['Input Sum'] = df[inputs].sum(axis=1)
        df = df[df['Input Sum'] > count_threshold]
        inputs = ['Input Sum'] * len(outputs)
    else:
        if len(outputs) != len(inputs):
            raise Exception(
                "outputs and inputs must be equal length when input_type is 'matched'"
            )
    output2input_dict = {}
    for output_col, input_col in zip(outputs, inputs):
        df[output_col] = np.where(df[output_col] > count_threshold,
                                  df[output_col], 0)
        df[input_col] = np.where(df[input_col] > count_threshold,
                                 df[input_col], 0)
        if verbose:
            print("Comparing input", input_col, "with output", output_col)
        df[output_col + '_odds_ratio'] = np.divide(df[output_col],
                                                   df[input_col])
        output2input_dict[output_col] = input_col

    constructs_ids = []
    cless = []
    output_sites = []
    input_sites = []
    mws = []
    if estimate_cells:
        from faust.utilities import estimate_cell_input
        estimated_cells_input = []
        estimated_cells_output = []
    for construct_id in tqdm(df['Target Gene'].unique(),
                             desc='Looping through target genes'):
        experiment = df[df['Target Gene'] == construct_id]
        control = df[df['Target Gene'].isin(controls)]
        for odds_ratio in [
                x for x in experiment.columns if x.endswith('_odds_ratio')
        ]:
            a = experiment[odds_ratio].values
            b = control[odds_ratio].values
            if input_type == 'matched':
                a = [x for x in a if not np.isinf(x) and not np.isnan(x)]
                b = [x for x in b if not np.isinf(x) and not np.isnan(x)]
            if len(a) == 0 or len(b) == 0:
                cless.append(np.nan)
                mws.append(np.nan)

            else:
                mw = mannwhitneyu(a, b, alternative='two-sided')
                cless.append(mw.statistic / len(a) / len(b))
                mws.append(mw.pvalue)
            output_site = odds_ratio.split('_odds_ratio')[0]
            output_sites.append(output_site)
            input_site = output2input_dict[odds_ratio.split('_odds_ratio')[0]]
            input_sites.append(input_site)
            constructs_ids.append(construct_id)
            if estimate_cells:
                estimated_cells_input.append(
                    estimate_cell_input(df, input_site, construct_id,
                                        count_threshold))
                estimated_cells_output.append(
                    estimate_cell_input(df, output_site, construct_id,
                                        count_threshold))

    summary_df = pd.DataFrame(constructs_ids, columns=['gene'])
    summary_df['output_site'] = output_sites
    summary_df['input_site'] = input_sites
    summary_df['CommonLanguageEffectSize'] = cless
    summary_df['MannWhitneyP'] = mws
    summary_df['BH_q'] = nan_fdrcorrection_q(summary_df['MannWhitneyP'].values)
    if estimate_cells:
        summary_df['input_estimated_cell_count'] = estimated_cells_input
        summary_df['output_estimated_cell_count'] = estimated_cells_output
    return summary_df


def nan_fdrcorrection_q(pvalues):
    from statsmodels.stats.multitest import fdrcorrection
    pvalues = np.array(pvalues)
    realmask = ~np.isnan(pvalues)
    qvalues = np.array([np.nan] * len(pvalues))
    qvalues[realmask] = fdrcorrection(pvalues[realmask])[1]
    return qvalues


def get_replicate_aggregated_statistics(summary_df,
                                        aggregation_column=None,
                                        inplace=False):
    """

    Parameters
    ----------
    summary_df :
        
    aggregation_column :
         (Default value = None)
    inplace :
         (Default value = False)

    Returns
    -------

    """
    if aggregation_column not in summary_df.columns:
        raise Exception("aggregation_column must be column in summary_df")

    from statsmodels.stats.multitest import fdrcorrection
    from scipy.stats import combine_pvalues
    from faust.utilities import nan_fdrcorrection_q

    if not inplace:
        summary_df = summary_df.copy()
    for aggregation in summary_df[aggregation_column].unique():

        aggregate_gene_max_p = summary_df[
            summary_df[aggregation_column] == aggregation].groupby(
                'gene')['MannWhitneyP'].apply(max)
        gene2aggregate_union_q = {
            gene: q
            for gene, q in zip(
                aggregate_gene_max_p.index.values,
                nan_fdrcorrection_q(aggregate_gene_max_p.values))
        }
        summary_df['{}_union_BH_q'.format(aggregation)] = summary_df[
            'gene'].apply(lambda x: gene2aggregate_union_q[x])

        aggregate_gene_fisher_combined_p = summary_df[
            summary_df[aggregation_column] == aggregation].groupby('gene')[
                'MannWhitneyP'].apply(combine_pvalues).apply(lambda x: x[1])
        gene2aggregate_intersection_q = {
            gene: q
            for gene, q in zip(
                aggregate_gene_fisher_combined_p.index.values,
                nan_fdrcorrection_q(aggregate_gene_fisher_combined_p.values))
        }
        summary_df['{}_intersection_BH_q'.format(aggregation)] = summary_df[
            'gene'].apply(lambda x: gene2aggregate_intersection_q[x])
    if not inplace:
        return summary_df


def estimate_read_error_singlets(n_observed_unique_grna_umis, n_observed_zeros,
                                 log2_counts_sum):
    from scipy.optimize import bisect

    left_side = lambda false_zeros: np.log(n_observed_unique_grna_umis /
                                           (n_observed_zeros - false_zeros))
    right_side = lambda false_zeros: log2_counts_sum / (
        n_observed_unique_grna_umis - false_zeros)
    f = lambda false_zeros: left_side(false_zeros) - right_side(false_zeros)
    estimate = bisect(f, 0, n_observed_zeros)
    return estimate / n_observed_zeros


def predicted_input(n_possible_grna_umi, n_detected_grna_umi):
    return n_possible_grna_umi * np.log(
        n_possible_grna_umi / (n_possible_grna_umi - n_detected_grna_umi))


def estimate_cell_input(df,
                        sample,
                        target_gene,
                        count_threshold,
                        target_gene_col='Target Gene'):
    if type(target_gene) is str:
        target_gene = [target_gene]
    from faust.utilities import predicted_input
    n_detected_grna_umi = (df[(df[target_gene_col].isin(target_gene)) & \
                             (df[sample]>count_threshold)][sample]>0).sum()
    n_possible_grna_umi = df[df[target_gene_col].isin(target_gene)].shape[0]
    if n_possible_grna_umi == n_detected_grna_umi:
        return np.nan
    else:
        return predicted_input(n_possible_grna_umi, n_detected_grna_umi)


def get_mageck_compatible_df(df,
                             sgRNA_col='Construct Barcode',
                             gene_col='Construct IDs',
                             UMI_col='UMI',
                             append_UMI=True,
                             output=None):
    df = df.copy()
    if append_UMI:
        df['sgRNA'] = df[sgRNA_col] + '_' + df[UMI_col]
    else:
        df['sgRNA'] = df[sgRNA_col]
    df['gene'] = df[gene_col]
    df = df[['sgRNA', 'gene'] + list(df.columns[df.dtypes == int])]
    if output is not None:
        if output.endswith('.csv'):
            df.to_csv(output, index=False)
        else:
            df.to_csv(output, index=False, sep='\t')
    return df


def get_mageck_ibar_compatible_df(df,
                                  sgRNA_col='Construct Barcode',
                                  gene_col='Construct IDs',
                                  UMI_col='UMI',
                                  output=None):
    df = df.copy()

    df['guide'] = df[sgRNA_col]
    df['gene'] = df[gene_col]
    df['barcode'] = df[UMI_col]
    df = df[['gene', 'guide', 'barcode'] + list(df.columns[df.dtypes == int])]
    if output is not None:
        if not output.endswith('.csv'):
            output = output + '.csv'
        df.to_csv(output, index=False)
    return df


def get_zfc_compatible_df(df,
                          sgRNA_col='Construct Barcode',
                          gene_col='Construct IDs',
                          UMI_col='UMI',
                          ctrl_col=None,
                          exp_col=None,
                          output=None):
    df = df.copy()

    df['guide'] = df[sgRNA_col]
    df['gene'] = df[gene_col]
    df['barcode'] = df[UMI_col]
    df['ctrl'] = df[ctrl_col]
    df['exp'] = df[exp_col]
    df = df[['gene', 'guide', 'barcode', 'ctrl', 'exp']]
    if output is not None:
        df.to_csv(output, index=False, sep='\t')
    return df


def get_riger_compatible_df(df,
                            sgRNA_col='Construct Barcode',
                            UMI_col='UMI',
                            append_UMI=True,
                            gene_col='Construct IDs',
                            score_col=None,
                            rank_col=None):
    df = df.copy()
    if append_UMI:
        df['Construct'] = df[sgRNA_col] + '_' + df[UMI_col]
    else:
        df['Construct'] = df[sgRNA_col]
    df['GeneSymbol'] = df[gene_col]
    df['NormalizedScore'] = df[score_col]
    df = df.sort_values(rank_col)
    df['Construct Rank'] = df.reset_index().index.values
    df['HairpinWeight'] = 1.0
    # Construct
    #    The name of the hairpin.
    # GeneSymbol
    #    A unique name for the gene.
    # NormalizedScore
    #   The hairpin score.
    # Construct Rank
    #    The hairpin rank.
    # HairpinWeight
    #  0 to 1--0 is unweighted
    return df[[
        'Construct', 'GeneSymbol', 'NormalizedScore', 'Construct Rank',
        'HairpinWeight'
    ]]


def run_alternative_test(df,
                         test=None,
                         exp_col=None,
                         ctrl_col=None,
                         sgRNA_col='Construct Barcode',
                         gene_col='Construct IDs',
                         UMI_col='UMI',
                         output=''):
    from faust.utilities import get_mageck_compatible_df, get_mageck_ibar_compatible_df, get_zfc_compatible_df
    implemented_tests = ['mageck', 'mageck-ibar', 'zfc']
    if test not in implemented_tests:  #+['riger']:
        raise Exception("test must be one of {}".format(implemented_tests))
    if output == '' or type(output) != str:
        raise Exception("output must be a non-null string")
    if test == 'mageck':
        transformed_df = get_mageck_compatible_df(df,
                                                  output=output,
                                                  gene_col=gene_col,
                                                  sgRNA_col=sgRNA_col,
                                                  UMI_col=UMI_col)
        command = 'mageck test -k {0} -t "{1}" -c "{2}" -n {0}'.format(
            output, exp_col, ctrl_col)
        primary_output = "{}.gene_summary.txt".format(output)
        secondary_output = None
    elif test == 'mageck-ibar':
        transformed_df = get_mageck_ibar_compatible_df(df,
                                                       output=output,
                                                       gene_col=gene_col,
                                                       sgRNA_col=sgRNA_col,
                                                       UMI_col=UMI_col)
        maxreps = transformed_df['gene'].value_counts().max() + 1
        if not output.endswith('.csv'):
            output = output + '.csv'
        command = 'mageck-ibar --RRApath "RRA --max-sgrnapergene-permutation {3}" -i {0} -t "{1}" -c "{2}" -o {0}'.format(
            output, exp_col, ctrl_col, maxreps)
        primary_output = output + ".gene.high.txt"
        secondary_output = output + ".gene.low.txt"
    elif test == 'zfc':
        transformed_df = get_zfc_compatible_df(df,
                                               exp_col=exp_col,
                                               ctrl_col=ctrl_col,
                                               output=output,
                                               gene_col=gene_col,
                                               sgRNA_col=sgRNA_col,
                                               UMI_col=UMI_col)
        command = "zfc --input {0} -o {1}".format(output, output + '_zfc')
        primary_output = "{}_zfc_gene.txt".format(output)
        secondary_output = None
    os.system(command)
    if secondary_output is None:
        return pd.read_table(primary_output)
    else:
        return pd.read_table(primary_output), pd.read_table(secondary_output)


def count_gpp_output(sgRNA_input, barcode_input, valid_constructs, valid_umis,
                     conditions, output):
    from tqdm import tqdm
    import pyfastx

    constructs = np.genfromtxt(valid_constructs, dtype=str)
    umis = np.genfromtxt(valid_umis, dtype=str)
    conditions = pd.read_csv(conditions, header=None).set_index(0)[1].to_dict()

    #    construct2counts = {construct:{umi:0 for umi in umis} for construct in constructs}
    construct2counts = {
        condition:
        {construct: {umi: 0
                     for umi in umis}
         for construct in constructs}
        for condition in conditions.keys()
        if type(conditions[condition]) == str
    }

    sgrna = pyfastx.Fastx(sgRNA_input)
    index = pyfastx.Fastx(barcode_input)
    for (name_sgrna, seq_sgrna, qual_sgrna,
         comment_sgrna), (name_barcode, seq_barcode, qual_barcode,
                          comment_barcode) in zip(
                              tqdm(sgrna, desc='looping through fastq'),
                              index):
        try:
            construct, umi = seq_sgrna.split(prefix)
            umi = umi[0:6]
            construct2counts[construct][umi] += 1
        except:
            pass

# df = pd.DataFrame.from_dict(construct2counts, 'index')#.T.reset_index()
# df = pd.DataFrame(df.stack()).reset_index()
# df.columns = ['Construct','UMI','count']
    df = pd.concat([
        pd.DataFrame.from_dict(construct2counts[key]).stack().rename(
            conditions[key]) for key in construct2counts.keys()
    ],
                   axis=1)
    df.to_csv(output, index=False)
