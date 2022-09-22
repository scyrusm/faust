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
                   verbose=True):
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
    if input_type not in ['single', 'matched']:
        raise Exception("input_type must be one of 'single', 'matched'")
    df = df.copy()
    if input_type == 'single':
        df['Input Sum'] = df[inputs].sum(axis=1)
        df['Input Mean'] = df[inputs].mean(axis=1)
        inputs = ['Input Sum'] * len(outputs)
        df = df[df['Input Mean'] > 0]
    else:
        if len(outputs) != len(inputs):
            raise Exception(
                "outputs and inputs must be equal length when input_type is 'matched'"
            )
    for output_col, input_col in zip(outputs, inputs):
        if verbose:
            print("Comparing input", input_col, "with output", output_col)
        df[output_col + '_odds_ratio'] = np.divide(df[output_col],
                                                   df[input_col])

    from scipy.stats import mannwhitneyu
    constructs_ids = []
    cless = []
    #tissues = []
    output_sites = []
    mws = []
    for construct_id in tqdm(df['Target Gene'].unique(),
                             desc='Looping through target genes'):
        experiment = df[df['Target Gene'] == construct_id]
        control = df[df['Target Gene'].isin(controls)]
        for odds_ratio in [
                x for x in experiment.columns if x.endswith('_odds_ratio')
        ]:
            a = experiment[odds_ratio].values
            b = control[odds_ratio].values
            if input_type == 'matched': # This must be done before calculating mannwhitneyu
                a = [x for x in a if not np.isinf(x) and not np.isnan(x)]
                b = [x for x in b if not np.isinf(x) and not np.isnan(x)]

            mw = mannwhitneyu(a, b, alternative='two-sided')
            cles = mw.statistic / len(experiment) / len(control)
            constructs_ids.append(construct_id)
            cless.append(cles)
            #  tissues.append(odds_ratio.split(' ')[0])
            output_sites.append(odds_ratio.split('_')[0])
            mws.append(mw.pvalue)
    summary_df = pd.DataFrame(constructs_ids, columns=['gene'])
    #  summary_df['tissue'] = tissues
    summary_df['output_site'] = output_sites
    summary_df['CommonLanguageEffectSize'] = cless
    summary_df['MannWhitneyP'] = mws
    from statsmodels.stats.multitest import fdrcorrection

    summary_df['BH_q'] = fdrcorrection(summary_df['MannWhitneyP'].values, )[1]
    return summary_df


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
    if not inplace:
        summary_df = summary_df.copy()
    for aggregation in summary_df[aggregation_column].unique():

        aggregate_gene_max_p = summary_df[
            summary_df[aggregation_column] == aggregation].groupby(
                'gene')['MannWhitneyP'].apply(max)
        gene2aggregate_union_q = {
            gene: q
            for gene, q in zip(aggregate_gene_max_p.index.values,
                               fdrcorrection(aggregate_gene_max_p.values)[1])
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
                fdrcorrection(aggregate_gene_fisher_combined_p.values)[1])
        }
        summary_df['{}_intersection_BH_q'.format(aggregation)] = summary_df[
            'gene'].apply(lambda x: gene2aggregate_intersection_q[x])
    if not inplace:
        return summary_df
