import pandas as pd
import os


def read_gpp_output(
        folders,
        chipfile=None,
        barcode2gene_dict=None,
        indices=['Construct Barcode', 'Construct IDs', 'UMI', 'Target Gene'],
        dropcols=[]):
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


def get_summary_df(df, controls, inputs, outputs):
    from tqdm import tqdm
    import pandas as pd
    df = df.copy()

    df['Input Sum'] = df[inputs].sum(axis=1)
    df['Input Mean'] = df[inputs].mean(axis=1)
    df = df[df['Input Sum'] > 0]
    odds_ratio_dict = {}
    for output in outputs:
        odds_ratio_dict[output] = []
        #for i, row in df.iterrows():
        #    odds_ratio_dict[output].append(row[output]/row['Input Sum'])
        odds_ratio_dict[output] = df[output] / df['Input Mean']
    for key in odds_ratio_dict.keys():
        df[key + '_odds_ratio'] = odds_ratio_dict[key]

    from scipy.stats import mannwhitneyu
    constructs_ids = []
    cless = []
    #tissues = []
    mws = []
    for construct_id in tqdm(df['Target Gene'].unique(),
                             desc='Looping through target genes'):
        experiment = df[df['Target Gene'] == construct_id]
        control = df[df['Target Gene'].isin(controls)]
        for odds_ratio in [
                x for x in experiment.columns if x.endswith('_odds_ratio')
        ]:
            mw = mannwhitneyu(experiment[odds_ratio].values,
                              control[odds_ratio].values,
                              alternative='two-sided')
            cles = mw.statistic / len(experiment) / len(control)
            constructs_ids.append(construct_id)
            cless.append(cles)
            #  tissues.append(odds_ratio.split(' ')[0])
            mws.append(mw.pvalue)
    summary_df = pd.DataFrame(constructs_ids, columns=['gene'])
    #  summary_df['tissue'] = tissues
    summary_df['CommonLanguageEffectSize'] = cless
    summary_df['MannWhitneyP'] = mws
    from statsmodels.stats.multitest import fdrcorrection

    summary_df['BH_q'] = fdrcorrection(summary_df['MannWhitneyP'].values, )[1]
    return summary_df
