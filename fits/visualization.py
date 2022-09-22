def plot_top_hits(summary_df, ascending=False,nhits=10,figsize=(9,5),fig=None, ax=None, 
                  sort_criterion='CommonLanguageEffectSize',hue=None):
    """

    Parameters
    ----------
    summary_df :
        
    ascending :
         (Default value = False)
    nhits :
         (Default value = 10)
    figsize :
         (Default value = (9)
    5) :
        
    fig :
         (Default value = None)
    ax :
         (Default value = None)
    sort_criterion :
         (Default value = 'CommonLanguageEffectSize')
    hue :
         (Default value = None)

    Returns
    -------

    """
    import matplotlib.pyplot as plt
    from panopticon.visualization import swarmviolin
    
    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    elif (fig is None and ax is not None) or (fig is not None and ax is None):
        raise Exception("either both or neither of fig, ax must be None")
    summary_df = summary_df.copy()
    hits = summary_df.groupby('gene')[sort_criterion].mean().sort_values(ascending=False).head(nhits).index.values

    swarmviolin(summary_df[summary_df['gene'].isin(hits)].sort_values(sort_criterion), y='CommonLanguageEffectSize',
                x='gene',ax=ax,hue=hue,
             #   annotate_hue_pvalue_fmt_str='p: {0:.2g}', pvalue='ttest',split=True,
             #   custom_annotation_dict=gene2tumor_intersection_q, violinplot_kwargs={'color':'b'},
                # swarmplot_kwargs={'color':'k'}

           )    
    return fig, ax