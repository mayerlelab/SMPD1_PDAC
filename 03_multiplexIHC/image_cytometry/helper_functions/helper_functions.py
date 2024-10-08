import pandas as pd
import re
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns



def get_cluster_proportions(adata,
                            cluster_key="cluster_final",
                            sample_key="replicate",
                            drop_values=None):
    """
    Input
    =====
    adata : AnnData object
    cluster_key : key of `adata.obs` storing cluster info
    sample_key : key of `adata.obs` storing sample/replicate info
    drop_values : list/iterable of possible values of `sample_key` that you don't want
    
    Returns
    =======
    pd.DataFrame with samples as the index and clusters as the columns and 0-100 floats
    as values
    """
    
    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([cluster_key, sample_key]).size().reset_index()
    props = sizes.set_index(cluster_key).groupby(sample_key).apply(lambda x: 100 * x / x.sum()).reset_index()
    props = props.pivot(columns=sample_key, index=cluster_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)
    
    if drop_values is not None:
        for drop_value in drop_values:
            props.drop(drop_value, axis=0, inplace=True)
    return props


def plot_cluster_proportions(cluster_props, 
                             cluster_palette=None,
                             xlabel_rotation=0, dpi=100): 
    fig, ax = plt.subplots(dpi=dpi)
    fig.patch.set_facecolor("white")
    
    cmap = None
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, 
            n_colors=len(cluster_palette), 
            as_cmap=True)
   
    cluster_props.plot(
        kind="bar", 
        stacked=True, 
        ax=ax, 
        legend=None,
        edgecolor = "black",
        colormap=cmap
    )
    
    ax.legend(bbox_to_anchor=(1, 1), frameon=False)
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=xlabel_rotation)
    ax.set_xlabel(cluster_props.index.name.capitalize())
    ax.set_ylabel("Proportion of cells")
    fig.tight_layout()
    ax.grid(False)
    
    return None

def multiple_testing_correction(ps, alpha=0.05, method='benjamini-hochberg', **kwargs):
    """ correct pvalues for multiple testing and add corrected `q` value
    :param ps: list of pvalues
    :param alpha: significance level default : 0.05
    :param method: multiple testing correction method [bonferroni|benjamini-hochberg]
    :returns (q, rej): two lists of q-values and rejected nodes
    """
    # Implement copy from GOATools.
    _p = np.array(ps)
    q = _p.copy()
    rej = _p.copy()
    mask = ~np.isnan(_p)
    p = _p[mask]
    if method == 'bonferroni':
        q[mask] = p * len(p)
        rej[mask] = q[mask] < alpha
    elif method == 'benjamini-hochberg':
        _rej, _q = fdrcorrection(p, alpha)
        rej[mask] = _rej
        q[mask] = _q
    else:
        raise ValueError(method)
    return q, rej


def fdrcorrection(pvals, alpha=0.05):
    """ benjamini hocheberg fdr correction. inspired by statsmodels 
    """
    # Implement copy from GOATools.
    pvals = np.asarray(pvals)
    pvals_sortind = np.argsort(pvals)
    pvals_sorted = np.take(pvals, pvals_sortind)

    ecdffactor = _ecdf(pvals_sorted)
    reject = pvals_sorted <= ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True
    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1
    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[pvals_sortind] = pvals_corrected
    del pvals_corrected
    reject_ = np.empty_like(reject)
    reject_[pvals_sortind] = reject
    return reject_, pvals_corrected_

def _ecdf(x):
    nobs = len(x)
    return np.arange(1,nobs+1)/float(nobs)