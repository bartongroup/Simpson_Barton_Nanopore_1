import itertools as it
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def read_all_chimeric_counts(h5_fns, sample_names, normalise=True):
    chimeric_counts = {}
    all_gene_non_chimeric_counts = {}
    for sample, h5_fn in zip(sample_names, h5_fns):
        chimeric_counts[sample] = pd.read_hdf(h5_fn, key='chimera_counts')
        all_gene_non_chimeric_counts[sample] = pd.read_hdf(h5_fn, key='non_chimeric_counts')
        norm_factors = pd.read_hdf(h5_fn, key='norm_factors')
        if normalise:
            chimeric_counts[sample] /= norm_factors
            all_gene_non_chimeric_counts[sample] /= norm_factors
    chimeric_counts = pd.concat(
        chimeric_counts, axis=1,
        sort=True, names=['sample', 'boot'])
    all_gene_non_chimeric_counts = pd.concat(
        all_gene_non_chimeric_counts, axis=1,
        sort=True, names=['sample', 'boot'])
    downstream_genes = {(chimera, strand): downstream for
                        chimera, downstream, strand
                        in chimeric_counts.index}
    chimeric_counts = chimeric_counts.groupby(level=(0, 2), axis=0).sum()
    non_chimeric_counts = all_gene_non_chimeric_counts.loc[chimeric_counts.index].copy()
    counts = pd.concat(
        {'chimeric': chimeric_counts,
         'nonchimeric': non_chimeric_counts},
        axis=1, sort=True, names=['readtype', 'sample', 'boot'],
    ).reorder_levels(['sample', 'readtype', 'boot'], axis=1).fillna(0)
    return counts, downstream_genes


def get_bootstrap_stats(bootstraps, cond_a, cond_b):
    bootstraps = bootstraps.copy() + 0.5
    cond_a_ratio = (
        bootstraps.loc[:, (cond_a, 'chimeric', pd.IndexSlice[:])].values /
        bootstraps.loc[:, (cond_a, 'nonchimeric', pd.IndexSlice[:])].values
    )
    cond_b_ratio = (
        bootstraps.loc[:, (cond_b, 'chimeric', pd.IndexSlice[:])].values /
        bootstraps.loc[:, (cond_b, 'nonchimeric', pd.IndexSlice[:])].values
    )
    ks_stat = []
    ks_p_val = []
    for i in range(len(bootstraps)):
        ks, p_val = stats.ks_2samp(cond_a_ratio[i], cond_b_ratio[i])
        ks_stat.append(ks)
        ks_p_val.append(p_val)
    ks_stat = np.array(ks_stat)
    ks_p_val = np.array(ks_p_val)
    n_boots = len(bootstraps.columns.unique(level=2))
    boot_lr = {}
    for n, (i, j) in enumerate(it.product(range(n_boots), repeat=2)):
        cond_a_data = bootstraps.loc[:, (cond_a, pd.IndexSlice[:], i)].copy()
        cond_a_data.columns = cond_a_data.columns.droplevel(0)
        cond_b_data = bootstraps.loc[:, (cond_b, pd.IndexSlice[:], j)].copy()
        cond_b_data.columns = cond_b_data.columns.droplevel(0)
        r = ((cond_a_data['chimeric'].values / cond_a_data['nonchimeric'].values) /
             (cond_b_data['chimeric'].values / cond_b_data['nonchimeric'].values))
        boot_lr[n] = np.log2(r).ravel()
    boot_lr = pd.DataFrame.from_dict(boot_lr)
    boot_lr.index = bootstraps.index
    boot_lr_res = boot_lr.quantile([0.5, 0.025, 0.975], axis=1).T
    boot_lr_res.columns = ['logodds_median', 'logodds_lower_ci95', 'logodds_upper_ci95']
    boot_lr_res['logodds_mean'] = boot_lr.mean(axis=1)
    boot_lr_res['ks_stat'] = ks_stat
    boot_lr_res['ks_p_val'] = ks_p_val
    _, boot_lr_res['ks_fdr'], *_ = multipletests(boot_lr_res.ks_p_val, method='bonferroni')
    return boot_lr_res


def generate_bootstrapped_logodds(h5_fns, cond_a_sample_name, cond_b_sample_name):
    counts, downstream_genes = read_all_chimeric_counts(
        h5_fns, [cond_a_sample_name, cond_b_sample_name], normalise=False)
    median_counts = counts.groupby(level=['sample', 'readtype'], axis=1).median()
    median_counts = counts.groupby(level=['sample', 'readtype'], axis=1).median()
    median_counts.columns = (median_counts.columns.get_level_values(0) + '_' +
                             median_counts.columns.get_level_values(1))
    logodds_ratios = get_bootstrap_stats(
        counts, cond_a_sample_name, cond_b_sample_name)
    logodds_ratios['downstream_genes'] = pd.Series(downstream_genes)
    logodds_ratios = logodds_ratios.join(median_counts)
    return logodds_ratios