from collections import defaultdict 
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy import stats

import click


BED12_COLUMNS = [
    'chrom', 'start', 'end', 'srna_id',
    'allen_score', 'strand', 'cut_site',
    'cut_site_end', 'color',
    'nblocks', 'block_lengths', 'block_starts',
    'log_ratio'
]


OUTPUT_COLUMNS = [
    'chrom', 'start', 'end', 'srna_id',
    'score', 'strand', 'cut_site',
    'cut_site_end', 'color',
    'nblocks', 'block_lengths', 'block_starts',
    'log_odds', 'p_val', 'fdr',
    'log_ratio', 'median_null_log_ratio',
    'allen_score',
    'log_ratio_p_val', 'allen_score_p_val',
]

TEST_COLUMNS = ['srna_id', 'allen_score', 'log_ratio']


def read_bed_file(bed_fn):
    bed_records = pd.read_csv(
        bed_fn,
        sep='\t',
        names=BED12_COLUMNS,
        dtype={'block_lengths': str, 'block_starts': str}
    )
    return bed_records


def make_null_distributions(gstar_shuffled):
    distributions = {}
    for srna_id, group in gstar_shuffled.groupby('srna_id'):
        distributions[srna_id] = (
            group.log_ratio.values, group.allen_score.values
        )
    return distributions
    

def test_significance(gstar_srnas, gstar_shuffled):
    null_distributions = make_null_distributions(gstar_shuffled)
    median_null_dists = []
    log_ratio_p_vals = []
    allen_score_p_vals = []
    combined_p_vals = []
    for _, srna_id, allen_score, log_ratio in gstar_srnas[TEST_COLUMNS].itertuples():
        log_ratio_dist, allen_score_dist = null_distributions[srna_id]
        lr_p = ((log_ratio_dist > log_ratio).sum() + 1) / (len(log_ratio_dist) + 1)
        as_p = ((allen_score_dist < allen_score).sum() + 1) / (len(allen_score_dist) + 1)
        _, p_val = stats.combine_pvalues([lr_p, as_p], method='fisher')
        median_null_dists.append(np.median(log_ratio_dist))
        log_ratio_p_vals.append(lr_p)
        allen_score_p_vals.append(as_p)
        combined_p_vals.append(p_val)

    gstar_srnas['median_null_log_ratio'] = median_null_dists
    gstar_srnas['log_odds'] = (
        gstar_srnas.log_ratio - gstar_srnas.median_null_log_ratio
    )

    gstar_srnas['log_ratio_p_val'] = log_ratio_p_vals
    gstar_srnas['allen_score_p_val'] = allen_score_p_vals
    gstar_srnas['p_val'] = combined_p_vals
    _, gstar_srnas['fdr'], *_ = multipletests(
        gstar_srnas.p_val, method='fdr_bh'
    )
    gstar_srnas['score'] = np.round(- np.log10(gstar_srnas.fdr))
    return gstar_srnas


@click.command()
@click.option('-s', '--GSTAr-sRNAs', required=True)
@click.option('-c', '--GSTAr-shuffled', required=True)
@click.option('-o', '--output-fn', required=True)
def nanopore_endcut(gstar_srnas, gstar_shuffled, output_fn):
    gstar_srnas = read_bed_file(gstar_srnas)
    gstar_shuffled = read_bed_file(gstar_shuffled)
    res = test_significance(gstar_srnas, gstar_shuffled)
    res[OUTPUT_COLUMNS].to_csv(output_fn, sep='\t', index=False, header=False)


if __name__ == '__main__':
    nanopore_endcut()