import numpy as np
import pandas as pd

import pyBigWig as pybw
import click


BED12_COLUMNS = [
    'chrom', 'start', 'end', 'srna_id',
    'allen_score', 'strand', 'cut_site',
    'cut_site_end', 'color',
    'nblocks', 'block_lengths', 'block_starts'
]

ITERCOLS = [
    'chrom', 'start', 'cut_site', 'strand', 'block_lengths', 'block_starts'
]


def read_bed_file(bed_fn):
    bed_records = pd.read_csv(
        bed_fn,
        sep='\t',
        names=BED12_COLUMNS,
        dtype={'block_lengths': str, 'block_starts': str}
    )
    return bed_records


def parse_intervals(start, split_point, inv_lengths, inv_starts, strand):
    inv_lengths = np.fromstring(inv_lengths, sep=',', dtype=np.int)
    inv_starts = np.fromstring(inv_starts, sep=',', dtype=np.int) + start
    inv_ends = inv_starts + inv_lengths
    invs = np.sort(np.concatenate([inv_starts, inv_ends]))
    idx = np.searchsorted(invs, split_point)
    invs = np.insert(invs, idx, [split_point, split_point])
    upslic, downslic = (invs[:idx + 1].reshape(-1, 2),
                        invs[idx + 1:].reshape(-1, 2))
    if strand == '-':
        upslic, downslic = downslic, upslic
    return upslic, downslic


def cleavage_lr(bw, chrom, upslic, downslic):
    '''
    Generate a log odds ratio similar to endCut but for nanopore data:
    instead of the ratio of cut site 5' tags to 
    '''
    upstream_cut = sum([
        np.nansum(bw.values(str(chrom), i, j)) for i, j in upslic
    ])
    downstream_cut = sum([
        np.nansum(bw.values(str(chrom), i, j)) for i, j in downslic
    ])
    return np.log2((downstream_cut + 1) / (upstream_cut + 1))


def add_cleavage_log_ratio_to_gstar(bed_fn, fwd_bigwig_fn, rev_bigwig_fn):
    fwd_bw = pybw.open(fwd_bigwig_fn)
    rev_bw = pybw.open(rev_bigwig_fn)
    log_ratios = []
    srna_alns = read_bed_file(bed_fn)
    for (_, chrom, start, cut_site, strand,
          block_lengths, block_starts) in srna_alns[ITERCOLS].itertuples():
        bw = fwd_bw if strand == '+' else rev_bw
        upslic, downslic = parse_intervals(
            start, cut_site, block_lengths, block_starts, strand
        )
        log_ratios.append(cleavage_lr(bw, chrom, upslic, downslic))
    srna_alns['log_ratio'] = log_ratios
    fwd_bw.close()
    rev_bw.close()
    return srna_alns


@click.command()
@click.option('-f', '--fwd-bigwig', required=True)
@click.option('-r', '--rev-bigwig', required=True)
@click.option('-s', '--genomic-mapped-srnas', required=True)
@click.option('-o', '--output-fn', required=True)
def srna_cleavage_logratio(fwd_bigwig, rev_bigwig, genomic_mapped_srnas, output_fn):
    genomic_mapped_srnas = add_cleavage_log_ratio_to_gstar(
        genomic_mapped_srnas, fwd_bigwig, rev_bigwig
    )
    genomic_mapped_srnas.to_csv(output_fn, sep='\t', index=False, header=False)


if __name__ == '__main__':
    srna_cleavage_logratio()