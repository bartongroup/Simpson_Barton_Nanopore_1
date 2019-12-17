import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests

import pysam
import click

from chimerID.io import parse_gtf_flat_exon_invs, parse_pysam_aln
from chimerID.intervals import intersect_spliced_invs


COLUMNS = ['chrom', 'start', 'end', 'gene_id',
           'strand', 'nreads_a', 'nreads_b',
           'median_a', 'ci_lower_a', 'ci_upper_a',
           'median_b', 'ci_lower_b', 'ci_upper_b',
           'ks', 'ks_p_val', 'mwu', 'mwu_p_val']

COL_ORDER = ['chrom', 'start', 'end', 'gene_id',
             'score', 'strand', 'nreads_a', 'nreads_b',
             'median_a', 'ci_lower_a', 'ci_upper_a',
             'median_b', 'ci_lower_b', 'ci_upper_b',
             'ks', 'ks_p_val', 'ks_fdr',
             'mwu', 'mwu_p_val', 'mwu_fdr']


def get_polya_dist(chrom, start, end, strand, bam, gene_invs=None,
                   overlap_thresh=200, gene_frac_thresh=0.2,
                   read_frac_thresh=0.5):
    if gene_invs is None:
        gene_invs = [[start, end]]
    gene_ln = sum([e - s for s, e in gene_invs])
    polya_lengths = []
    for aln in bam.fetch(chrom, start, end):
        _, _, _, read_id, read_strand, read_invs, _, _  = parse_pysam_aln(aln)
        if strand != read_strand:
            continue
        read_ln = sum([e - s for s, e in read_invs])
        abs_overlap = intersect_spliced_invs(gene_invs, read_invs)
        read_frac = abs_overlap / read_ln
        gene_frac = abs_overlap / gene_ln
        if abs_overlap >= overlap_thresh and \
               read_frac >= read_frac_thresh and \
               gene_frac >= gene_frac_thresh:
            pa = aln.get_tag('pA')
            polya_lengths.append(pa)
    return np.array(polya_lengths)


def test_polya_distrib_changes(bam_a_fn, bam_b_fn, gtf_fn, nreads=10):
    res = []
    with pysam.AlignmentFile(bam_a_fn) as bam_a, pysam.AlignmentFile(bam_b_fn) as bam_b:
        for gene_id, chrom, strand, gene_invs in parse_gtf_flat_exon_invs(gtf_fn):
            start, end = gene_invs[0][0], gene_invs[-1][1]
            polya_a = get_polya_dist(chrom, start, end, strand, bam_a, gene_invs)
            polya_b = get_polya_dist(chrom, start, end, strand, bam_b, gene_invs)
            if len(polya_a) > nreads and len(polya_b) > nreads:
                ks, ks_p_val = ks_2samp(polya_a, polya_b)
                mwu, mwu_p_val = mannwhitneyu(polya_a, polya_b)
                res.append(
                    [chrom, start, end, gene_id, strand,
                     len(polya_a), len(polya_b),
                     np.median(polya_a), np.percentile(polya_a, 2.5), np.percentile(polya_a, 97.5),
                     np.median(polya_b), np.percentile(polya_b, 2.5), np.percentile(polya_b, 97.5),
                     ks, ks_p_val, mwu, mwu_p_val]
                )
    res = pd.DataFrame(res, columns=COLUMNS)
    _, ks_fdr, *_ = multipletests(res.ks_p_val, method='fdr_bh')
    res['ks_fdr'] = ks_fdr
    res['score'] = np.negative(np.log10(res.ks_p_val))
    _, mwu_fdr, *_ = multipletests(res.mwu_p_val, method='fdr_bh')
    res['mwu_fdr'] = mwu_fdr
    res = res[COL_ORDER]
    return res


def get_median_polya_lengths(gtf_fn, bam_fns, nreads=10):
    res = []
    bams = [pysam.AlignmentFile(b) for b in bam_fns]
    for gene_id, chrom, strand, gene_invs in parse_gtf_flat_exon_invs(gtf_fn):
        start, end = gene_invs[0][0], gene_invs[-1][1]
        index = [chrom, start, end, strand, gene_id]
        row = []
        exprs = []
        for bam in bams:
            polya_dist = get_polya_dist(chrom, start, end, strand, bam, gene_invs)
            row.append(np.median(polya_dist))
            exprs.append(len(polya_dist))
        exprs = np.array(exprs)
        if (exprs > nreads).all():
            res.append(index + row)
    res = pd.DataFrame(
        res, 
        columns=['chrom', 'start', 'end', 'strand', 'gene_id'] + list(bam_fns),
    )
    return res


@click.group()
def cli():
    pass


@cli.command()
@click.option('-a', '--bam-a-fn', required=True)
@click.option('-b', '--bam-b-fn', required=True)
@click.option('-g', '--gtf-fn', required=True)
@click.option('-o', '--output-bed', required=True)
@click.option('-n', '--nreads', required=False, default=10)
def test_differential_polya(bam_a_fn, bam_b_fn, gtf_fn, output_bed, nreads):
    results = test_polya_distrib_changes(bam_a_fn, bam_b_fn, gtf_fn, nreads)
    results.to_csv(
        output_bed, sep='\t',
        na_rep='.', header=False,
        index=False,
        float_format='%.10g'
    )


@cli.command()
@click.option('-g', '--gtf-fn', required=True)
@click.option('-o', '--output-tsv', required=True)
@click.argument('bams', nargs=-1)
def median_polya_lengths(gtf_fn, output_tsv, bams):
    results = get_median_polya_lengths(gtf_fn, bams)
    results.to_csv(
        output_tsv, sep='\t',
        na_rep='.', header=True,
        index=False,
        float_format='%.10g'
    )


if __name__ == '__main__':
    cli()