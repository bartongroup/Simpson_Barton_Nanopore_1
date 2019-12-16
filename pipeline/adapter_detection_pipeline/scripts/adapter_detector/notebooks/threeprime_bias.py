import re
import itertools as it
import numpy as np

import pysam

def flatten(bundle):
    flattened = []
    all_invs = iter(sorted(it.chain(*bundle)))
    inv_start, inv_end = next(all_invs)
    for start, end in all_invs:
        if start <= inv_end:
            inv_end = max(inv_end, end)
        else:
            flattened.append([inv_start, inv_end])
            inv_start, inv_end = start, end
    if not flattened or flattened[-1] != [inv_start, inv_end]:
        flattened.append([inv_start, inv_end])
    return flattened


def get_gtf_gene_id(attrs):
    return re.search('gene_id \"(.*?)\";', attrs).group(1)


def get_gtf_exons(gtf_fn):
    with open(gtf_fn) as gtf:
        for record in gtf:
            record = record.strip().split('\t')
            if record[2] == 'exon':
                gene_id = get_gtf_gene_id(record[8])
                yield record[0], int(record[3]) - 1, int(record[4]), gene_id, record[6]


def parse_gtf_flat_exon_invs(gtf_fn):
    gene_cluster = []
    gtf_iter = get_gtf_exons(gtf_fn)
    curr_chrom, start, end, curr_gene_id, curr_strand = next(gtf_iter)
    gene_cluster.append([[start, end]])
    for chrom, start, end, gene_id, strand in gtf_iter:
        if gene_id != curr_gene_id:
            yield curr_gene_id, curr_chrom, curr_strand, flatten(gene_cluster)
            curr_gene_id, curr_chrom, curr_strand = gene_id, chrom, strand
            gene_cluster = []
            gene_cluster.append([[start, end]])
        else:
            gene_cluster.append([[start, end]])
    if gene_cluster:
        yield curr_gene_id, curr_chrom, curr_strand, flatten(gene_cluster)


def get_coverage(gtf_fn, bam_fn):
    gene_coverage = {}
    with pysam.AlignmentFile(bam_fn) as bam:
        for gene_id, chrom, strand, invs in parse_gtf_flat_exon_invs(gtf_fn):
            ln = sum([e - s for s, e in invs])
            cov = np.zeros(ln, dtype=np.uint)
            offset = 0
            for start, end in invs:
                for aln in bam.fetch(chrom, start, end):
                    if ['+', '-'][aln.is_reverse] == strand:
                        cov_start = max(aln.reference_start, start) - start
                        cov_end = min(aln.reference_end, end) - start
                        cov[cov_start + offset: cov_end + offset] += 1
                offset += end - start
            gene_coverage[gene_id] = cov
    return gene_coverage


def three_prime_bias_metric(arr):
    arr = np.trim_zeros(arr)
    med = np.median(arr)
    mad = np.abs(arr - med)
    return np.median(mad) / np.max(mad)