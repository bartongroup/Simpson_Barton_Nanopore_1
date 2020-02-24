from collections import Counter
import pysam

from .io import (
    bed12_interval_iterator,
    parse_gtf_flat_exon_invs,
    parse_pysam_aln,
    to_bed12, write_to_tabix
)
from .intervals import (
    fraction_of_shorter, flatten,
    intersect_spliced_invs
)


def format_merged_genes(chrom, invs, gene_ids, strand):
    gene_id = ','.join(gene_ids)
    return to_bed12(chrom, flatten(invs), strand=strand, name=gene_id)


def merge_overlapping_genes(gtf_fn, overlap_threshold=0.75):
    inv_cluster = {'+': [], '-': []}
    gene_id_cluster = {'+': [], '-': []}
    curr_chrom = {}
    curr_start = {}
    curr_end = {}
    gene_iter = parse_gtf_flat_exon_invs(gtf_fn)
    gene_id, chrom, strand, invs = next(gene_iter)
    curr_chrom[strand] = chrom
    inv_cluster[strand].append(invs)
    gene_id_cluster[strand].append(gene_id)
    curr_start[strand], curr_end[strand] = invs[0][0], invs[-1][1]
    for gene_id, chrom, strand, invs in gene_iter:
        if not curr_chrom.get(strand, False):
            curr_chrom[strand] = chrom
            inv_cluster[strand] = [invs,]
            gene_id_cluster[strand] = [gene_id,]
            curr_start[strand], curr_end[strand] = invs[0][0], invs[-1][1]
        if curr_chrom[strand] != chrom:
            yield format_merged_genes(curr_chrom[strand], inv_cluster[strand], gene_id_cluster[strand], strand)
            curr_chrom[strand] = chrom
            inv_cluster[strand] = [invs,]
            gene_id_cluster[strand] = [gene_id,]
            curr_start[strand], curr_end[strand] = invs[0][0], invs[-1][1]
        else:
            start, end = invs[0][0], invs[-1][1]
            f = fraction_of_shorter((start, end), (curr_start[strand], curr_end[strand]))
            if f < overlap_threshold:
                yield format_merged_genes(
                    curr_chrom[strand],
                    inv_cluster[strand],
                    gene_id_cluster[strand],
                    strand
                )
                curr_chrom[strand] = chrom
                inv_cluster[strand] = [invs,]
                gene_id_cluster[strand] = [gene_id,]
                curr_start[strand], curr_end[strand] = invs[0][0], invs[-1][1]
            else:
                inv_cluster[strand].append(invs)
                gene_id_cluster[strand].append(gene_id)
                curr_start[strand] = min(start, curr_start[strand])
                curr_end[strand] = max(end, curr_end[strand])
    for strand in '+-':
        if inv_cluster[strand]:
            yield format_merged_genes(
                curr_chrom[strand],
                inv_cluster[strand],
                gene_id_cluster[strand],
                strand
            )


def create_flattened_reference_bed(gtf_fn, output_prefix, overlap_threshold=0.75):
    bed_iter = merge_overlapping_genes(gtf_fn, overlap_threshold)
    output_bed = write_to_tabix(bed_iter, output_prefix)
    return output_bed
