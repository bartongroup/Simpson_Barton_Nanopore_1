import pandas as pd

from itertools import combinations
from operator import itemgetter

from .intervals import intersect_spliced_invs
from .io import parse_pysam_aln, bed12_parse_exons


def find_overlapping_genes(read_invs, genes_tabix, blacklist_genes,
                           chrom, start, end, strand,
                           ovrlp_thresh=0.2,
                           abs_ovrlp_thresh=200):
    overlapping_gene_invs = []
    overlapping_gene_ids = []
    for gene_record in genes_tabix.fetch(chrom, start, end):
        gene_record = gene_record.split()
        gene_id = gene_record[3]
        gene_strand = gene_record[5]
        if gene_strand != strand or gene_id in blacklist_genes:
            continue
        gene_start, gene_end, gene_invs = bed12_parse_exons(gene_record)
        gene_length = sum([e - s for s, e in gene_invs])
        abs_overlap = intersect_spliced_invs(read_invs, gene_invs)
        overlap = abs_overlap / gene_length
        if overlap > ovrlp_thresh and abs_overlap > abs_ovrlp_thresh:
            overlapping_gene_invs.append([gene_start, gene_end])
            overlapping_gene_ids.append(gene_id)
    return overlapping_gene_invs, overlapping_gene_ids


def read_wholly_contained(r_start, r_end, gene_invs, gene_ids, tol):
    for (g_start, g_end), g_id in zip(gene_invs, gene_ids):
        if r_start > (g_start - tol) and r_end < (g_end + tol):
            return True, g_id
    else:
        return False, None


def identify_gene_of_origin(read_strand, gene_invs, gene_ids):
    if read_strand == '+':
        order = sorted(range(len(gene_ids)),
                       key=lambda i: gene_invs[i][0])
    elif read_strand == '-':
        order = sorted(range(len(gene_ids)),
                       key=lambda i: gene_invs[i][1],
                       reverse=True)
    sorted_gene_ids = [gene_ids[i] for i in order]
    gene_of_origin = sorted_gene_ids[0]
    downstream_genes = '|'.join(sorted_gene_ids[1:])
    return gene_of_origin, downstream_genes


def identify_chimeras(bam_iter, tabix, blacklist_genes,
                      ovrlp_thresh=0.2,
                      abs_ovrlp_thresh=200):
    for n, aln in enumerate(bam_iter):
        chrom, start, end, read_id, strand, invs, is_secondary, mapq = aln
        if mapq == 0 or is_secondary:
            # low quality mappings on repetitive regions
            # are likely to generate spurious chimeras
            continue
        record = (invs, chrom, strand)
        overlap_invs, gene_ids = find_overlapping_genes(
            invs, tabix, blacklist_genes,
            chrom, start, end, strand, 
            ovrlp_thresh, abs_ovrlp_thresh
        )
        if not gene_ids:
            continue
        rwc, g = read_wholly_contained(
            start, end,
            overlap_invs, gene_ids,
            tol=abs_ovrlp_thresh
        )
        if len(gene_ids) > 1 and not rwc:
            gene_of_origin, downstream_genes = identify_gene_of_origin(
                strand, overlap_invs, gene_ids)
            yield True, (gene_of_origin, downstream_genes, strand), read_id, record, n
        elif rwc:
            yield False, (g, strand), read_id, record, n
        elif gene_ids:
            yield False, (gene_ids[0], strand), read_id, record, n