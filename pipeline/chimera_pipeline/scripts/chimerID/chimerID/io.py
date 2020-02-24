import re
import os
import gzip
import subprocess
from tempfile import mkstemp
import numpy as np

import pysam

from .intervals import flatten

def bed12_parse_exons(record):
    start = int(record[1])
    end = int(record[2])
    exstarts = np.fromstring(record[11], sep=',') + start
    exends = exstarts + np.fromstring(record[10], sep=',')
    exons = np.dstack([exstarts, exends])[0]
    return start, end, exons


def bed12_interval_iterator(bed_fn, query=None):
    if query is None:
        query_contig = None
        query_start = None
        query_end = None
    else:
        query_contig, query_start, query_end = query
    with pysam.TabixFile(bed_fn) as bed:
        for record in bed.fetch(query_contig, query_start, query_end):
            record = record.split()
            chrom = record[0]
            strand = record[5]
            name = record[3]
            start, end, exons = bed12_parse_exons(record)
            yield chrom, start, end, name, strand, exons


def bam_cigar_to_invs(aln):
    invs = []
    start = aln.reference_start
    end = aln.reference_end
    strand = '-' if aln.is_reverse else '+'
    left = start
    right = left
    aln_length = 0
    for op, ln in aln.cigar:
        if op in (1, 4, 5):
            # does not consume reference
            continue
        elif op in (0, 2, 7, 8):
            # consume reference but do not add to invs yet
            right += ln
        elif op == 3:
            invs.append([left, right])
            aln_length += right - left
            left = right + ln
            right = left
    if right > left:
        invs.append([left, right])
    assert invs[0][0] == start
    assert invs[-1][1] == end
    return invs, start, end, strand


def parse_pysam_aln(aln):
    chrom = aln.reference_name
    read_id = aln.query_name
    invs, start, end, strand = bam_cigar_to_invs(aln)
    is_secondary = aln.is_secondary
    mapq = aln.mapping_quality
    return chrom, start, end, read_id, strand, invs, is_secondary, mapq


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


BED12_RECORD = (
    '{chrom}\t{start:d}\t{end:d}\t'
    '{name}\t60\t{strand}\t'
    '{start:d}\t{end:d}\t'
    '66,146,244\t{nblocks:d}\t'
    '{blocksizes}\t{blockstarts}\n'
)


def to_bed12(chrom, invs, strand, name='locus'):
    start = invs[0][0]
    end = invs[-1][1]
    nblocks = len(invs)
    lengths = ','.join(['{:d}'.format(e - s) for s, e in invs])
    starts = ','.join(['{:d}'.format(s - start) for s, e in invs])
    return BED12_RECORD.format(
        chrom=chrom, start=start, end=end,
        name=name, strand=strand, nblocks=nblocks,
        blocksizes=lengths, blockstarts=starts
    )


def write_to_tabix(bed_iter, output_prefix):
    handle, fn = mkstemp(suffix='.bed')
    with open(fn, 'w') as bed:
        for record in bed_iter:
            bed.write(record)
    output_bed = f'{output_prefix}.reference_loci.bed.gz'
    subprocess.check_call(
        f'sort -k1,1 -k2,2n {fn} | bgzip > {output_bed}; tabix -p bed {output_bed}',
        shell=True
    )
    os.close(handle)
    return output_bed