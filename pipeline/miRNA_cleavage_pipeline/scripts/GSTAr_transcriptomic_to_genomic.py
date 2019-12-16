from collections import defaultdict 
import numpy as np
import pandas as pd

import click

GSTAR_TRANSCRIPTOMIC_COLUMNS = [
    'srna_id', 'transcript_id',
    'start', 'end', 'cut_site', 'allen_score'
]

GSTAR_GENOMIC_COLUMNS = [
    'chrom', 'start', 'end',
    'srna_id', 'transcript_id',
    'cut_site','allen_score',
    'strand', 'invs',
]


def parse_exons(record):
    start = int(record[1])
    exstarts = np.fromstring(record[11], sep=',', dtype=np.int) + start
    exlengths = np.fromstring(record[10], sep=',', dtype=np.int)
    return exstarts, exlengths


def read_annotation(bed_fn):
    annotation = {}
    with open(bed_fn) as bed:
        for record in bed:
            record = record.split()
            transcript_id = record[3]
            chrom = record[0]
            strand = record[5]
            exstarts, exlengths = parse_exons(record)
            annotation[transcript_id] = (chrom, strand, exstarts, exlengths)
    return annotation


def find_genomic_pos(tx_cut_site, exstarts, exlengths, strand):
    if strand == '-':
        tot_len = exlengths.sum()
        tx_cut_site = tot_len - tx_cut_site
    cumlengths = np.cumsum(exlengths)
    exidx = np.searchsorted(
        cumlengths,
        tx_cut_site, 
    )
    offset = exstarts[exidx] - exlengths[:exidx].sum()
    return tx_cut_site + offset


def find_genomic_inv(tx_start, tx_end, exstarts, exlengths, strand):
    exends = exstarts + exlengths
    start = find_genomic_pos(tx_start, exstarts, exlengths, strand)
    end = find_genomic_pos(tx_end, exstarts, exlengths, strand)
    if strand == '-':
        start, end = end, start
    start_idx = np.searchsorted(exstarts, start, side='right')
    end_idx = np.searchsorted(exstarts, end, side='left')
    if start_idx == end_idx:
        return np.array([[start, end]])
    else:
        invs = [start]
        i = start_idx
        while i != end_idx:
            invs.append(exends[i - 1])
            invs.append(exstarts[i])
            i += 1
        invs.append(end)
        invs = np.array(invs).reshape(-1, 2)
        invs = np.array([i for i in invs if not i[0] == i[1]])
        return invs


def convert_to_genome_coordinates(gstar_alns, annotation, win_size):
    gstar_genomic = []
    for (_, srna_id, transcript_id,
         tx_start, tx_end, tx_cut_site,
         allen_score) in gstar_alns.itertuples():
        chrom, strand, exstarts, exlengths = annotation[transcript_id]
        tlen = sum(exlengths)
        if tx_cut_site > tlen or tx_end > tlen:
            raise ValueError('transcript length in GSTAr and annotation are incompatible')
        cut_site = find_genomic_pos(tx_cut_site, exstarts, exlengths, strand)
        invs = find_genomic_inv(
            max(tx_cut_site - win_size, 0), min(tx_cut_site + win_size, tlen - 1),
            exstarts, exlengths, strand
        )
        start, end = invs[0][0], invs[-1][1]
        gstar_genomic.append([
            chrom, start, end,
            srna_id, transcript_id,
            cut_site, allen_score,
            strand, invs
        ])
    gstar_genomic = pd.DataFrame(gstar_genomic, columns=GSTAR_GENOMIC_COLUMNS)
    return gstar_genomic


def read_gstar_output(gstar_fn, annotation, win_size):
    gstar_alns = pd.read_csv(
        gstar_fn,
        sep='\t',
        names=GSTAR_TRANSCRIPTOMIC_COLUMNS,
    )
    gstar_alns = convert_to_genome_coordinates(gstar_alns, annotation, win_size)
    return gstar_alns


def invs_to_bed12(invs):
    offset = invs[0][0]
    nblocks = len(invs)
    block_starts = ','.join([f'{s - offset:d}' for s in invs[:, 0]])
    block_lengths = ','.join([f'{e - s:d}' for s, e in invs])
    return nblocks, block_starts, block_lengths


def write_bed12(output_fn, gstar_genomic):
    with open(output_fn, 'w') as f:
        for (_, chrom, start, end, srna_id, transcript_id,
             cut_site, allen_score, strand, invs) in  gstar_genomic.itertuples():
            nblocks, block_starts, block_lengths = invs_to_bed12(invs)
            f.write(
                f'{chrom}\t{start}\t{end}\t{srna_id}\t'
                f'{allen_score}\t{strand}\t{cut_site}\t{cut_site}\t'
                f'0,114,178\t{nblocks}\t{block_lengths}\t{block_starts}\n'
            )


@click.command()
@click.option('-a', '--annotation-bed12', required=True)
@click.option('-s', '--GSTAr-sRNAs', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-w', '--window-size', required=False, default=20)
def nanopore_endcut(annotation_bed12, gstar_srnas, output_fn, window_size):
    annotation = read_annotation(annotation_bed12)
    gstar_genomic = read_gstar_output(gstar_srnas, annotation, window_size)
    write_bed12(output_fn, gstar_genomic)


if __name__ == '__main__':
    nanopore_endcut()