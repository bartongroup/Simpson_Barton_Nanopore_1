import sys
import re
from collections import defaultdict
import itertools as it

import numpy as np
import click


ALTERNATIVE_JUNCTION_LABELS = {
    'left': {'+': 'alternate_donor', '-': 'alternate_acceptor'},
    'right': {'+': 'alternate_acceptor', '-': 'alternate_donor'},
}

ALTERNATIVE_END_LABELS = {
    'left': {'+': 'five_p_var', '-': 'three_p_var'},
    'right': {'+': 'three_p_var', '-': 'five_p_var'},
}


BED_RECORD = (
    '{chrom}\t{start:d}\t{end:d}\t'
    '{name}_{ctype}\t60\t{strand}\n'
)


def exons_to_introns(exon_invs):
    if len(exon_invs) == 1:
        return []
    else:
        exon_invs = sorted(exon_invs)
        intron_invs = list(zip(
            [ex[1] for ex in exon_invs[:-1]],
            [ex[0] for ex in exon_invs[1:]]
        ))
        return intron_invs


def fix_complex(chunk_labels, strand):
    ln = len(chunk_labels)
    for i, label in enumerate(chunk_labels):
        if label == 'complex':
            if i == 0:
                chunk_labels[i] = ALTERNATIVE_END_LABELS['left'][strand]
            elif i == ln - 1:
                chunk_labels[i] = ALTERNATIVE_END_LABELS['right'][strand]
            elif chunk_labels[i - 1] == 'intron' and chunk_labels[i + 1] == 'intron':
                chunk_labels[i] = 'partial_exon'
            elif chunk_labels[i - 1] == 'exon' and chunk_labels[i + 1] == 'exon':
                chunk_labels[i] = 'partial_intron'
            elif chunk_labels[i - 1] == 'alternate_donor' and chunk_labels[i + 1] == 'alternate_acceptor':
                chunk_labels[i] = 'partial_intron' if strand == '+' else 'partial_exon'
            elif chunk_labels[i - 1] == 'alternate_acceptor' and chunk_labels[i + 1] == 'alternate_donor':
                chunk_labels[i] = 'partial_exon' if strand == '+' else 'partial_intron'
            elif chunk_labels[i - 1] == 'exon' and (
                    chunk_labels[i + 1].startswith('alternate') or chunk_labels[i + 1].endswith('var')):
                chunk_labels[i] = 'partial_intron'
            elif chunk_labels[i + 1] == 'exon' and (
                    chunk_labels[i - 1].startswith('alternate') or chunk_labels[i - 1].endswith('var')):
                chunk_labels[i] = 'partial_intron'
            elif chunk_labels[i - 1] == 'intron' and (
                    chunk_labels[i + 1].startswith('alternate') or chunk_labels[i + 1].endswith('var')):
                chunk_labels[i] = 'partial_exon'
            elif chunk_labels[i + 1] == 'intron' and (
                    chunk_labels[i - 1].startswith('alternate') or chunk_labels[i - 1].endswith('var')):
                chunk_labels[i] = 'partial_exon'
    return chunk_labels


def get_chunks(exon_invs, strand):
    starts = set([invs[0][0] for invs in exon_invs.values()])
    ends = set([invs[-1][1] for invs in exon_invs.values()])
    first_exons = set([invs[0] for invs in exon_invs.values()])
    last_exons = set([invs[-1] for invs in exon_invs.values()])
    all_invs = list(it.chain(*exon_invs.values()))
    chunks = np.unique(np.array(all_invs).ravel())
    chunks = np.stack([chunks[:-1], chunks[1:]], axis=1)
    exons = set(it.chain(*exon_invs.values()))
    exon_starts = set([inv[0] for inv in exons])
    exon_ends = set([inv[1] for inv in exons])
    intron_invs = {t_id: exons_to_introns(ex) for t_id, ex in exon_invs.items()}
    introns = set(it.chain(*intron_invs.values()))
    chunk_labels = []
    for chunk in chunks:
        chunk = tuple(chunk)
        if chunk in introns:
            chunk_labels.append('intron')
        elif chunk in exons:
            if chunk in first_exons:
                chunk_labels.append('first_exon' if strand == '+' else 'terminal_exon')
            elif chunk in last_exons:
                chunk_labels.append('first_exon' if strand == '-' else 'terminal_exon')
            else:
                for intron in introns:
                    if intron[0] < chunk[1] and intron[1] > chunk[1]:
                        chunk_labels.append('cassette_exon')
                        break
                else:
                    chunk_labels.append('constitutive_exon')
        elif starts.intersection(chunk):
            chunk_labels.append(ALTERNATIVE_END_LABELS['left'][strand])
        elif ends.intersection(chunk):
            chunk_labels.append(ALTERNATIVE_END_LABELS['right'][strand])
        elif exon_ends.issuperset(chunk):
            chunk_labels.append(ALTERNATIVE_JUNCTION_LABELS['left'][strand])
        elif exon_starts.issuperset(chunk):
            chunk_labels.append(ALTERNATIVE_JUNCTION_LABELS['right'][strand])
        else:
            chunk_labels.append('complex')
    chunk_labels = fix_complex(chunk_labels, strand)
    return chunks, chunk_labels


def get_gtf_gene_id(attrs):
    return re.search('gene_id \"(.*?)\";', attrs).group(1)


def get_gtf_transcript_id(attrs):
    return re.search('transcript_id \"(.*?)\";', attrs).group(1)


def get_gtf_exons(gtf_fn):
    with open(gtf_fn) as gtf:
        for record in gtf:
            record = record.strip().split('\t')
            if record[2] == 'exon':
                gene_id = get_gtf_gene_id(record[8])
                transcript_id = get_gtf_transcript_id(record[8])
                yield record[0], int(record[3]) - 1, int(record[4]), gene_id, transcript_id, record[6]


def parse_gtf_chunk_invs(gtf_fn):
    gene_cluster = defaultdict(list)
    gtf_iter = get_gtf_exons(gtf_fn)
    curr_chrom, start, end, curr_gene_id, transcript_id, curr_strand = next(gtf_iter)
    gene_cluster[transcript_id].append((start, end))
    for chrom, start, end, gene_id, transcript_id, strand in gtf_iter:
        if gene_id != curr_gene_id:
            chunks, chunk_labels = get_chunks(gene_cluster, curr_strand)
            yield curr_gene_id, curr_chrom, curr_strand, chunks, chunk_labels
            curr_gene_id, curr_chrom, curr_strand = gene_id, chrom, strand
            gene_cluster = defaultdict(list)
            gene_cluster[transcript_id].append((start, end))
        else:
            gene_cluster[transcript_id].append((start, end))
    if gene_cluster:
        chunks, chunk_labels = get_chunks(gene_cluster, curr_strand)
        yield curr_gene_id, curr_chrom, curr_strand, chunks, chunk_labels


def to_bed(chrom, start, end, strand, name, chunk_type):
    return BED_RECORD.format(
        chrom=chrom, start=start, end=end,
        name=name, ctype=chunk_type, strand=strand
    )


@click.command()
@click.option('-g', '--gtf', required=True)
def cli(gtf):
    for gene_id, chrom, strand, chunks, chunk_labels in parse_gtf_chunk_invs(gtf):
        for (start, end), chunk_type in zip(chunks, chunk_labels):
            sys.stdout.write(to_bed(chrom, start, end, strand, gene_id, chunk_type))


if __name__ == '__main__':
    cli()