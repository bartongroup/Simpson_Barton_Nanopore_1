import pysam
import click


RC = str.maketrans('ACGTN', 'TGCAN')


def rev_comp(seq):
    return seq.translate(RC)[::-1]


def get_five_prime_softclipped(aln):
    if not aln.is_reverse:
        op, ln = aln.cigar[0]
        seq = aln.query_sequence
    else:
        op, ln = aln.cigar[-1]
        seq = rev_comp(aln.query_sequence)
    if op == 4:
        sc = seq[:ln]
    else:
        sc = ''
    return sc


def filter_bam_by_5p_g_softclip(inbam, capped_bam, uncapped_bam):
    for aln in inbam.fetch():
        softclipped_5p = get_five_prime_softclipped(aln)
        if softclipped_5p and softclipped_5p[-1] == 'G':
            capped_bam.write(aln)
        else:
            uncapped_bam.write(aln)


@click.command()
@click.option('-i', '--input-bam-fn')
@click.option('-c', '--capped-bam-fn')
@click.option('-u', '--uncapped-bam-fn')
def filter_capped(input_bam_fn, capped_bam_fn, uncapped_bam_fn):
    with pysam.AlignmentFile(input_bam_fn) as inbam:
        with pysam.AlignmentFile(capped_bam_fn, mode='wb', template=inbam) as capped_bam, \
             pysam.AlignmentFile(uncapped_bam_fn, mode='wb', template=inbam) as uncapped_bam:
            filter_bam_by_5p_g_softclip(inbam, capped_bam, uncapped_bam)


if __name__ == '__main__':
    filter_capped()