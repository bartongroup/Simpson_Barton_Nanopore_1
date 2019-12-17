import numpy as np
import pandas as pd

import pysam
import click

BINS = [0, 10, 25, 50, 75, 100, 150, 200, 500, np.inf]


def add_pa_tag(bam_fn, output_bam_fn, polya_lengths_fn, bins):
    polya_lengths = pd.read_table(
        polya_lengths_fn, sep='\t',
        usecols=['readname', 'contig', 'position', 'polya_length', 'qc_tag'],
        na_values=['-1'],
        dtype={'readname': str, 'contig': str, 'position': int, 'polya_length': float, 'qc_tag': 'category'},
        index_col='readname'
    )
    polya_lengths = polya_lengths[polya_lengths.qc_tag == 'PASS'].dropna()

    bins_str = np.array([f'{i:04.0f}-{j:04.0f}' for i, j in zip(bins[:-1], bins[1:])])
    polya_lengths['pa_binned'] = bins_str[np.digitize(polya_lengths.polya_length, bins=bins) - 1]
    polya_lengths_dict = polya_lengths.polya_length.to_dict()
    polya_lengths_bins = polya_lengths.pa_binned.to_dict()
    with pysam.AlignmentFile(bam_fn) as inbam, pysam.AlignmentFile(output_bam_fn, 'wb', template=inbam) as outbam:
        for aln in inbam.fetch():
            try:
                pa = polya_lengths_dict[aln.query_name]
                ba = polya_lengths_bins[aln.query_name]
            except KeyError:
                continue
            aln.set_tag('pA', pa, value_type='f')
            aln.set_tag('bA', ba, value_type='Z')
            outbam.write(aln)


@click.command()
@click.option('-b', '--bam-fn', required=True)
@click.option('-o', '--output-bam-fn', required=True)
@click.option('-p', '--nanopolish-polya-tsv-fn', required=True)
def label_bam(bam_fn, output_bam_fn, nanopolish_polya_tsv_fn):
    '''
    Take TSV output from Nanopolish polya and use it to annotate a BAM file
    using pA tag (float, original polyA length from nanopolish) and bA tag
    (str, binned polyA length). Reads which are unmapped or QC filtered by
    nanopolish are removed.
    '''
    add_pa_tag(bam_fn, output_bam_fn, nanopolish_polya_tsv_fn, bins=BINS)


if __name__ == '__main__':
    label_bam()