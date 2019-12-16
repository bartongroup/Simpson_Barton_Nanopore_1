import numpy as np
import pandas as pd
import pysam
import click


def get_read_mapping_locs(bam_fn):
    mapping = []
    with pysam.AlignmentFile(bam_fn) as bam:
        for aln in bam.fetch():
            mapping.append([
                aln.query_name,
                aln.reference_name,
                aln.reference_start,
                aln.reference_end,
                ['+', '-'][aln.is_reverse],
            ])
    return pd.DataFrame(mapping, columns=['read_id', 'chrom', 'genomic_start', 'genomic_end', 'strand'])


def get_dist_to_prev(read):
    if not read.strand_same_as_prev or not read.chrom_same_as_prev:
        return np.nan
    if read.strand == '+':
        return read.genomic_start_prev_read - read.genomic_end
    elif read.strand == '-':
        return read.genomic_start - read.genomic_end_prev_read


def read_sequencing_summary(ss_fn, bam_fn):
    ss = []
    for ss_fn in ss_fn:
        ss.append(pd.read_csv(ss_fn, sep='\t'))
    ss = pd.concat(ss, axis=0)
    ss = ss.sort_values(['channel', 'start_time'])
    ss['end_time'] = ss['start_time'] + ss['duration']
    ss['prev_read_id'] = ss.groupby('channel').read_id.shift(1)
    ss['prev_end_time'] = ss.groupby('channel').end_time.shift(1)
    ss['missing_signal_time'] = ss.start_time - ss.prev_end_time
    mapped_loc = get_read_mapping_locs(bam_fn)
    ss = ss.merge(mapped_loc, on='read_id', how='left')
    ss = ss.merge(
        mapped_loc,
        left_on='prev_read_id',
        right_on='read_id',
        suffixes=('', '_prev_read'),
        how='left'
    ).drop('read_id_prev_read', axis=1)
    ss['chrom_same_as_prev'] = ss.chrom == ss.chrom_prev_read
    ss['strand_same_as_prev'] = ss.strand == ss.strand_prev_read
    ss['genomic_dist_to_prev'] = ss.apply(get_dist_to_prev, axis=1)
    return ss


@click.command()
@click.option('-b', '--bam-fn', required=True)
@click.option('-s', '--sequencing-summary-fn', required=True, multiple=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-d', '--max-distance', default=1000)
@click.option('-d', '--min-distance', default=-10)
def cli(bam_fn, sequencing_summary_fn, output_fn,
        max_distance, min_distance):
    ss = read_sequencing_summary(sequencing_summary_fn, bam_fn)
    ss['oversplit'] = (
        ss_mapped_only.chrom_same_as_prev & 
        ss_mapped_only.strand_same_as_prev & 
        (ss_mapped_only.genomic_dist_to_prev <= max_distance).fillna(False) &
        (ss_mapped_only.genomic_dist_to_prev >= min_distance).fillna(False)
    )
    ss.to_csv(output_fn, sep='\t')

if __name__ == '__main__':
    cli()