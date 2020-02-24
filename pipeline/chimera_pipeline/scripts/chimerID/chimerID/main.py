import re
import click

from .reference import create_flattened_reference_bed
from .bootstrap import parallel_bootstrap
from .io import to_bed12
from .logodds import generate_bootstrapped_logodds

def parse_query(query):
    try:
        chrom, start, end = re.search('^(.+):(\d+)-(\d+)$', query).groups()
    except ValueError:
        raise ValueError('Could not parse query')
    start, end = int(start), int(end)
    return chrom, start, end


@click.group()
def chimerID():
    pass


@chimerID.command()
@click.option('-b', '--bam', required=True)
@click.option('-g', '--gtf', required=True)
@click.option('-o', '--output-prefix', required=True)
@click.option('-bl', '--blacklisted-genes', required=False, default=None)
@click.option('-p', '--processes', default=1)
@click.option('-n', '--n-boots', default=100)
@click.option('-s', '--sample-frac', default=0.75)
@click.option('-ggo', '--gene-gene-merge-overlap-threshold', type=float, default=0.75)
@click.option('-rgo', '--read-gene-overlap-threshold', type=float, default=0.2)
@click.option('-argo', '--abs-read-gene-overlap-threshold', type=int, default=200)
def detect_chimeras(bam, gtf, output_prefix, blacklisted_genes,
                    processes, n_boots, sample_frac,
                    gene_gene_merge_overlap_threshold,
                    read_gene_overlap_threshold,
                    abs_read_gene_overlap_threshold):
    bed_fn = create_flattened_reference_bed(
        gtf, output_prefix,
        overlap_threshold=gene_gene_merge_overlap_threshold
    )
    if blacklisted_genes is not None:
        with open(blacklisted_genes) as bl:
            blacklisted_genes = set([b.strip() for b in bl])
    else:
        blacklisted_genes = set()
    chimeric_counts, non_chimeric_counts, reads, norm_factors = parallel_bootstrap(
        n_boots=n_boots, n_proc=processes,
        bam_fn=bam,
        bed_fn=bed_fn,
        blacklisted_genes=blacklisted_genes,
        sample_frac=sample_frac,
        ovrlp_thresh=read_gene_overlap_threshold,
        abs_ovrlp_thresh=abs_read_gene_overlap_threshold,
    )
    chimeric_counts.to_hdf(f'{output_prefix}.boot_gp_counts.h5', key='chimera_counts')
    non_chimeric_counts.to_hdf(f'{output_prefix}.boot_gp_counts.h5', key='non_chimeric_counts', mode='r+')
    norm_factors.to_hdf(f'{output_prefix}.boot_gp_counts.h5', key='norm_factors', mode='r+')
    with open(f'{output_prefix}.chimeric_reads.bed', 'w') as f:
        for read_id, (invs, chrom, strand) in reads.items():
            f.write(to_bed12(chrom, invs, strand, read_id))


@chimerID.command()
@click.option('-a', '--cond-a-h5-fn', required=True)
@click.option('-b', '--cond-b-h5-fn', required=True)
@click.option('-as', '--cond-a-name', required=True)
@click.option('-bs', '--cond-b-name', required=True)
@click.option('-o', '--output-fn', required=True)
def bootstrap_logodds(cond_a_h5_fn, cond_b_h5_fn, cond_a_name, cond_b_name, output_fn):
    h5_fns = [cond_a_h5_fn, cond_b_h5_fn]
    logodds_data = generate_bootstrapped_logodds(h5_fns, cond_a_name, cond_b_name)
    logodds_data.to_csv(output_fn, sep='\t')
    

if __name__ == '__main__':
    cli()