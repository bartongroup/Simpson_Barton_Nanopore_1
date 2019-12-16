import re

import numpy as np
import pandas as pd
import pysam
import pybedtools

import click


DER_COLS = ['chrom', 'start', 'end', 'comp', 'logFC', 'strand',
            'meanCoverage', 'mean_cntrl', 'mean_treat', 'p_val', 'fdr']
PRIORITY = ['three_prime_UTR', 'five_prime_UTR', 'CDS',
            'miRNA_primary_transcript', 'rRNA', 'tRNA',
            'snoRNA', 'snRNA', 'lnc_RNA', 'antisense_lncRNA', 'ncRNA',
            'pseudogenic_exon', 'transposable_element', 'transcript_region']


def annotate_inv(genes, chrom, start, end, strand):
    feats = set()
    parents = set()
    antisense_feats = set()
    antisense_parents = set()
    overlapping_feats = list(genes.fetch(chrom, start, end, parser=pysam.asGTF()))
    if not len(overlapping_feats):
        return None
    else:
        for record in overlapping_feats:
            if record.strand == strand:
                feats.add(record.feature)
                try:
                    parents.add(re.search(r'Parent=(.*?)[;:\.]', record.attributes).group(1))
                except AttributeError:
                    parents.add(re.search(r'ID=(.*?)[;:\\.]', record.attributes).group(1))
            else:
                antisense_feats.add(record.feature)
                try:
                    antisense_parents.add(re.search(r'Parent=(.*?)[;:\.]', record.attributes).group(1))
                except AttributeError:
                    antisense_parents.add(re.search(r'ID=(.*?)[;:\\.]', record.attributes).group(1))
    if not len(feats) and len(antisense_feats):
        feats = 'antisense'
        parents = antisense_parents
    else:
        feats = feats.difference({'gene', 'exon', 'protein', 'mRNA'})
        if not feats:
            feats = 'intron'
        else:
            for ftype in PRIORITY:
                if ftype in feats:
                    feats = ftype
                    break
            else:
                feats = '|'.join(feats)
    return [feats, '|'.join(parents)]


def annotate_derfinder(derfinder_tsv_fn, genes_fn):
    df_res = pd.read_csv(
        derfinder_tsv_fn,
        sep='\t',
        names=DER_COLS
    )
    genic_der_records = []
    intergenic_der_records = []
    with pysam.TabixFile(genes_fn) as genes:
        for record in df_res.itertuples(index=False):
            chrom, start, end, _, _, strand, *_ = record
            annot = annotate_inv(
                genes, chrom, start, end, strand)
            if annot is None:
                intergenic_der_records.append(record)
            else:
                record = list(record) + annot
                genic_der_records.append(record)
    genic_der_records = pd.DataFrame(
        genic_der_records, columns=DER_COLS + ['feature_type', 'gene_ids'])
    genic_der_records['upstream'] = np.nan
    genic_der_records['downstream'] = np.nan
    # Use pybedtools to find closest genes for intergenic DERs
    genes = pybedtools.BedTool(genes_fn).filter(lambda record: "Parent" not in record[8]).saveas()
    intergenic_bt = pybedtools.BedTool.from_dataframe(pd.DataFrame(intergenic_der_records))
    upstream = intergenic_bt.closest(genes, D='a', id=True).to_dataframe(header=None)
    downstream = intergenic_bt.closest(genes, D='a', iu=True).to_dataframe(header=None)
    upstream_genes = upstream.iloc[:, -2].str.extract('ID=(.*?)[;:\.]')
    downstream_genes = downstream.iloc[:, -2].str.extract('ID=(.*?)[;:\.]')

    intergenic_der_records = pd.DataFrame(intergenic_der_records, columns=DER_COLS)
    intergenic_der_records['feature_type'] = 'intergenic'
    intergenic_der_records['gene_ids'] = np.nan
    intergenic_der_records['upstream'] = upstream_genes
    intergenic_der_records['downstream'] = downstream_genes
    der_records_annot = pd.concat(
        [genic_der_records, intergenic_der_records]
    ).reset_index(drop=True)
    return der_records_annot


def format_gtf(record, name):
    GTF_RECORD = ('{chrom}\t{source}\t{feature}\t'
                  '{start:d}\t{end:d}\t.\t{strand}\t.\t'
                  'gene_id {name}; transcript_id {name}; feature_type {feature}; '
                  'overlaps {overlaps}; upstream {upstream}; downstream {downstream}; '
                  'logFC {logFC:.3f}; fdr {fdr:.3g}; mean_cov {cov:.2f};\n')
    fmtted = GTF_RECORD.format(
        chrom=record.chrom, source=record.comp, feature=record.feature_type,
        start=record.start, end=record.end, strand=record.strand,
        name=name, overlaps=record.gene_ids,
        upstream=record.upstream, downstream=record.downstream,
        logFC=record.logFC, fdr=record.fdr, cov=record.meanCoverage
    )
    return fmtted


@click.command()
@click.option('-b', '--derfinder-bed-fn', required=True)
@click.option('-g', '--gff-fn', required=True, help='Araport11 GFF')
@click.option('-o', '--output-fn', required=True)
def cli(derfinder_bed_fn, gff_fn, output_fn):
    annotated = annotate_derfinder(derfinder_bed_fn, gff_fn)
    annotated = annotated.sort_values(['chrom', 'start'])
    with open(output_fn, 'w') as f:
        for i, record in enumerate(annotated.itertuples()):
            f.write(format_gtf(record, f'derfinder_region_{i}'))


if __name__ == '__main__':
    cli()