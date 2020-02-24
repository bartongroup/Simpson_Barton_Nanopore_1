import random
from collections import Counter
import numpy as np
import pandas as pd

import pysam
from joblib import Parallel, delayed

from .io import parse_pysam_aln
from .chimera import (
    identify_chimeras,
)


def iter_bam(bam_fn, query=None):
    if query is None:
        contig = None
        start = None
        end = None
    else:
        contig, start, end = query
    with pysam.AlignmentFile(bam_fn) as bam:
        for aln in bam.fetch(contig, start, end):
            # low quality alignments on repetitive regions
            # can produce dodgy chimeras
            if aln.is_secondary or aln.mapq == 0:
                continue
            yield parse_pysam_aln(aln)


def run_chimerid(*, bam_fn, bed_fn, query=None,
                 blacklisted_genes=None,
                 ovrlp_thresh=0.2,
                 abs_ovrlp_thresh=200):
    if blacklisted_genes is None:
        blacklisted_genes = set()
    chimeric_reads = {}
    with pysam.TabixFile(bed_fn) as bed:
        chimeric_gene_counts = Counter()
        non_chimeric_gene_counts = Counter()
        bam_iter = iter_bam(bam_fn, query=query)
        for res in identify_chimeras(
                bam_iter, bed, blacklisted_genes,
                ovrlp_thresh, abs_ovrlp_thresh):
            is_chimera, gene_id, read_id, record, _ = res
            if is_chimera:
                chimeric_reads[read_id] = record
                chimeric_gene_counts[gene_id] += 1
            else:
                non_chimeric_gene_counts[gene_id] += 1
    chimeric_gene_counts = pd.Series(chimeric_gene_counts)
    non_chimeric_gene_counts = pd.Series(non_chimeric_gene_counts)
    return (chimeric_gene_counts,
            non_chimeric_gene_counts,
            chimeric_reads)


def get_bam_references(bam_fn):
    querys = []
    with pysam.AlignmentFile(bam_fn) as bam:
        for ref in bam.references:
            querys.append((ref, None, None))
    return querys


def parallel_run(n_proc, **chimerid_kwargs):
    chimeric_reads = {}
    res_chimeric = []
    res_non_chimeric = []
    querys = get_bam_references(chimerid_kwargs['bam_fn'])
    n_proc = min(n_proc, len(querys))
    with Parallel(n_jobs=n_proc) as pool:
        res = pool(
            delayed(run_chimerid)(query=q, **chimerid_kwargs)
            for q in querys
        )
        for r in res:
            (chimeric_gene_counts,
             non_chimeric_gene_counts,
             reads) = r
            chimeric_reads.update(reads)
            res_chimeric.append(chimeric_gene_counts)
            res_non_chimeric.append(non_chimeric_gene_counts)
    chimeric_gene_counts = pd.concat(res_chimeric, axis=0)
    non_chimeric_gene_counts = pd.concat(res_non_chimeric, axis=0)
    return chimeric_gene_counts, non_chimeric_gene_counts, chimeric_reads


def bootstrap(chimeric, nonchimeric, sample_frac):
    norm_factor = chimeric.sum() + nonchimeric.sum()
    sample_size = int(norm_factor * sample_frac)
    idx = np.arange(len(chimeric) + len(nonchimeric))
    weights = np.concatenate([chimeric, nonchimeric]) / norm_factor
    s = np.random.choice(idx, size=sample_size, p=weights, replace=True)
    s = np.bincount(s, minlength=len(idx))
    chimeric_sample = pd.Series(s[:len(chimeric)], index=chimeric.index)
    nonchimeric_sample = pd.Series(s[len(chimeric):], index=nonchimeric.index)
    return sample_size, chimeric_sample, nonchimeric_sample


def parallel_bootstrap(n_boots, n_proc, sample_frac=1,
                       **chimerid_kwargs):
    (chimeric_gene_counts, non_chimeric_gene_counts,
     chimeric_reads) = parallel_run(n_proc, **chimerid_kwargs)
    
    boot_chimeric_gene_counts = {}
    boot_nonchimeric_gene_counts = {}
    boot_norm_factors = {}
    norm_factors = {}
    for i in range(n_boots):
        sample_size, c_sample, nc_sample = bootstrap(
            chimeric_gene_counts,
            non_chimeric_gene_counts,
            sample_frac,
        )
        boot_chimeric_gene_counts[i] = c_sample
        boot_nonchimeric_gene_counts[i] = nc_sample
        norm_factors[i] = sample_size
    boot_chimeric_gene_counts = pd.concat(
        boot_chimeric_gene_counts, axis=1, sort=True
    )
    boot_nonchimeric_gene_counts = pd.concat(
        boot_nonchimeric_gene_counts, axis=1, sort=True
    )
    norm_factors = pd.Series(norm_factors)
    return (boot_chimeric_gene_counts, boot_nonchimeric_gene_counts,
            chimeric_reads, norm_factors)