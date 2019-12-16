import sys
import os
from glob import glob
import random
import tempfile

import numpy as np
import pandas as pd
from scipy import stats

import pyBigWig
import click

from rpy2 import robjects as robj
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri


pandas2ri.activate()


derfinder = importr('derfinder')
robj.r['options'](species='arabidopsis_thaliana')
robj.r['options'](chrsStyle=robj.r("NULL"))


def get_chroms_list(bw_fn):
    with pyBigWig.open(bw_fn) as bw:
        chroms = list(bw.chroms().keys())
    return chroms


def load_data(file_names, sample_names, conditions,
              chroms=None, cutoff=10, **fullcov_kws):
    chroms = robj.StrVector(chroms)
    fn = robj.StrVector(file_names)
    sn = robj.StrVector(sample_names)
    cn = robj.StrVector(conditions)
    cond_data = robj.DataFrame({'row.names': sn, 'cond': cn}, stringsasfactor=True)
    fn.names = sn
    full_cov = derfinder.fullCoverage(
        files=fn, chrs=chroms,
        cutoff=cutoff, **fullcov_kws
    )
    return full_cov, cond_data


def make_models(full_cov, cond_data):
    sample_depths = derfinder.sampleDepth(
        derfinder.collapseFullCoverage(full_cov), 1)
    models = derfinder.makeModels(sample_depths, testvars=cond_data.rx('cond')[0])
    return models


def run_base_level_derfinder(full_cov, cond_data, models, chrom, n_permutations, seed):
    seeds = robj.IntVector([seed + i for i in range(n_permutations)])
    tmpdir = tempfile.mkdtemp(dir='.')
    res = derfinder.analyzeChr(
        chr=robj.StrVector([chrom]), coverageInfo=full_cov.rx(chrom)[0],
        models=models, groupInfo=cond_data.rx('cond')[0], writeOutput=False,
        cutoffFstat=5e-02, nPermute=n_permutations, seeds=seeds,
        returnOutput=True, runAnnotation=False, lowMemDir=os.path.join(tmpdir, chrom, 'chunksDir')
    )
    res = res.rx('regions')[0].rx('regions')[0]
    res = robj.r['as.data.frame'](res)
    res = robj.r['type.convert'](res, **{'as.is':True})
    res = pandas2ri.ri2py(res)
    return res


def run_derfinder(cntrl_bws, treat_bws,
                  cntrl_name, treat_name,
                  chroms=None, expression_cutoff=10,
                  n_permutations=50, seed=10123):
    assert cntrl_name != treat_name
    sample_names = ['{}_{}'.format(cntrl_name, i) for i in range(len(cntrl_bws))] + \
                   ['{}_{}'.format(treat_name, i) for i in range(len(treat_bws))] 
    conds = [s.split('_')[0] for s in sample_names]
    bws = cntrl_bws + treat_bws
    if chroms is None:
        # find ALL the chroms:
        chroms = get_chroms_list(bws[0])
    elif isinstance(chroms, str):
        chroms = [chroms,]
    all_chrom_results = []
    full_cov, cond_data = load_data(
        bws, sample_names, conds,
        chroms=chroms,
        cutoff=expression_cutoff
    )
    models = make_models(full_cov, cond_data)
    for chrom in chroms:
        chrom_res = run_base_level_derfinder(
            full_cov, cond_data, models, chrom, n_permutations, seed
        )
        all_chrom_results.append(chrom_res)
    results = pd.concat(all_chrom_results, axis=0)
    return results


def write_bed_of_sig_res(results, output_file, cntrl_name, treat_name, effect_size_threshold=1):
    name = '{}_vs_{}'.format(cntrl_name, treat_name)
    fch_col = results.columns[results.columns.str.startswith('log2FoldChange')][0]
    results_filt = results.copy()
    results_filt['name'] = name
    results_filt = results_filt.loc[
        results.significantQval.astype(bool) & (np.abs(results[fch_col]) > effect_size_threshold),
        ['seqnames', 'start', 'end', 'name', fch_col, 'strand', 'meanCoverage',
         f'mean{cntrl_name}', f'mean{treat_name}', 'pvalues', 'qvalues']
    ]
    results_filt = results_filt.sort_values(['seqnames', 'start'])
    results_filt['pvalues'] = -np.log10(results_filt.pvalues)
    results_filt['qvalues'] = -np.log10(results_filt.qvalues)
    results_filt.to_csv(output_file, sep='\t', header=False, index=False, float_format='%5g')


@click.command()
@click.option('-cb', '--cntrl-bws', multiple=True, required=True)
@click.option('-tb', '--treat-bws', multiple=True, required=True)
@click.option('-cn', '--cntrl-name', required=True)
@click.option('-tn', '--treat-name', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-b', '--bed-output-fn', required=True)
@click.option('-s', '--strand', required=False, type=click.Choice(['+', '-', '.']))
@click.option('-est', '--effect-size-threshold', required=False, default=1, type=float)
@click.option('-c', '--expression-cutoff', required=False, default=10, type=float)
@click.option('-np', '--n-permutations', required=False, default=50, type=int)
def cli(cntrl_bws, treat_bws,
        cntrl_name, treat_name,
        output_fn, bed_output_fn, strand,
        effect_size_threshold,
        expression_cutoff,
        n_permutations):
    res = run_derfinder(
        cntrl_bws, treat_bws,
        cntrl_name, treat_name,
        n_permutations=n_permutations,
        expression_cutoff=expression_cutoff
    )
    res['strand'] = strand
    res.to_csv(output_fn, sep='\t')
    write_bed_of_sig_res(res, bed_output_fn,
                         cntrl_name, treat_name,
                         effect_size_threshold)


if __name__ == '__main__':
    cli()