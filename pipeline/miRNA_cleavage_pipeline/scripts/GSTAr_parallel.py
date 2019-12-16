import os
import shutil
import subprocess
import tempfile
from io import StringIO

import pandas as pd
import joblib
import click


GSTAR_PATH = os.path.join(
    os.path.split(
        os.path.abspath(__file__)
    )[0],
    'GSTAr.pl'
)

GSTAR_COLUMNS = [
    'srna_id', 'transcript_id',
    'start', 'end', 'cut_site',
    'mfe_perfect_match', 'mfe_actual', 'mfe_ratio',
    'allen_score', 'paired', 'unpaired',
    'structure', 'sequence'
]

GSTAR_USECOLS = [
    'srna_id', 'transcript_id',
    'start', 'end', 'cut_site', 'allen_score'
]


def chunk_fasta(fasta_fn, n_chunks):
    # only works for single line fasta
    n_records = 0
    with open(fasta_fn) as fasta:
        for record in fasta:
            if record[0] == '>':
                n_records += 1
    n_seqs_per_chunk, r = divmod(n_records, n_chunks)
    with open(fasta_fn) as fasta:
        chunk = []
        n_records = 0
        chunks_created = 0
        while True:
            try:
                name = next(fasta)
            except StopIteration:
                break
            seq = next(fasta)
            chunk.append(name)
            chunk.append(seq)
            n_records += 1
            if chunks_created < r and n_records == n_seqs_per_chunk + 1:
                yield chunk
                chunk = []
                n_records = 0
                chunks_created += 1
            elif chunks_created >= r and n_records == n_seqs_per_chunk:
                yield chunk
                chunk = []
                n_records = 0
                chunks_created += 1


def read_gstar_output(gstar_fn):
    '''
    read the tsv file format produced by GSTAr in -t mode
    '''
    gstar = pd.read_csv(
        gstar_fn,
        sep='\s+',
        names=GSTAR_COLUMNS,
        usecols=GSTAR_USECOLS,
        skiprows=8,
        comment='#'
    )
    gstar['start'] = gstar['start'] - 1
    gstar['cut_site'] = gstar['cut_site'] - 1
    return gstar


def parallel_GSTAr(fasta_chunk, transcriptome_fn, tmpdir):
    # parallel GSTAr instances have to be run in their own separate directories
    with tempfile.TemporaryDirectory(dir=tmpdir) as chunk_dir:
        # write the fasta chunk to the temporary directory
        srna_chunk_fn = os.path.join(chunk_dir, 'sRNAs.fa')
        with open(srna_chunk_fn, 'w') as fasta:
            for line in fasta_chunk:
                fasta.write(line)
        # run GSTAr in the temporary directory
        proc = subprocess.Popen(
            ['perl', GSTAR_PATH, '-t', srna_chunk_fn, transcriptome_fn],
            cwd=chunk_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = proc.communicate()
        return read_gstar_output(StringIO(stdout.decode()))
        

@click.command()
@click.option('-s', '--srna-fn', required=True)
@click.option('-t', '--transcriptome-fn', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-n', '--n-jobs', required=False, default=4)
def run_GSTAr(srna_fn, transcriptome_fn, output_fn, n_jobs):
    tmpdir = os.environ['TMPDIR']
    transcriptome_base_fn = os.path.split(transcriptome_fn)[1]
    transcriptome_copy = os.path.join(tmpdir, transcriptome_base_fn)
    shutil.copyfile(transcriptome_fn, transcriptome_copy)

    with joblib.Parallel(n_jobs) as pool:
        res = pool(
            joblib.delayed(parallel_GSTAr)(fasta_chunk, transcriptome_copy, tmpdir)
            for fasta_chunk in chunk_fasta(srna_fn, n_jobs)
        )
    res = pd.concat(res, axis=0)
    res.to_csv(output_fn, sep='\t', index=False, header=False)


if __name__ == '__main__':
    run_GSTAr()