import os
import numpy as np
from multiprocessing import cpu_count
import logging

from joblib import Parallel, delayed
import pysam
import click
import click_log

from .utils import get_fast5_read_id_mapping, get_output_filenames
from .model import load_adapter_model
from .signal import get_fast5_fiveprime
from .trim import Trimmer


logger = logging.getLogger(__name__)
click_log.basic_config(logger)


def predict_bam_batch(batch, reads, 
                      model, threshold,
                      batch_size, squiggle_size):
    '''
    Make a prediction for a batch of signals
    '''
    read_ids, sig_len, batch = zip(*batch)
    sig_len = np.array(sig_len)
    batch = np.array(batch).reshape(batch_size, squiggle_size, 1)
    batch = batch[~np.isnan(sig_len)]
    read_ids = [r_id for r_id, sl in zip(read_ids, sig_len) if not np.isnan(sl)]
    preds = model.predict(batch).squeeze()
    logger.info('INFO: Made prediction for batch')
    for score, read_id, read in zip(preds, read_ids, reads):
        assert read_id == read.query_name
        read.set_tag('ac', score, value_type='f')
        yield score > threshold, read


def predict_fq_batch(batch, reads, 
                     model, trimmer, threshold,
                     batch_size, squiggle_size):
    '''
    Make a prediction for a batch of signals
    '''
    read_ids, _, batch = zip(*batch)
    batch = np.array(batch).reshape(batch_size, squiggle_size, 1)
    preds = model.predict(batch).squeeze()
    reads = trimmer(reads, preds > threshold)
    logger.info('INFO: Made prediction for batch')
    for score, read_id, read in zip(preds, read_ids, reads):
        assert read_id == read.name
        yield score > threshold, read


def bam_batch_iter(bam, read_id_fast5_filemap,
                   batch_size, keep_unmapped):
    batch = []
    batch_reads = []
    for read in bam.fetch(until_eof=keep_unmapped):
        read_id = read.query_name
        try:
            f5_fn = read_id_fast5_filemap[read_id]
        except KeyError:
            logger.warning(
                'WARNING: Fast5 file not found for read {} in bam file'.format(read_id))
            continue
        batch.append((read_id, f5_fn))
        batch_reads.append(read)
        if len(batch) == batch_size:
            logger.info('INFO: Identified batch of read fast5 files')
            yield batch, batch_reads, batch_size
            batch = []
            batch_reads = []
    if len(batch):
        yield batch, batch_reads, len(batch)


def fq_batch_iter(fq, read_id_fast5_filemap, batch_size):
    batch = []
    batch_reads = []
    for read in fq:
        read_id = read.name
        try:
            f5_fn = read_id_fast5_filemap[read_id]
        except KeyError:
            logger.warning(
                'WARNING: Fast5 file not found for read {} in bam file'.format(read_id))
            continue
        batch.append((read_id, f5_fn))
        batch_reads.append(read)
        if len(batch) == batch_size:
            logger.info('INFO: Identified batch of read fast5 files')
            yield batch, batch_reads, batch_size
            batch = []
            batch_reads = []
    if len(batch):
        yield batch, batch_reads, len(batch)


def bam_filter(bam_fn, read_id_fast5_filemap, model,
               pass_output_bam_fn, fail_output_bam_fn, threshold=0.5,
               batch_size=5000, squiggle_size=2000,
               keep_unmapped=False, processes=8):
    '''
    Filter a bam file using a model trained to detect adapters in
    signal from 5' end.
    
    Parameters:
    ----------
    
        bam_fn: str, required
          Path to input bam file

        read_id_fast5_filemap: dict, required
          Dict of read_id: fast5 filepath key value pairs

        model: keras.models.Model, required
          Model used to make predictions. Must accept input of
          shape (batch_size, squiggle_size, 1)

        pass_output_bam_fn: str, required
          Bam file path to write reads passing filter to.
          Will be overwritten

        fail_output_bam_fn: str, required
          Bam file path to write reads failing filter to.
          Will be overwritten

        threshold: float, optional, default: 0.5
          Threshold at which to cutoff pass/fail

        batch_size, int, optional, default: 5000
          Number of signals to pass to model prediction
          at a time

        squiggle_size: int, optional, default: 2000
          Length of signal to take from 5' end of RNA fast5 signal,
          to pass to model.

        keep_unmapped: bool, optional, default: False
          Whether to process unmapped reads in the bam file. If false,
          unmapped reads are filtered out (i.e. are not present in
          either pass or fail bam files). This speeds up computation.

        processed: int, optional, default: 8
          Number of parallel processes to use to read signals from fast5

    Returns:
    --------
    
        None
    '''
    in_bam = pysam.AlignmentFile(bam_fn, mode='rb')
    # create two output bams
    out_bams = {
        True: pysam.AlignmentFile(pass_output_bam_fn,
                                  mode='wb', template=in_bam),
        False: pysam.AlignmentFile(fail_output_bam_fn,
                                   mode='wb', template=in_bam)
    }
    filter_counts = {True: 0, False: 0}
    with Parallel(n_jobs=processes) as pool:
        for batch, reads, bs in bam_batch_iter(in_bam, read_id_fast5_filemap,
                                               batch_size, keep_unmapped):
            # not sure if joblib guarantees the output order is the same as
            # input order, so pass the read_ids as well to be safe
            batch_signals = pool(
                delayed(get_fast5_fiveprime)(*f, squiggle_size) for f in batch)
            logger.info('INFO: Collected Five prime signals for batch')
            for decision, read in predict_bam_batch(batch_signals, reads, model,
                                                    threshold, bs, squiggle_size):
                out_bams[decision].write(read)
                filter_counts[decision] += 1
            logger.info('INFO: Wrote batch to bam files')
            logger.info('INFO: Processed {} reads'.format(sum(filter_counts.values())))
            logger.info('INFO: Of these:')
            logger.info('INFO:   {} passed the threshold'.format(filter_counts[True]))
            logger.info('INFO:   {} failed the threshold'.format(filter_counts[False]))

            
    in_bam.close()
    for bam in out_bams.values():
        bam.close()


def fastq_filter(fastq_fn, read_id_fast5_filemap, model,
                 pass_output_fq_fn, fail_output_fq_fn, threshold=0.5,
                 batch_size=5000, squiggle_size=2000,
                 trim=False, processes=8):
    '''
    Filter a bam file using a model trained to detect adapters in
    signal from 5' end.
    
    Parameters:
    ----------
    
        fastq_fn: str, required
          Path to input fastq file

        read_id_fast5_filemap: dict, required
          Dict of read_id: fast5 filepath key value pairs

        model: keras.models.Model, required
          Model used to make predictions. Must accept input of
          shape (batch_size, squiggle_size, 1)

        pass_output_fq_fn: str, required
          Fastq file path to write reads passing filter to.
          Will be overwritten

        fail_output_fq_fn: str, required
          Fastq file path to write reads failing filter to.
          Will be overwritten

        threshold: float, optional, default: 0.5
          Threshold at which to cutoff pass/fail

        batch_size, int, optional, default: 5000
          Number of signals to pass to model prediction
          at a time

        squiggle_size: int, optional, default: 2000
          Length of signal to take from 5' end of RNA fast5 signal,
          to pass to model.

        trim: bool, optional, default: False
          If True, a second model trained at the sequence level will
          be used to estimate the length of the adapter and cut it off

        processed: int, optional, default: 8
          Number of parallel processes to use to read signals from fast5

    Returns:
    --------
    
        None
    '''
    in_fq = pysam.FastxFile(fastq_fn)
    # create two output fastqs
    out_fqs = {
        True: open(pass_output_fq_fn, mode='w'),
        False: open(fail_output_fq_fn, mode='w')
    }
    if trim:
        trimmer = Trimmer().trim_adapters
    else:
        trimmer = lambda x, y: x
    filter_counts = {True: 0, False: 0}
    with Parallel(n_jobs=processes) as pool:
        for batch, reads, bs in fq_batch_iter(in_fq, read_id_fast5_filemap,
                                              batch_size):
            # not sure if joblib guarantees the output order is the same as
            # input order, so pass the read_ids as well to be safe
            batch_signals = pool(
                delayed(get_fast5_fiveprime)(*f, squiggle_size) for f in batch)
            logger.info('INFO: Collected Five prime signals for batch')
            for decision, read in predict_fq_batch(batch_signals, reads, model, trimmer,
                                                   threshold, bs, squiggle_size):
                out_fqs[decision].write('{:s}\n'.format(str(read)))
                filter_counts[decision] += 1
            logger.info('INFO: Wrote batch to Fastq files')
            logger.info('INFO: Processed {} reads'.format(sum(filter_counts.values())))
            logger.info('INFO: Of these:')
            logger.info('INFO:   {} passed the threshold'.format(filter_counts[True]))
            logger.info('INFO:   {} failed the threshold'.format(filter_counts[False]))

            
    in_fq.close()
    for fq in out_fqs.values():
        fq.close()
        
        
@click.group()
def cli():
    pass


@cli.command()
@click_log.simple_verbosity_option(logger)
@click.option('-p', '--pos-bam-file', required=True, type=str)
@click.option('-n', '--neg-bam-file', required=True, type=str)
@click.option('-o', '--output-bam-basename', required=True, type=str)
@click.option('-t', '--threshold', required=True, type=float)
@click.option('--keep-unmapped/--no-keep-unmapped', default=False)
def rethreshold_bam(pos_bam_file, neg_bam_file,
                    output_bam_basename,
                    threshold, keep_unmapped):
    logger.info('INFO: Rethresholding bams at {}'.format(threshold))
    input_bams = [
        pysam.AlignmentFile(bam_fn, mode='rb')
        for bam_fn in (pos_bam_file, neg_bam_file)
    ]
    pass_output_bam_fn, fail_output_bam_fn = get_output_filenames(
        output_bam_basename, 'bam'
    )
    out_bams = {
        True: pysam.AlignmentFile(pass_output_bam_fn,
                                  mode='wb', template=input_bams[0]),
        False: pysam.AlignmentFile(fail_output_bam_fn,
                                   mode='wb', template=input_bams[1])
    }
    filter_counts = {True: 0, False: 0}
    for in_bam in input_bams:
        for record in in_bam.fetch(until_eof=keep_unmapped):
            decision = record.get_tag('ac') > threshold
            out_bams[decision].write(record)
            filter_counts[decision] += 1
    logger.info('INFO: Processed {} reads'.format(sum(filter_counts.values())))
    logger.info('INFO: Of these:')
    logger.info('INFO:   {} passed the threshold'.format(filter_counts[True]))
    logger.info('INFO:   {} failed the threshold'.format(filter_counts[False]))

    for bam in [*input_bams, *out_bams.values()]:
        bam.close()


@cli.command()
@click_log.simple_verbosity_option(logger)
@click.option('-b', '--bam-file', required=True, type=str)
@click.option('-s', '--summary-files-basedir', required=True, type=str)
@click.option('-f', '--fast5-files-basedir', required=True, type=str)
@click.option('-o', '--output-bam-basename', required=True, type=str)
@click.option('-m', '--model-file', required=False, default=None)
@click.option('-k', '--rna-kit', required=False, type=click.Choice(['SQK-RNA001', 'SQK-RNA002']))
@click.option('-t', '--threshold', required=False, default=0.5, type=float)
@click.option('-n', '--batch-size', required=False, default=5000, type=int)
@click.option('-l', '--squiggle-size', required=False, default=2000, type=int)
@click.option('--keep-unmapped/--no-keep-unmapped', default=False)
@click.option('-p', '--processes', required=False, default=-1, type=int)
def filter_bam(bam_file,
               summary_files_basedir,
               fast5_files_basedir,
               output_bam_basename,
               model_file,
               rna_kit,
               threshold,
               batch_size,
               squiggle_size,
               keep_unmapped,
               processes):
    if processes == -1:
        processes = cpu_count()
    logger.info('INFO: {} cpus available'.format(processes))
    read_id_f5_filemap = get_fast5_read_id_mapping(
        summary_files_basedir,
        fast5_files_basedir,
    )
    logger.info(
        ('INFO: Found {} read_ids with associated Fast5'
        ' files in sequencing summaries').format(len(read_id_f5_filemap)))
    pass_output_bam_fn, fail_output_bam_fn = get_output_filenames(
        output_bam_basename
    )
    logger.info('INFO: Reads passing threshold will be written to {}'.format(
        pass_output_bam_fn))
    logger.info('INFO: Reads failing threshold will be written to {}'.format(
        fail_output_bam_fn))
    logger.info('INFO: Loading model')
    model = load_adapter_model(model_file, rna_kit=rna_kit)
    logger.debug('DEBUG: Model Summary:')
    model.summary(print_fn=logger.debug)
    bam_filter(
        bam_file, read_id_f5_filemap, model,
        pass_output_bam_fn, fail_output_bam_fn,
        threshold, batch_size, squiggle_size,
        keep_unmapped, processes
    )


@cli.command()
@click_log.simple_verbosity_option(logger)
@click.option('-q', '--fastq-file', required=True, type=str)
@click.option('-s', '--summary-files-basedir', required=True, type=str)
@click.option('-f', '--fast5-files-basedir', required=True, type=str)
@click.option('-o', '--output-fastq-basename', required=True, type=str)
@click.option('-m', '--model-file', required=False, default=None)
@click.option('-t', '--threshold', required=False, default=0.5, type=float)
@click.option('-n', '--batch-size', required=False, default=5000, type=int)
@click.option('-l', '--squiggle-size', required=False, default=3000, type=int)
@click.option('-c', '--trim/--no-trim', required=False, default=False)
@click.option('-p', '--processes', required=False, default=-1, type=int)
def filter_fastq(fastq_file,
                 summary_files_basedir,
                 fast5_files_basedir,
                 output_fastq_basename,
                 model_file,
                 threshold,
                 batch_size,
                 squiggle_size,
                 trim,
                 processes):
    if processes == -1:
        processes = cpu_count()
    logger.info('INFO: {} cpus available'.format(processes))
    read_id_f5_filemap = get_fast5_read_id_mapping(
        summary_files_basedir,
        fast5_files_basedir,
    )
    logger.info(
        ('INFO: Found {} read_ids with associated Fast5'
        ' files in sequencing summaries').format(len(read_id_f5_filemap)))
    pass_output_fq_fn, fail_output_fq_fn = get_output_filenames(
        output_fastq_basename, 'fastq'
    )
    logger.info('INFO: Reads passing threshold will be written to {}'.format(
        pass_output_fq_fn))
    logger.info('INFO: Reads failing threshold will be written to {}'.format(
        fail_output_fq_fn))
    logger.info('INFO: Loading model')
    model = load_adapter_model(model_file)
    logger.debug('DEBUG: Model Summary:')
    model.summary(print_fn=logger.debug)
    fastq_filter(
        fastq_file, read_id_f5_filemap, model,
        pass_output_fq_fn, fail_output_fq_fn,
        threshold, batch_size, squiggle_size,
        trim, processes
    )
    

if __name__ == '__main__':
    run_bam_filter()