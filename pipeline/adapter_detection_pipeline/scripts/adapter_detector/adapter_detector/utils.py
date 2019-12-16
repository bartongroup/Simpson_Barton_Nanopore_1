import os


def find_summary_files(basedir):
    '''
    Search basedir recursively to find files named
    sequencing_summary.txt, will follow symlinks
    '''
    i = 0
    for currdir, _, files in os.walk(basedir, followlinks=True):
        for fn in files:
            if fn == 'sequencing_summary.txt':
                i += 1
                yield os.path.join(currdir, fn)
    if not i:
        raise IOError('Could not find any sequencing_summary.txt files')


def parse_sequencing_summary(seq_summary_fn):
    '''
    Parse a sequencing summary file to create a
    dictionary of read id to fast5 file name mappings.
    Fast5 file names do not include full path yet.
    '''
    read_id_filemap = {}
    with open(seq_summary_fn) as ss:
        # skip header
        _ = next(ss)
        for record in ss:
            filename, read_id, *_ = record.split('\t')
            read_id_filemap[read_id] = filename
    return read_id_filemap


def find_fast5s(basedir):
    '''
    Search basedir recursively to find fast5 files 
    will follow symlinks
    '''
    i = 0
    for currdir, _, files in os.walk(basedir, followlinks=True):
        for fn in files:
            if fn.endswith('.fast5'):
                i += 1
                yield os.path.join(currdir, fn)
    if not i:
        raise IOError('Could not find any Fast5 files')


def get_full_filepath_for_read_id_mappings(read_id_filemap, fast5_basedir):
    '''
    Extend the filepath for fast5s in a read id to fast5 mapping.
    Files or read ids that cannot be found in fast5_basedir
    will be removed.
    '''
    fullpath_filemap = {}
    read_id_map_with_fullpaths = {}
    for fullpath in find_fast5s(fast5_basedir):
        _, basename = os.path.split(fullpath)
        if basename not in fullpath_filemap:
            fullpath_filemap[basename] = [fullpath,]
        else:
            fullpath_filemap[basename].append(fullpath)
    read_id_map_with_fullpaths = {}
    for read_id, basename in read_id_filemap.items():
        try:
            read_id_map_with_fullpaths[read_id] = fullpath_filemap[basename]
        except KeyError:
            continue
    return read_id_map_with_fullpaths


def get_fast5_read_id_mapping(summaries_basedir, fast5_basedir):
    '''
    Creates a dictionary mapping read_ids to the location of
    their raw fast5 files on disk.
    
    Parameters:
    -----------
    
        summaries_basedir: str, required
          Directory which is recursively searched for sequencing_summary.txt
          files produced by albacore or other basecaller
        fast5_basedir: str, required
          Directory which is recursively searched for FAST5 files.

    Returns:
    --------

        read_id_filemap: dict
          Dict of read_id: fast5 filepath key value pairs
    '''
    read_id_filemap = {}
    for summary_fn in find_summary_files(summaries_basedir):
        read_id_filemap.update(parse_sequencing_summary(summary_fn))
    if not read_id_filemap:
        raise ValueError('Read ID filemap is empty')
    read_id_filemap = get_full_filepath_for_read_id_mappings(
        read_id_filemap, fast5_basedir
    )
    if not read_id_filemap:
        raise ValueError('Read ID filemap is empty after extending paths')
    return read_id_filemap


def get_output_filenames(output_basename, filetype='bam'):
    if filetype == 'bam':
        ext = '.bam'
    elif filetype == 'fastq':
        ext = '.fq'
    else:
        raise TypeError('Not FastQ or BAM')
    if os.path.isdir(output_basename):
        raise IOError('output_bam_basename is a directory')
    else:
        return (output_basename + '.passes' + ext,
                output_basename + '.fails' + ext)