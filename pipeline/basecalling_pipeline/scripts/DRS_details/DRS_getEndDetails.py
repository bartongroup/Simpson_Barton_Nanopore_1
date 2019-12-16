#!/usr/bin/env python
'''
----------------------------------
ONTdrstools.DRS_getEndAlignment.py
----------------------------------

This script parses an aligned set of ONT DRS data, generating several sets of
files focussed on the 5' and 3' ends of the reads. Specifically, it can output
strand specific bigwig (or wig) tracks for the 3' & 5' ends of the reads, and
also write fasta format files containing the 3' and 5' soft-clipped sequences
of the reads. These may contain valuable information such as the poly-A tail 
adapter sequence (likely to be poorly called) and possibly the 5' methyl cap. 

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.1
:created_on: 2018-02-27

Command-line Arguments
======================

**usage\:** 
    DRS_summary.py
    :param: <input bam file>
    :option:`-l|--log` *<file>*
    [:option:`-p|--prefix` <str>]
    [:option:`-v|--verbose`] 
    [:option:`--version`] 
    [:option:`--help`]

Required Parameters
-------------------

:para: <input bam file>

  The input bam file

:option:`--logfile|-l`        

  The name (inc. path) of the log file from the wrapper.

Optional Parameter
------------------

:option: `--prefix|-p` <str> (default: None)

  Prefix string for output filenames. Can optionally include a full path. 
  Defaults to the input filename.

:option: `--fpwig|-f`

  Output 5' ends wig file? (default: False)

:option: `--tpwig|-t` (default: False)

  Output 3' ends wig file?

:option: `--bigwig|-b` (default: False)

  Output wig files as bigwigs?

:option: `--fpclipseq|-g` (default: False)

  Output 5' ends soft-clipped sequence?

:option: `--tpclipseq|-y` (default: False)

  Output 3' ends soft-clipped sequence?

:option:`--help|-h`

  Print a basic description of the tool and its options to STDOUT.

:option:`--version`    

  Show program's version number and exit.
    
:option:`--verbose|-v`     

  Turn on verbose logging (recommended).

Output
======

Undefined, as yet :D
'''

ver=1.30

__scriptname__= "DRS_summary"
__version__ = str(ver)
__usage__ = "\n\t%s <input bam file> -l|--logfile [-p|--prefix <str>]\n\t" \
            "[--version][-v|--verbose][--help]"
__progdesc__ = '''
This script parses an aligned set of ONT DRS data, generating several sets of
summary statistics, making several plots, and generating a wig file for the
pileup of 3' ends of the reads..
'''

__progepi__ = '''
--------------------------
ONTdrstools.DRS_summary.py
--------------------------
'''

import os, sys, pysam, numpy 
import script_options.standard_parsers as sp
import script_options.custom_callables as cc
import script_logging.standard_logging as sl
from Bio import SeqIO, Seq, SeqRecord
from Bio.Alphabet import generic_dna
from parsing_routines.wig_tools import wigData, wigTrack

def addScriptOptions(parser, pos_args, kw_args):
    
    """ add script-specific script options """
    
    script_options_group = parser.add_argument_group('Options')
    
    hlpstr = "Output 3' ends soft-clipped sequence? Output sequences " \
             "are written from nearest ase upstream/downstream of the " \
             "aligned sequence outwars. So for example; a forward " \
             "strand alignment, (5')A..T<alignment>G..C(3'), and a " \
             "reverse strand alignment, (3')T..G<alignment>C..A(5'), " \
             "will be written as: (5' file)> T..A & C..A, (3' file)> " \
             "G..C & G..T"
    option_short_name = "y"
    option_name = "tpclipseq"
    options_default = False
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store_true',
                                      dest = option_name,
                                      help = hlpstr,
                                      default=options_default
                                      )
    
    kw_args.append((option_name, option_name, options_default))

    hlpstr = "Output 5' ends soft-clipped sequence? Output sequences " \
             "are written from nearest ase upstream/downstream of the " \
             "aligned sequence outwars. So for example; a forward " \
             "strand alignment, (5')A..T<alignment>G..C(3'), and a " \
             "reverse strand alignment, (3')T..G<alignment>C..A(5'), " \
             "will be written as: (5' file)> T..A & C..A, (3' file)> " \
             "G..C & G..T"
    option_short_name = "g"
    option_name = "fpclipseq"
    options_default = False
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store_true',
                                      dest = option_name,
                                      help = hlpstr,
                                      default=options_default
                                      )
    
    kw_args.append((option_name, option_name, options_default))
    
    hlpstr = "Minimum length of output soft-clipped sequence"
    option_short_name = "L"
    option_name = "minlen"
    options_default = 0
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type = int,
                                      help = hlpstr,
                                      default=options_default
                                      )
    
    kw_args.append((option_name, option_name, options_default))

    hlpstr = "This defines a genomic region that us used to subset " \
             "the reads that are output. Reads that start or end " \
             "within this region will be reported. The format of " \
             "the string is <chr>:<start>:<stop>. Default is None."
    option_short_name = "r"
    option_name = "region"
    options_default = "None"
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type = str,
                                      help = hlpstr,
                                      default=options_default
                                      )
    
    kw_args.append((option_name, option_name, options_default))
    
    hlpstr = "Padding in bp around the genomic region specified. " \
             "Reads will e reported if they start/end within the " \
             "specified region +/- the padding bp. Default is 0."
    option_short_name = "a"
    option_name = "regionpad"
    options_default = 0
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type = int,
                                      help = hlpstr,
                                      default=options_default
                                      )
    
    kw_args.append((option_name, option_name, options_default))
    
    hlpstr = "Output 3' ends wig file?"
    option_short_name = "t"
    option_name = "tpwig"
    options_default = False
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store_true',
                                      dest = option_name,
                                      help = hlpstr,
                                      default=options_default
                                      )
    
    kw_args.append((option_name, option_name, options_default))

    hlpstr = "Output 5' ends wig file?"
    option_short_name = "f"
    option_name = "fpwig"
    options_default = False
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store_true',
                                      dest = option_name,
                                      help = hlpstr,
                                      default=options_default
                                      )
    
    kw_args.append((option_name, option_name, options_default))
    
    hlpstr = "Output wig files as Bigwigs? Add this flag to output as a wig."
    option_short_name = "b"
    option_name = "wig"
    options_default = True
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store_false',
                                      dest = option_name,
                                      help = hlpstr,
                                      default=options_default
                                      )
    
    kw_args.append((option_name, option_name, options_default))
    
    hlpstr = "Prefix string for output filenames. Can optionally include a " \
             "full path. Defaults to the input filename."
    option_short_name = "p"
    option_name = "prefix"
    options_default = ""
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type = str,
                                      help = hlpstr,
                                      )
    kw_args.append((option_name, option_name, options_default))
    
    return(parser, pos_args, kw_args)

def getRegionPadSeqs(thisbam, region=None, minlength=0, region_pad=20, LOG_EVERY_N=1000,
                     logger=None):
    
    """parses a bam file gathering statistics
    
    Returns a dictionary of various gathered statistics."""
    
    def readEndsInRegion(read, region, pad=20):
        if str(read.reference_name)!=region[0]:
            return(False)
        else:
            if read.reference_start > region[1]-pad and read.reference_end < region[2]+pad:
                return("both")
            elif read.reference_start > region[1]-pad and read.reference_start < region[2]+pad:
                if read.is_reverse:
                    return("3p")
                else:
                    return("5p")
            elif read.reference_end > region[1]-pad and read.reference_end < region[2]+pad:
                if read.is_reverse:
                    return("5p")
                else:
                    return("3p")
            else:
                return(False)
    
    if type(thisbam) is str:
        thisbam = pysam.AlignmentFile(thisbam, "rb")
    
    if region is not None:
        region = region.split(":")
        region[0] = str(region[0])
        region[1] = int(region[1])
        region[2] = int(region[2])
        read_iterator = thisbam.fetch(region[0], region[1], region[2])
    else:
        read_iterator = thisbam.fetch()
    
    seq_records_5p=[]
    seq_records_3p=[]
    nlogs=1
    counter=0
    end_match_counter=0
    no_sequence_reads=0
    
    for read in read_iterator:
        if read.query_sequence is not None:
            doRead = True
            getEnd="both"
            
            if region is not None:
                getEnd = readEndsInRegion(read, region, pad=region_pad)
                if getEnd=="3p" or getEnd=="5p" or getEnd=="both":
                    end_match_counter+=1
                else:
                    doRead=False
            
            if doRead:
                fp_cigtup = read.cigartuples[0]
                tp_cigtup = read.cigartuples[-1]
                if read.is_reverse:
                    fp_cigtup = read.cigartuples[-1]
                    tp_cigtup = read.cigartuples[0]

                if fp_cigtup[0]==4 and fp_cigtup[1]>minlength and (getEnd=="both" or getEnd=="5p"):
                    seq = Seq.Seq(read.query_sequence[0:fp_cigtup[1]][::-1], generic_dna)
                    if read.is_reverse:
                        seq = Seq.Seq(read.query_sequence[read.query_length-fp_cigtup[1]:read.query_length], generic_dna)
                    rec = SeqRecord.SeqRecord(seq, id=read.query_name, name=read.query_name)
                    seq_records_5p.append(rec)

                if tp_cigtup[0]==4 and tp_cigtup[1]>minlength and (getEnd=="both" or getEnd=="3p"):
                    seq = Seq.Seq(read.query_sequence[read.query_length-tp_cigtup[1]:read.query_length], generic_dna)
                    if read.is_reverse:
                        seq = Seq.Seq(read.query_sequence[0:tp_cigtup[1]][::-1], generic_dna)
                    rec = SeqRecord.SeqRecord(seq, id=read.query_name, name=read.query_name)
                    seq_records_3p.append(rec)
        else:
            no_sequence_reads+=1
        
        counter+=1
       
        if (counter % LOG_EVERY_N)==0:
            msg="processed {these} reads...".format(these=(nlogs*LOG_EVERY_N))
            if logger is not None:
                logger.info(msg)
            else:
                print(msg)
            nlogs+=1

    msg = []
    msg.append("Processed {counted} reads...".format(counted=counter))
    msg.append("\tSkipped {counted} reads with no sequence(??)...".format(counted=no_sequence_reads))
    if region is not None:
        msg.append("\t\tidentified {counted} reads ending within the region...".format(counted=end_match_counter))
    msg.append("\t\tidentified {counted} reads with 5-prime skips...".format(counted=len(seq_records_5p)))
    msg.append("\t\tidentified {counted} reads with 3-prime skips...".format(counted=len(seq_records_3p)))
    
    if logger is not None:
        logger.info("\n".join(msg))
    else:
        print("\n".join(msg))
    
    
    return(seq_records_5p, seq_records_3p)

def buildStartAndEndWigData(thisbam, LOG_EVERY_N=1000, logger=None):
    
    """parses a bam file for 3' and 5' ends and builds these into wig-track data
    
    Returns a dictionary of various gathered statistics."""
    
    def formatToWig(wigdata):
        
        """ take in the read position data and output in wigTrack format"""
        
        this_wigTracks = {}
        for key in wigdata.keys():
            track = wigTrack()
            track.wigtype = "fixedStep"
            track.chr = key
            track.start = 1
            track.step = 1
            track.position = numpy.arange(len(wigdata[key]))+track.start
            this_wigTracks[key] = track
            this_wigTracks[key].data = wigdata[key]

        this_wigData = wigData()
        this_wigData.tracks = this_wigTracks
        return(this_wigData)
    
    if type(thisbam) is str:
        thisbam = pysam.AlignmentFile(thisbam, "rb")
    
    
    all_wigdata={
        "fwd":{
            "five_prime":{},
            "three_prime":{}
        },
        "rev":{
            "five_prime":{},
            "three_prime":{}            
        }
    }
    
    chromSizes = dict(zip(thisbam.references, thisbam.lengths))
    for key in chromSizes.keys():
        for strand in all_wigdata.keys():
            for end in all_wigdata[strand].keys():
                all_wigdata[strand][end][key] = numpy.zeros(chromSizes[key])

    counter=0
    nlogs=0
    for read in thisbam.fetch():
        if read.is_reverse:
            all_wigdata["rev"]["five_prime"][read.reference_name][read.reference_end-1]+=1
            all_wigdata["rev"]["three_prime"][read.reference_name][read.reference_start]+=1
        else:
            all_wigdata["fwd"]["five_prime"][read.reference_name][read.reference_start]+=1
            all_wigdata["fwd"]["three_prime"][read.reference_name][read.reference_end-1]+=1
        counter+=1
        if (counter % LOG_EVERY_N)==0:
            msg = "processed {these} reads...".format(these=(nlogs*LOG_EVERY_N))
            if logger is not None:
                logger.info(msg)
            else:
                print(msg)
            nlogs+=1
    
    msg = "Processed {counted} reads...".format(counted=counter)
    if logger is not None:
        logger.info(msg)
    else:
        print(msg)
    msg = "Formatting wig tracks..."
    if logger is not None:
        logger.info(msg)
    else:
        print(msg)

    for strand in all_wigdata.keys():
        for end in all_wigdata[strand].keys():
            all_wigdata[strand][end] = formatToWig(all_wigdata[strand][end])
    
    return(all_wigdata, chromSizes)

if __name__ == '__main__':

    # parse command line options
    # Set standard parser
    parser, pos_args, kw_args = sp.standard_parser(__version__,
                                                   prog = __scriptname__, 
                                                   usage = __usage__,
                                                   description = __progdesc__,
                                                   epilog = __progepi__,
                                                   infile = True,
                                                   outfile = False,
                                                   tmpdir = False)
        
    parser, pos_args, kw_args = addScriptOptions(parser, pos_args, kw_args)
    
    args = parser.parse_args()
           
    # setup standard logging information
    script_logger = sl.standard_logger(__version__, sys.argv, args, pos_args, 
                                       kw_args, script_name=__scriptname__)
    
    if args.tpwig or args.fpwig:
        script_logger.info("Parsing read-end data from %s..." % args.infile)
        wigdata, csizes = buildStartAndEndWigData(args.infile, logger=script_logger)
        if args.tpwig:
            for strand in wigdata.keys():
                fileout = "{}{}_three-prime.bigwig".format(args.prefix, strand)
                wigdata[strand]["three_prime"].writeWig(fileout, isBigWig=args.wig,
                                                        fulldata=True, chromSizes=csizes)
        if args.fpwig:
            for strand in wigdata.keys():
                fileout = "{}{}_five-prime.bigwig".format(args.prefix, strand)
                wigdata[strand]["five_prime"].writeWig(fileout, isBigWig=args.wig,
                                                       fulldata=True, chromSizes=csizes)
    
    if args.tpclipseq or args.fpclipseq:
        script_logger.info("Extracting soft-clipped sequences from %s..." % args.infile)
        region = None
        if args.region!="None":
            region=args.region
        fp, tp = getRegionPadSeqs(args.infile, minlength=args.minlen, region=region, 
                                  region_pad=args.regionpad, logger=script_logger)
        if args.tpclipseq:
            fileout = "{}three-prime_softclipped.fa".format(args.prefix)
            with open(fileout, "w") as out_handle:
                SeqIO.write(tp, out_handle, "fasta")
        if args.fpclipseq:
            fileout = "{}five-prime_softclipped.fa".format(args.prefix)
            with open(fileout, "w") as out_handle:
                SeqIO.write(fp, out_handle, "fasta")

    script_logger.info("Finished. Have a nice day! ;)")
