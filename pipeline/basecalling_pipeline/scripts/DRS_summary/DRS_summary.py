#!/usr/bin/python
'''
--------------------------
ONTdrstools.DRS_summary.py
--------------------------

This script parses an aligned set of ONT DRS data, generating several sets of
summary statistics, making several plots, and generating a wig file for the
pileup of 3' ends of the reads.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.0
:created_on: 2017-12-18

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

:option: `--prefix|-p` <str>

  Prefix string for output filenames. Can optionally include a full path. 
  Defaults to the input filename.

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

import sys, pysam, numpy, json, matplotlib
matplotlib.use('Agg')
from json_tricks import dump
import script_options.standard_parsers as sp
import script_options.custom_callables as cc
import script_logging.standard_logging as sl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
from parsing_routines.gff_gtf_tools import annotation

def addScriptOptions(parser, pos_args, kw_args):
    
    """ add script-specific script options """
    
    script_options_group = parser.add_argument_group('Options')
    
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
    
    hlpstr = "Annotation (gtf format) filename. Default is None."
    option_short_name = "a"
    option_name = "annotation"
    options_default = ""
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type = cc.input_file,
                                      help = hlpstr,
                                      )
    
    kw_args.append((option_name, option_name, options_default))
    
    return(parser, pos_args, kw_args)

def parseBam(thisbam, region=None, LOG_EVERY_N=10000):
    
    """parses a bam file gathering statistics
    
    Returns a dictionary of various gathered statistics."""
    
    if type(thisbam) is str:
        thisbam = pysam.AlignmentFile(args.infile, "rb")
    
    refs = thisbam.references
    read_lengths={}
    read_lens=[]
    read_aligned_lengths={}
    aligned_len=[]
    read_end_pos={}
    tprime_skipped=[]
    threeprime_skipped={}
    fprime_skipped=[]
    fiveprime_skipped={}
    
    if region is None:
        read_iterator = thisbam.fetch()
    else:
        read_iterator = thisbam.fetch(str(region[0]), int(region[1]), int(region[2]))
    
    counter=0
    nlogs=1
    for read in read_iterator:
        this_read_len = read.query_length
        if this_read_len==0:
            this_read_len = read.infer_query_length()
        read_lens.append(this_read_len)
        try: 
            read_lengths[read.query_length].append(read.query_name)
        except KeyError:
            read_lengths[read.query_length]=[read.query_name]
        
        aligned_len.append(read.query_alignment_length)
        try: 
            read_aligned_lengths[read.query_alignment_length].append(read.query_name)
        except KeyError:
            read_aligned_lengths[read.query_alignment_length]=[read.query_name]
        
        read_ref = refs[read.reference_id]
        if read_ref not in read_end_pos.keys():
            read_end_pos[read_ref]={}
        
        tp_skipped = 0
        fp_skipped = 0
        if read.is_reverse:
            read_end = read.reference_start
            if read.cigar[0][0]==4:
                tp_skipped = read.cigar[0][1]
            if read.cigar[-1][0]==4:
                fp_skipped = read.cigar[-1][1]  
        else:
            read_end = read.reference_end
            if read.cigar[-1][0]==4:
                tp_skipped = read.cigar[-1][1]
            if read.cigar[0][0]==4:
                fp_skipped = read.cigar[0][1]
        
        try:
            read_end_pos[read_ref][read_end].append(read.query_name)
        except KeyError:
            read_end_pos[read_ref][read_end]=[read.query_name]
        
        tprime_skipped.append(tp_skipped)
        try:
            threeprime_skipped[tp_skipped].append(read.query_name)
        except KeyError:
            threeprime_skipped[tp_skipped]=[read.query_name]
        
        fprime_skipped.append(fp_skipped)
        try:
            fiveprime_skipped[fp_skipped].append(read.query_name)
        except KeyError:
            fiveprime_skipped[fp_skipped]=[read.query_name]
        
        counter+=1
        if (counter % LOG_EVERY_N)==0:
            print "processed %i reads..." % (nlogs*LOG_EVERY_N)
            nlogs+=1
    
    print "Finished. Processed %i reads." % counter
        
    return({"read_lengths": numpy.array(read_lens),
            "read_lengths_detail": read_lengths,
            "aligned_lengths": numpy.array(aligned_len),
            "aligned_lengths_detail": read_aligned_lengths,
            "3p_ends": read_end_pos,
            "3p_end_skips": numpy.array(tprime_skipped),
            "3p_end_skips_detail":threeprime_skipped,
            "5p_end_skips": numpy.array(fprime_skipped),
            "5p_end_skips_details": fiveprime_skipped})

def plotReadLengthHist(stats_data, binwidth=20, alpha=0.5, figsize=(12,8), fileout=None):
    
    """ plot the read length histograms """
        
    fig = plt.figure(figsize=figsize)
    
    common_params = dict(
                        bins=numpy.arange(0, max(stats_data["aligned_lengths"])+binwidth, binwidth), 
                        log=True,
                        alpha=alpha
                        )
    
    plt.hist(stats_data["read_lengths"],
             color="darkblue", label="Read length", **common_params)
    plt.hist(stats_data["aligned_lengths"],
             color="darkgreen", label="aligned length", **common_params)
    plt.legend(loc=1)
    plt.xlabel("length (bp)")
    plt.ylabel("number of reads")
    if fileout is not None:
        plt.savefig(fileout, dpi=300, transparent=True)
    
    return(fig)

def plotLengthDensity(stats_data, gridsize=100, cmap="Greens", figsize=(16,12), fileout=None):
    
    """ plot the read length histograms """
        
    fig = plt.figure(figsize=figsize)
    
    plt.hexbin(stats_data["aligned_lengths"], stats_data["read_lengths"],
               gridsize=gridsize, xscale="log", yscale="log", bins="log",
               cmap=cmap)
    
    axes = plt.gca()
    mins = [axes.get_xlim()[0], axes.get_ylim()[0]]
    maxs = [axes.get_xlim()[1], axes.get_ylim()[1]]
    plt.plot([min(mins),min(maxs)],[min(mins),min(maxs)], color="black")
    
    cbar = plt.colorbar(shrink=0.4)
    cbar.set_label(r'$log_{10}(no. of reads)$', rotation=270, labelpad=15)
    plt.xlabel("aligned length (bp)")
    plt.ylabel("read length (bp)")
    if fileout is not None:
        plt.savefig(fileout, dpi=300, transparent=True)
    
    return(fig)

def plotSkipLengthHist(stats_data, binwidth=20, alpha=1.0, figsize=(12,8), fileout=None):
    
    """ plot the read length histograms """
        
    fig = plt.figure(figsize=(12,8))
    
    common_params = dict(
                        bins=numpy.arange(0, max(stats_data["3p_end_skips"])+binwidth, binwidth), 
                        log=True,
                        alpha=alpha
                        )
    
    plt.subplot(2,1,1)
    plt.hist(stats_data["3p_end_skips"],
             color="darkblue", label="3' end skipped bases", **common_params)
    plt.ylabel("number of reads")
    plt.legend(loc=1)
    plt.subplot(2,1,2)
    plt.hist(stats_data["5p_end_skips"],
             color="darkgreen", label="5' end skipped bases", **common_params)
    plt.ylabel("number of reads")
    plt.xlabel("read length (bp)")
    plt.legend(loc=1)
    if fileout is not None:
        plt.savefig(fileout, dpi=300, transparent=True)
    
    return(fig)

def annotCount(annot_filename, bam_file, feature="genes", annot_fmt="gtf", LOG_EVERY_N=10000):
    
    """ count the number of reads mapping to each feature """
    
    # load annotation
    print "loading annotation..."
    annot = annotation(annot_filename, filetype=annot_fmt)
    annot.set_feature(feature)
    afeats = annot.get_selection()
    
    print "counting reads..."
    # pointer to the bam file (again)
    if type(bam_file) is str:
        thisbam = pysam.AlignmentFile(bam_file, "rb")
    else:
        thisbam=bam_file
    
    def getRevReads(bamfile, chrid, start, stop):
        
        """ Make a generator object that returns only rev strand reads"""
        
        return(read for read in bamfile.fetch(chrid, start, stop) if read.is_reverse)
    
    def getFwdReads(bamfile, chrid, start, stop):
        
        """ Make a generator object that returns only rev strand reads"""
        
        return(read for read in bamfile.fetch(chrid, start, stop) if not read.is_reverse)
    
    counts = numpy.zeros(len(afeats),
                         dtype=[("name","|S20"),("count","int"),("frac_read_coverage","float")])
    
    i=0
    nlogs=1
    for afeat in afeats:
        if afeat.strand=="+":
            reads = getFwdReads(thisbam, afeat.chrid, afeat.start, afeat.stop)
        elif afeat.strand=="-":
            reads = getRevReads(thisbam, afeat.chrid, afeat.start, afeat.stop)
        
        count=0
        bases_covered=0
        for read in reads:
            count+=1
            start = afeat.start
            stop = afeat.stop
            if read.reference_start>start:
                start = read.reference_start
            if read.reference_end<stop:
                stop = read.reference_end
            bases_covered = bases_covered+(stop-start)
        
        frac_coverage = 0
        if count>0:
            mean_coverage = float(bases_covered)/(afeat.stop-afeat.start)
            frac_coverage = float(mean_coverage)/count
                    
        counts["name"][i] = afeat.desc["gene_id"]
        counts["count"][i] = count
        counts["frac_read_coverage"][i] = frac_coverage
                
        i+=1
        if (i % LOG_EVERY_N) == 0:
            print "processed %i features..." % (nlogs*LOG_EVERY_N)
            nlogs+=1
    
    return(counts)

def plotAnnotCountvCov(counts, log=True, bins=50, feature_label="genes", cmap="Greens",
                       figsize=(12,8), fileout=None):
    
    """plots the counts vs the fractional coverage for a dataset
    
    uses output from annotCount. 
    """
    
    fig = plt.figure(figsize=figsize)
    
    ind = numpy.where(counts["count"]>0)[0]
    xscale="linear"
    if log:
        xscale="log"
    
    plt.hexbin(counts["count"][ind], counts["frac_read_coverage"][ind],
               gridsize=bins, bins="log", xscale=xscale, marginals=False, cmap=cmap)
    cbar = plt.colorbar(shrink=0.4)
    cbar.set_label(r'$log_{10}(no. of %s)$' % feature_label, rotation=270, labelpad=15)
     
    plt.ylabel(r"mean fractional read coverage")
    plt.xlabel(r"read counts")
    if fileout is not None:
        plt.savefig(fileout, dpi=300, transparent=True)
    
    return(fig)
     
def plotAnnotCountScatter(counts1, counts2, c1label="", c2label="",
                          log=True, xylog=True, bins=50, feature_label="genes",
                          figsize=(10,8), fileout=None):
    
    """plots the read counts for two datasets for the same annotation
    
    uses output from annotCount. Make sure the same annotation was used
    to generate both sets of counts - that way things will all be in the right order
    """
    
    fig = plt.figure(figsize=figsize)
    xscale="linear"
    yscale="linear"
    x=counts1["count"]
    y=counts2["count"]
    if xylog:
        ind1 = counts1["count"]>0
        ind2 = counts2["count"]>0
        ind=ind1*ind2
        x=counts1["count"][ind]
        y=counts2["count"][ind]        
        xscale="log"
        yscale="log"
            
    plt.hexbin(x, y, gridsize=bins, bins="log", xscale=xscale, 
               yscale=yscale, marginals=False,cmap='Greens')
    
    axes = plt.gca()
    mins = [axes.get_xlim()[0], axes.get_ylim()[0]]
    maxs = [axes.get_xlim()[1], axes.get_ylim()[1]]
    plt.plot([min(mins),min(maxs)*2],[min(mins),min(maxs)*2], color="black")

    cbar = plt.colorbar(shrink=0.4)
    cbar.set_label('# of %s' % feature_label, rotation=270, labelpad=15)
     
    plt.ylabel(r"%s read counts" % c1label)
    plt.xlabel(r"%s read counts" % c2label)
    if fileout is not None:
        plt.savefig(fileout, dpi=300, transparent=True)
    
    return(fig)

def plotAnnotCovScatter(counts1, counts2, c1label="", c2label="",
                        log=True, bins=50, feature_label="genes",
                        figsize=(10,8), fileout=None):
    
    """plots the fractional read coverage for two datasets for the same annotation
    
    uses output from annotCount. Make sure the same annotation was used
    to generate both sets of counts - that way things will all be in the right order
    """
    
    fig = plt.figure(figsize=figsize)
        
    plt.hexbin(counts1["frac_read_coverage"], counts2["frac_read_coverage"],
               gridsize=bins, bins="log", xscale="linear", marginals=False,
               cmap='Greens')
    plt.plot([0,1],[0,1], color="black")
    cbar = plt.colorbar(shrink=0.4)
    cbar.set_label('# of %s' % feature_label, rotation=270, labelpad=15)
     
    plt.ylabel(r"%s mean fractional read coverage" % c1label)
    plt.xlabel(r"%s mean fractional read coverage" % c2label)
    if fileout is not None:
        plt.savefig(fileout, dpi=300, transparent=True)
    
    return(fig)  

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
        
    script_logger.info("Parsing alignment data from %s..." % args.infile)
    
    print "Reading bamfile..."
    thisbam = pysam.AlignmentFile(args.infile, "rb")
    thisbamstats = parseBam(thisbam)
    thiscounts=[]
    
    if len(thisbamstats["read_lengths"])>0:
    
        figs = []
        script_logger.info("Plotting read length histograms...")
        if args.prefix=="":
            fileout = None
        else:
            fileout = "%s_readLenHist.png" % args.prefix
            script_logger.info("\t... writing out to file %s" % fileout)        
        figs.append(plotReadLengthHist(thisbamstats,fileout=fileout))
        
        script_logger.info("Plotting read length vs alignment length...")
        if args.prefix=="":
            fileout = None
        else:
            fileout = "%s_readLen_v_alignedLen.png" % args.prefix
            script_logger.info("\t... writing out to file %s" % fileout)   
        figs.append(plotLengthDensity(thisbamstats,fileout=fileout))
        
        script_logger.info("Plotting 5' and 3' skipped sequence lengths...")
        if args.prefix=="":
            fileout = None
        else:
            fileout = "%s_skipLenHist.png" % args.prefix
            script_logger.info("\t... writing out to file %s" % fileout)
        figs.append(plotSkipLengthHist(thisbamstats,fileout=fileout))
        
        if args.annotation is not None:
            thiscounts = annotCount(args.annotation, thisbam)
            script_logger.info("Plotting gene counts...")
            if args.prefix=="":
                fileout = None
            else:
                fileout = "%s_geneexpression.png" % args.prefix
                script_logger.info("\t... writing out to file %s" % fileout)
            figs.append(plotAnnotCountvCov(thiscounts,fileout=fileout))
    else:
        script_logger.info("No aligned reads found to plot.")
    
    script_logger.info("Writing stats to json file %s ..." % ("%s_plotstats.json" % args.prefix))
    dump({"stats": thisbamstats, "counts": thiscounts}, "%s_plotstats.json" % args.prefix,
          sort_keys=True, indent=4, separators=(',', ': '))
    script_logger.info("Finished. Have a nice day! ;)")