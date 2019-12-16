#!/usr/bin/env python
'''
----------------------------
ONTdrstools.DRS_seqErrors.py
----------------------------

This script parses an aligned set of ONT DRS data, looking at the sequencing
errors in the data. The errors distributions are summarized by type and they
are broken down by the base relative to the proportions of each base in the
reference sequence underlying the read alignment.

Note that this script requires that the cs tags exist for read alignments in
the input bam file. See: https://github.com/lh3/minimap2

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.0
:created_on: 2018-05-09

Command-line Arguments
======================

**usage\:** 
    DRS_seqErrors.py
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

__scriptname__= "DRS_seqErrors"
__version__ = str(ver)
__usage__ = "\n\t%s <input bam file> -l|--logfile [-p|--prefix <str>]\n\t" \
            "[--version][-v|--verbose][--help]"
__progdesc__ = '''
This script parses an aligned set of ONT DRS data, looking at the sequencing
errors in the data. The errors distributions are summarized by type and they
are broken down by the base relative to the proportions of each base in the
reference sequence underlying the read alignment.

Note that this script requires that the cs tags exist for read alignments in
the input bam file. See: https://github.com/lh3/minimap2
'''

__progepi__ = '''
----------------------------
ONTdrstools.DRS_seqErrors.py
----------------------------
'''

import os, sys, pysam, numpy, re, matplotlib
matplotlib.use('Agg')
import script_options.standard_parsers as sp
import script_options.custom_callables as cc
import script_logging.standard_logging as sl
import matplotlib.pyplot as plt
from json_tricks import dump

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
    
    return(parser, pos_args, kw_args)

def getLongestAlignments(bamfile, logger=None):
    
    """Get the best alignment of each read - where best == longest"""
    
    logger.info("Parsing alignments, keeping longest alignments for " \
                "multiple mapping reads...")
    
    best_alns={}
    rej_alns={}
    
    for readaln in bamfile.fetch():
        if readaln.query_name not in best_alns.keys():
            best_alns[readaln.query_name] = readaln
        elif readaln.alen > best_alns[readaln.query_name].alen:
            best_alns[readaln.query_name] = readaln
        else:
            rej_alns[readaln.query_name] = readaln
    
    return(best_alns, rej_alns)

def countBaseInstances(thisstr, updatedic):
    
    """ for a string count the a, t, g,& c's and update the input dictionary """

    bases = ["A","T","G","C"]
    for base in bases:
        updatedic[base]+=thisstr.count(base)

    return(updatedic)

def parseCStag(cstag, readseq, logger=None, debug=False):
    
    """Parses and extracts the information stored in the 'cs' flag of bamfile alignments.
    
    The information we're looking for with this are identity matches, deletions in reads
    relative to the reference, insertions in reads relative to the reference and finally
    substitutions of the reference base for other bases. 
    
    See https://github.com/lh3/minimap2"""
    
    # define the regex that matches the cs flag components
    r = re.compile(":[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+")
    
    # match vs the sequence
    csarr = numpy.array(r.findall(cstag))
    
    # build a dictionary of the results including recording the reference bases and
    # including the by-bp stats.
    
    cs_stats = {"identity":0,
                "insertion":0,
                "deletion":0,
                "substitution":0,
                "refbases":""}
    
    cs_bpstats = {"insertion":{"A":0, "T":0, "C":0, "G":0},
                  "deletion":{"A":0, "T":0, "C":0, "G":0},
                  "substitution":{"A":[], "T":[], "C":[], "G":[]},
                  "identity":{"A":0, "T":0, "C":0, "G":0}}
    
    pos_in_read = 0
    for block in csarr:
        # first parse identities
        if block.startswith(":"): 
            ilen = int(block.split(":")[1])
            # add identity match length
            cs_stats["identity"]+=ilen
            # list sequence matching the block
            bases = readseq[pos_in_read:pos_in_read+ilen]
            # update current position in the read
            pos_in_read+=ilen
            # append current list of reference bases with the matchign bases
            cs_stats["refbases"] = "{}{}".format(cs_stats["refbases"], bases)
            # count the ATCG entries in the sequence and update the bp
            cs_bpstats["identity"] = countBaseInstances(bases, cs_bpstats["identity"])
        # then parse substitutions
        elif block.startswith("*"):
            bases = block.split("*")[1].upper()
            if len(bases)==2:
                # substitutions are recorded individually so increment substitutions
                cs_stats["substitution"]+=1
                # bases[0] is the reference base, bases[1] is the substitution
                cs_bpstats["substitution"][bases[0]].append(bases[1])
                # append current list of reference bases with the reference base
                cs_stats["refbases"] = "{}{}".format(cs_stats["refbases"], bases[0])
                # update current position in the read
                pos_in_read+=1
            else:
                # this should never happen. If it does something is screwed!
                raise ValueError("Error in substitution block: {}".format(block))
        # then parse insertions
        elif block.startswith("+"):
            bases = block.split("+")[1].upper()
            # insertions are recorded in blocks so increment the insertion count with 
            # block length. There are no relevant reference bases to add for insertions
            # so we'ell skip that bit.
            cs_stats["insertion"]+=len(bases)
            # count the ATCG entries in the sequence and update the bp
            cs_bpstats["insertion"] = countBaseInstances(bases, cs_bpstats["insertion"])
            # update current position in the read
            pos_in_read+=len(bases)
        # finally parse deletions
        elif block.startswith("-"):
            bases = block.split("-")[1].upper()
            # deletions are recorded in blocks so increment the deletion count with 
            # block length.
            cs_stats["deletion"]+=len(bases)
            # count the ATCG entries in the sequence and update the bp
            cs_bpstats["deletion"] = countBaseInstances(bases, cs_bpstats["deletion"])
            # append current list of reference bases with the reference base
            cs_stats["refbases"] = "{}{}".format(cs_stats["refbases"], bases)
        else:
            raise ValueError("Malformed CS tag: {}".format(block))
    
    # convert substitutions distionary entries to numpy arrays not lists
    for key in cs_bpstats["substitution"].keys():
        cs_bpstats["substitution"][key] = numpy.array(cs_bpstats["substitution"][key])
        
    # return stats
    return(cs_stats, cs_bpstats)

def getGlobalAlignmentStats(reads, parseCS=True, logger=None):
    
    """Get a summary of the alignment stats for the reads"""
    
    logger.info("Building global alignments error stats...")
    
    stats = {"matches":[],
             "insertion":[],
             "deletion":[],
             "skip":[],
             "softclip":[],
             "hardclip":[],
             "padding":[],
             "seqmatch":[],
             "seqmismatch":[],
             "back":[],
             "EditDist":[],
             "nbases":[],
             "nalignedbases":[]
            }
    
    if parseCS:
        stats["refbases"]={"A":0, "T":0, "C":0, "G":0}
        stats["identity"]=[]
        stats["substitution"]=[]
        stats["bp_stats"]={"insertion":{"A":0, "T":0, "C":0, "G":0},
                           "deletion":{"A":0, "T":0, "C":0, "G":0},
                           "substitution":{"A":[], "T":[], "C":[], "G":[]},
                           "identity":{"A":0, "T":0, "C":0, "G":0}}
    
    for read in reads:
        # basic info
        stats["nbases"].append(read.query_length)
        stats["nalignedbases"].append(read.query_alignment_length)
        
        # sam cigar information
        read_cigar_stats = read.get_cigar_stats()[0]
        stats["matches"].append(read_cigar_stats[0])
        stats["insertion"].append(read_cigar_stats[1])
        stats["deletion"].append(read_cigar_stats[2])
        stats["skip"].append(read_cigar_stats[3])
        stats["softclip"].append(read_cigar_stats[4])
        stats["hardclip"].append(read_cigar_stats[5])
        stats["padding"].append(read_cigar_stats[6])
        stats["seqmatch"].append(read_cigar_stats[7])
        stats["seqmismatch"].append(read_cigar_stats[8])
        stats["back"].append(read_cigar_stats[9])
        stats["EditDist"].append(read_cigar_stats[10])
        
        # additional cs flag information
        if parseCS:
            cs_stats, bp_stats = parseCStag(read.get_tag('cs'), read.seq)
            
            # sanity checks:
            if cs_stats["insertion"]!=read_cigar_stats[1] or cs_stats["deletion"]!=read_cigar_stats[2] or (cs_stats["identity"]+cs_stats["substitution"])!=read_cigar_stats[0]:
                print(read.query_name)
                print("cs stats\n", cs_stats)
                print("cigar stats\n", read_cigar_stats)
                raise ValueError("cs flag information does not tally with sam cigar string information")
            else:
                stats["refbases"] = countBaseInstances(cs_stats["refbases"], stats["refbases"])
                stats["identity"].append(cs_stats["identity"])
                stats["substitution"].append(cs_stats["substitution"])
                for key in bp_stats.keys():
                    for base in bp_stats[key]:
                        if key=="substitution":
                            stats["bp_stats"][key][base] = numpy.append(stats["bp_stats"][key][base],bp_stats[key][base])
                        else:
                            stats["bp_stats"][key][base]+=bp_stats[key][base]
    
    for key in stats.keys():
        if key!="bp_stats" and key!="refbases":
            stats[key] = numpy.array(stats[key])
    
    return(stats)

def plotErrorDistributions(aln_stats, nbins=100, saveas=None, logger=None):
    
    """Plot sthe distributions of each error type in the data."""
    
    fig = plt.figure(figsize=(8,18), dpi=150)
    
    plt.subplot(411)
    x = plt.hist(aln_stats["nalignedbases"]/aln_stats["nbases"], bins=nbins, alpha=0.5, label="aligned")
    plt.xlabel("fraction of bases in read")
    plt.legend(loc=2)
    plt.ylabel("count")
    
    plt.subplot(412)
    y = plt.hist(aln_stats["identity"]/aln_stats["nalignedbases"], bins=nbins, alpha=0.5, label="identity match")
    plt.xlabel("fraction of aligned bases in read")
    plt.legend(loc=2)
    plt.ylabel("count")
    
    plt.subplot(413)
    x = plt.hist(aln_stats["insertion"]/aln_stats["nalignedbases"], bins=nbins, alpha=0.3, label="insertions")
    y = plt.hist(aln_stats["deletion"]/aln_stats["nalignedbases"], bins=x[1], alpha=0.3, label="deletions")
    z = plt.hist(aln_stats["substitution"]/aln_stats["nalignedbases"], bins=x[1], alpha=0.3, label="substitutions")
    plt.xlabel("fraction of aligned bases in read")
    plt.legend(loc=1)
    plt.ylabel("count")
    
    plt.subplot(414)
    y = aln_stats["skip"]/aln_stats["nalignedbases"]
    y = plt.hist(y[numpy.where(y!=0)[0]], bins=x[1], alpha=0.3, label="skips")
    plt.legend(loc=1)
    plt.ylabel("count")
    plt.xlabel("fraction of aligned bases in read")
    
    if saveas is not None:
        if logger is not None:
            logger.info("Saving distribution plots to {}...".format(saveas))
        plt.savefig(saveas, dip=600)
    
    return(fig)

def getProportionStats(aln_stat, logger=None):
    
    """ convert the alignment stats to proportions so we can compare the occurance of
    each error type vs the reference bp proportions """
    
    logger.info("Fractions of each base in the reference sequence underlying each read:")
    proportions={"refbases":{},
                 "bp_stats":{}}
    
    
    
    for base in aln_stats["refbases"]:
        proportion = aln_stats["refbases"][base]/sum(aln_stats["refbases"].values())
        SE = numpy.sqrt((proportion*(1-proportion))/sum(aln_stats["refbases"].values()))
        CI = (1.96*SE) + (0.5/sum(aln_stats["refbases"].values()))
        proportions["refbases"][base] = {"proportion": proportion, "SE": SE, "95CI": CI}
        logger.info("{}: {:.2f} +/-{:.2f} % (95% CI)".format(base, proportion*100, CI*100))
        
    for key in aln_stats["bp_stats"]:
        proportions["bp_stats"][key] = {}
        if key!="substitution":
            logger.info("{} fractions relative to all {}s by (reference) base:".format(key, key))
            for base in aln_stats["bp_stats"][key].keys():
                proportion = aln_stats["bp_stats"][key][base]/aln_stats[key].sum()
                SE = numpy.sqrt((proportion*(1-proportion))/aln_stats[key].sum())
                CI = (1.96*SE) + (0.5/aln_stats[key].sum())
                proportions["bp_stats"][key][base] = {"proportion": proportion, "SE": SE, "95CI": CI}
                logger.info("{}: {:.2f} +/-{:.2f} % (95% CI)".format(base, proportion*100, CI*100))
    
                
    logger.info("\nSubstitution fractions relative to all substitutions by reference base:")
    for base in aln_stats["bp_stats"]["substitution"].keys():
        proportion = len(aln_stats["bp_stats"]["substitution"][base])/aln_stats["substitution"].sum()
        SE = numpy.sqrt((proportion*(1-proportion))/aln_stats["substitution"].sum())
        CI = (1.96*SE) + (0.5/aln_stats["substitution"].sum())
        proportions["bp_stats"]["substitution"][base] = {"proportion": proportion, "SE": SE, "95CI": CI, "breakdown":{}}
        logger.info("{}({:.2f} +/-{:.2f} % 95% CI):".format(base, proportion*100, CI*100))
        baseto_unique, baseto_counts = numpy.unique(aln_stats["bp_stats"]["substitution"][base], return_counts=True)
        baseto_dict = dict(zip(baseto_unique, baseto_counts))
        proportions["bp_stats"]["substitution"][base]["breakdown"] = {}
        logger.info("\tSubstitution fractions relative to all substitutions of reference base {}, by target base:".format(base))
        for baseto in baseto_dict:
            proportion = baseto_dict[baseto]/baseto_counts.sum()
            SE = numpy.sqrt((proportion*(1-proportion))/baseto_counts.sum())
            CI = (1.96*SE) + (0.5/baseto_counts.sum())
            proportions["bp_stats"]["substitution"][base]["breakdown"][baseto] = {"proportion": proportion, "SE": SE, "95CI": CI}
            logger.info("\t{}: {:.2f} +/-{:.2f} % (95% CI)".format(baseto, proportion*100, CI*100))
    
    return(proportions)

def plotProportions(proportions, saveas=None, logger=None):
    
    """plot the error proportions fir each of the categories, and the details of the substitutions"""
    
    bases = ["A","T","G","C"]
    
    fig1 = plt.figure(figsize=(10,6), dpi=150)
    plotprops = []
    ploterrors = []
    for base in bases:
        plotprops.append(proportions['refbases'][base]["proportion"])
        ploterrors.append(proportions['refbases'][base]["95CI"])
    refline = plt.plot(bases, plotprops, linestyle='--', zorder=1)
    refpoints = plt.errorbar(bases, plotprops, ploterrors, fmt="o", markersize=2, label="Reference proportions", zorder=2, color=refline[-1].get_color())
    
    for key in proportions['bp_stats'].keys():
        plotprops = []
        ploterrors = []
        for base in bases:
            plotprops.append(proportions['bp_stats'][key][base]["proportion"])
            ploterrors.append(proportions['bp_stats'][key][base]["95CI"])
        thisline = plt.plot(bases, plotprops, linestyle='--', zorder=1)
        thispoints = plt.errorbar(bases, plotprops, ploterrors, fmt="o", markersize=2, label=key, zorder=2, color=thisline[-1].get_color())
    
    plt.legend(loc=3)
    plt.xlabel("Base Pair")
    plt.ylabel("Proportion of base in the set")
    
    if saveas is not None:
        if logger is not None:
            logger.info("Saving error proportion plots to {}_errprops.png ...".format(saveas))
        plt.savefig("{}_errprops.png".format(saveas), dip=600)
    
    fig2 = plt.figure(figsize=(10,6), dpi=150)
    x=141
    p=None
    coldict={"A":"red", "T":"green", "G":"blue", "C":"purple"}
    for base in bases:
        if p is None:
            p = plt.subplot(x)
            ax = plt.gca()
            plt.ylabel("Proportion of target base in substitutions from the reference base")
        else:
            ax = plt.subplot(x, sharey=p)
            plt.setp(ax.get_yticklabels(), visible=False)
        plotbases = []
        plotprops = []
        ploterrors = []
        for baseto in bases:
            if baseto in proportions['bp_stats']["substitution"][base]["breakdown"].keys():
                plotbases.append("{}->{}".format(base, baseto))
                plotprops.append(proportions['bp_stats']["substitution"][base]["breakdown"][baseto]["proportion"])
                ploterrors.append(proportions['bp_stats']["substitution"][base]["breakdown"][baseto]["95CI"])
        thisline = plt.plot(plotbases, plotprops, linestyle='--', zorder=1, color=coldict[base])
        thispoints = plt.errorbar(plotbases, plotprops, ploterrors, fmt="o", markersize=2, label="{} substitutions".format(base), zorder=2, color=thisline[-1].get_color())
        plt.ylim((0,1.0))
        x+=1
        plt.legend(loc=1)
    
    if saveas is not None:
        if logger is not None:
            logger.info("Saving substitution breakdown proportion plots to {}_subsprops.png ...".format(saveas))
        plt.savefig("{}_subsprops.png".format(saveas), dip=600)
    
    return(fig1,fig2)

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
    
    # point at the bamfile
    script_logger.info("loading bamfile....")
    thisbam = pysam.AlignmentFile(args.infile, "rb")
    
    if thisbam.mapped>0:
        # parse the alignments for the stats
        alns = getLongestAlignments(thisbam, logger=script_logger)
        aln_stats = getGlobalAlignmentStats(alns[0].values(), logger=script_logger)
        fig1 = plotErrorDistributions(aln_stats, saveas="{}_errDists.png".format(args.prefix), logger=script_logger)
        
        # log the output accuracy
        mean_identity = (aln_stats["identity"]/aln_stats["nalignedbases"]).mean()
        stddev_identity = (aln_stats["identity"]/aln_stats["nalignedbases"]).std()
        script_logger.info("Match percentage for this data (+-std dev): " \
                           "{:.2f}+-{:.3f}".format(mean_identity, stddev_identity))
        
        # convert the alignment stats to proportion stats
        proportions = getProportionStats(aln_stats, logger=script_logger)
        
        # plot the proportion breakdowns
        fig3, fig4 = plotProportions(proportions, saveas=args.prefix, logger=script_logger)
        
        # write the raw info to an output json file
        fh = open("{}_alignment_stats.json".format(args.prefix), "w")
        #fh.write(json.dumps
        dump({"alignment_stats":aln_stats, "proportion_stats":proportions}, fh,
             sort_keys=True, indent=4, separators=(',', ': '))
        #fh.close()
    else:
        raise ValueError("Bam file {} has no mapped reads. Are you sure there were spike-ins" \
                         "in this sample?".format(args.infile))
    
    script_logger.info("Finished. Have a nice day! ;)")
