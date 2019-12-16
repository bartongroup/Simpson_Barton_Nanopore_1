#!/usr/bin/env python
'''
--------------------------------------------
ONTdrstools.DRS_splitCanonicalSpliceReads.py
--------------------------------------------

This script parses an aligned set of ONT DRS data and examines the splicing 
characteristics of the dataset. The parses the alignment of each read in the
bam file looking for splices and, for each splice present, classifies the 
splice as canonical or non canonical if it has the GU-AG splicing di-nucleotide
motif at the terminals of the intron.

Optionally, it will also examine an annotation gtf/gff and will use this to 
classify the splices as annotated, or not and record which transcripts the splice
is found in (-a). Importantly, if splice sites are novel (i.e., no annotated), the
script examines a padded region (10bp by default) around the splice position for 
alternative annotated splice sites and novel alternative canonical splice sites.

Optionally it will also use the position weight matricies from Sheth et al 2006
(DOI:10.1093/nar/gkl556) to classify splice sites as either U2 or U2 splice sites
if the species designation matches one of the wive spiecies they have position
weight matrices for (human, mouse, fly, c elegans and arabidopsis).

Optionally, the reads can be read out into separate man files for the different
categories.

All the read, splice and transcript classification information is output as json 
files, reads are split into new bam files, and several summary pltos are made.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.3
:created_on: 2018-03-23

Command-line Arguments
======================

**usage\:** 
    DRS_splitCanonicalSpliceReads.py
    :param: <input bam file>
    :option:`-l|--log` *<file>*
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

:option:`--help|-h`

  Print a basic description of the tool and its options to STDOUT.

:option:`--version`    

  Show program's version number and exit.
    
:option:`--verbose|-v`     

  Turn on verbose logging (recommended).

Output
======
'''

ver=1.2

__scriptname__= "DRS_splitCanonicalSpliceReads"
__version__ = str(ver)
__usage__ = "\n\t%s <input bam file> <input genome file> -l|--logfile" \
            "[--version][-v|--verbose][--help]"
__progdesc__ = '''
This script parses an aligned set of ONT DRS data and examines the splicing 
characteristics of the dataset. The parses the alignment of each read in the
bam file looking for splices and, for each splice present, classifies the 
splice as canonical or non canonical if it has the GU-AG splicing di-nucleotide
motif at the terminals of the intron.

Optionally, it will also examine an annotation gtf/gff and will use this to 
classify the splices as annotated, or not and record which transcripts the splice
is found in (-a). Importantly, if splice sites are novel (i.e., no annotated), the
script examines a padded region (10bp by default) around the splice position for 
alternative canonical splice sites, classifying them as annotated or not.

Optionally it will also use the position weight matricies from Sheth et al 2006
(DOI:10.1093/nar/gkl556) to classify splice sites as either U2 or U2 splice sites
if the species designation matches one of the five species they have position
weight matrices for (human, mouse, fly, c elegans and arabidopsis).

Optionally, the reads can be read out into separate man files for the different
categories.

All the read, splice and transcript classification information is output as json 
files, reads are split into new bam files, and several summary pltos are made.

'''

__progepi__ = '''
--------------------------------
DRS_splitCanonicalSpliceReads.py
--------------------------------
'''

import os, sys, pysam, json, re, math, glob, matplotlib, numpy, itertools, time
matplotlib.use('Agg')
import script_options.standard_parsers as sp
import script_options.custom_callables as cc
import script_logging.standard_logging as sl
from Bio import SeqIO, motifs
from Bio.Seq import Seq
from Bio.motifs import matrix
from Bio.Alphabet.IUPAC import unambiguous_dna
from parsing_routines.gff_gtf_tools import annotation
from parsing_routines.general_classes_and_functions import computeIntrons
from bisect import bisect_left
import matplotlib.pyplot as plt

def addScriptOptions(parser, pos_args, kw_args):
    
    """ add script-specific script options """
    
    script_req_group = sp.get_std_req_group(parser)

    hlpstr = "Input bamfile"
    option_short_name = "b"
    option_name = "bamfile"
    
    script_req_group.add_argument('-%s' % option_short_name,
                                  '--%s' % option_name,
                                  action = 'store',
                                  dest = option_name,
                                  type = cc.input_file,
                                  help = hlpstr)
    kw_args.append((option_name, option_name, None))
    
    hlpstr = "Input genome fasta"
    option_short_name = "g"
    option_name = "genomefile"
    
    script_req_group.add_argument('-%s' % option_short_name,
                                  '--%s' % option_name,
                                  action = 'store',
                                  dest = option_name,
                                  type = cc.input_file,
                                  help = hlpstr)
    kw_args.append((option_name, option_name, None))

    
    script_options_group = parser.add_argument_group('Options')
    
    hlpstr = "Path to gtf or gff annotation. If provided, novel splices " \
             "will be identified by comparing the detected intron " \
             "coordinates with the give annotation. Default is None"
    option_short_name = "a"
    option_name = "annotation"
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store',
                                      type = cc.input_file,
                                      dest = option_name,
                                      help = hlpstr
                                      )
    kw_args.append((option_name, option_name, None))
    
    hlpstr = "Annotation format. Options are 'gff3' or 'gft'. Default is 'gff3'."
    option_short_name = ""
    option_name = "input_format"
    option_default = "gff3"
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type = str,
                                      help = hlpstr,
                                      default = option_default,
                                      choices=["gff3","gtf"]
                                      )
    kw_args.append((option_name, option_name, option_default))
    
    hlpstr = "Strip 'Chr' from the beginning of annotation reference chromosomes?"
    option_short_name = ""
    option_name = "stripchr"
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store_true',
                                      dest = option_name,
                                      help = hlpstr
                                      )
    kw_args.append((option_name, option_name, False))
    
    hlpstr = "Comma-separated list of chromosome synonyms, delimited by ':'. " \
             "For example: C:Pt,M:Mt"
    option_short_name = ""
    option_name = "chr_synonyms"
    option_default = ""
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type = str,
                                      help = hlpstr,
                                      default = option_default
                                      )
    kw_args.append((option_name, option_name, option_default))
    
    hlpstr = "Use position weight matrices (default is false)"
    option_short_name = ""
    option_name = "pwm"
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store_true',
                                      dest = option_name,
                                      help = hlpstr
                                      )
    kw_args.append((option_name, option_name, False))
    
    hlpstr = "pwm species to use. Options are: A_thaliana, D_melanogaster, " \
             "M_musculus, C_elegans, or H_sapiens. Default is A_thaliana."
    option_short_name = ""
    option_name = "pwm_species"
    option_default = "A_thaliana"
    option_choices = ["A_thaliana","C_elegans","D_melagonaster","M_musculus","H_sapiens"]
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type=str,
                                      help = hlpstr,
                                      default = option_default,
                                      choices = option_choices
                                      )
    kw_args.append((option_name, option_name, option_default))
    
    hlpstr = "pwm log-odds threshold to use. Default is 3.0"
    option_short_name = ""
    option_name = "pwm_thresh"
    option_default = 3.0
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      help = hlpstr,
                                      default = option_default
                                      )
    kw_args.append((option_name, option_name, option_default))
    
    hlpstr = "Path to pwm files. Default is ./position_weight_matricies"
    option_short_name = ""
    option_name = "pwm_path"
    option_default = os.path.join(os.path.dirname(os.path.realpath(__file__)),"position_weight_matricies")
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      help = hlpstr,
                                      type = cc.input_path,
                                      default = option_default
                                      )
    kw_args.append((option_name, option_name, option_default))
    
    hlpstr = "Alternative splicing search region size (+-bp). Default is 10"
    option_short_name = ""
    option_name = "altsplicepad"
    option_default = 10
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      help = hlpstr,
                                      default = option_default
                                      )
    kw_args.append((option_name, option_name, option_default))
    
    hlpstr = "Minimum skips in bam file to consider the skips and intron (bp)." \
             "Default is 20"
    option_short_name = ""
    option_name = "minintronsize"
    option_default = 20
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      help = hlpstr,
                                      default = option_default
                                      )
    kw_args.append((option_name, option_name, option_default))
    
    hlpstr = "Prefix for the output files. Default is the prefex to the input bam file"
    option_short_name = "p"
    option_name = "prefix"
    option_default = ""
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store',
                                      type = str,
                                      dest = option_name,
                                      help = hlpstr,
                                      default = option_default
                                      )
    kw_args.append((option_name, option_name, option_default))
    
    hlpstr = "Split reads into new bam files?"
    option_short_name = ""
    option_name = "splitreads"
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store_true',
                                      dest = option_name,
                                      help = hlpstr
                                      )
    kw_args.append((option_name, option_name, False))
    
    hlpstr = "Which classifications to subset reads by (comma-separated values). Options are " \
             "one or more of: annotated, canonical, U2. Default is annotated."
    option_short_name = ""
    option_name = "spliton"
    option_default = "annotated"
    
    script_options_group.add_argument('--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type=str,
                                      help = hlpstr,
                                      default = option_default
                                      )
    kw_args.append((option_name, option_name, option_default))
    
    return(parser, pos_args, kw_args)
    
def classifySplice(pwms, chrid, start, stop, strand, genome, donorpad=(3,10), acceptorpad=(14,3), minthresh=3.0,
                   species_prefix=""):
    
    donor_splice_site=None
    acceptor_splice_site=None
    if strand=="-":
        donor_splice_site = genome[chrid][stop-donorpad[1]:stop+donorpad[0]].reverse_complement()
        acceptor_splice_site = genome[chrid][start-acceptorpad[1]:start+acceptorpad[0]].reverse_complement()
    else:
        donor_splice_site = genome[chrid][start-donorpad[0]:start+donorpad[1]]
        acceptor_splice_site = genome[chrid][stop-acceptorpad[0]:stop+acceptorpad[1]]

    pwm_scores={"acceptor":{}, "donor":{}}
    
    for pwm in pwms.keys():
        thispwm = pwms[pwm]
        pssm = thispwm.log_odds()
        
        spliceend=None
        if "acceptor" in pwm:
            splice_site = Seq(str(acceptor_splice_site.seq), pssm.alphabet)
            spliceend = "acceptor"
        elif "donor" in pwm:
            splice_site = Seq(str(donor_splice_site.seq), pssm.alphabet)
            spliceend = "donor"
        
        splicetype = re.sub("_{}".format(spliceend),"",pwm)
        
        for position, score in pssm.search(splice_site, threshold=minthresh):
            try:
                pwm_scores[spliceend][splicetype].append((position, score))
            except KeyError:
                pwm_scores[spliceend][splicetype]=[(position, score)]
        
    matching_end_classifications=None
    for key in pwm_scores["acceptor"].keys():
        if key in pwm_scores["donor"]:
            thistuple = (key, pwm_scores["donor"][key], pwm_scores["acceptor"][key])
            if matching_end_classifications is None:
                matching_end_classifications = thistuple
            else:
                if (key.startswith("U12") and matching_end_classifications[0].startswith("U12")) or (key.startswith("U2") and matching_end_classifications[0].startswith("U2")):
                    ascore=None
                    dscore=None
                    for pos, score in thistuple[1]:
                        if ascore is None or score > ascore:
                                ascore = 2**score
                    for pos, score in thistuple[2]:
                        if dscore is None or score > dscore:
                                dscore = 2**score
                    thisscore = ascore+dscore
                    prevscore = 2**(matching_end_classifications[1][0][1])+2**(matching_end_classifications[2][0][1])
                    if prevscore<thisscore:
                        matching_end_classifications = thistuple
                else:
                    U12score = None
                    U2score=None
                    if key.startswith("U12"):
                        for pos, score in thistuple[1]:
                            if U12score is None or score > U12score:
                                U12score = 2**score
                        for pos, score in matching_end_classifications[1]:
                            if U2score is None or score > U2score:
                                U2score = 2**score                                
                        if U12score>25*U2score:
                            matching_end_classifications = thistuple
                    elif matching_end_classifications[0].startswith("U12"):
                        for pos, score in matching_end_classifications[1]:
                            if U12score is None or score > U12score:
                                U12score = 2**score
                        for pos, score in thistuple[1]:
                            if U2score is None or score > U2score:
                                U2score = 2**score
                        if U12score<25*U2score:
                            matching_end_classifications = thistuple
    
    return(matching_end_classifications, pwm_scores)

def charcterizeIntrons(read_iterator, genome, annotation_splices=None, splicepad=0,
                       min_intron_length=20, pwms=None, pwmscorethreshold=3.0, logger=None, LOG_EVERY_N=10000,
                       pwm_species=""):
    
    """Return a dictionary {readID:[(start, stop)...]}
       Listing the intronic sites for each read (identified by 'N' in the cigar strings).
       
       read_iterator can be the result of a .fetch(...) call.
       Or it can be a generator filtering such reads. Example
       samfile.find_introns((read for read in samfile.fetch(...) if read.is_reverse)
    """
    
    def locNearestCanonicals(donor, acceptor, padsize, refname, start, stop):
        
        """ given donor and acceptor regions, locate the nearest canonical splice and 
        return the details """
    
        import re, numpy, itertools
        
        # we're going to find all instances of GA in the donor region, and AG int he acceptor region and 
        # then work out the positions relative to the actual mapped splice
        donor_canonical_matches = numpy.array([m.start() for m in re.finditer('(?=GT)', str(donor))])-padsize
        acceptor_canonical_matches = numpy.array([m.start() for m in re.finditer('(?=AG)', str(acceptor))])-padsize
        possible_canonical_pairs = list(itertools.product(donor_canonical_matches, acceptor_canonical_matches))
        possible_canonical_splices = ["{}:{}-{}".format(refname, x[0]+start, x[1]+stop) for x in possible_canonical_pairs if x!=(0,0)]
        return(possible_canonical_splices)
    
    def locNearestAnnotated(keystr, annotation_splices_dict, padsize, regex_match = re.compile("^(.+?):(.+?)-(.+?)$")):
    
        """ locates any nearby annotated splice sites based on their positions using a 
        pre-seperated dictionary and binary search """
    
        keyvals = regex_match.match(keystr).groups()
        target_start = int(keyvals[1])
        target_stop = int(keyvals[2])
        
        nearest_splices=[]
        if keyvals[0] in annotation_splices_dict.keys():
            starts = annotation_splices_dict[keyvals[0]]["starts"]
            stops = annotation_splices_dict[keyvals[0]]["stops"]
        
            i = bisect_left(starts, target_start-padsize)
            
            while i<1E10:
                if i==len(starts):
                    break
                if starts[i]>target_start+padsize:
                    break
                if starts[i]>target_start-padsize and starts[i]<target_start+padsize and stops[i]>target_stop-padsize and stops[i]<target_stop+padsize:
                    nearest_splices.append("{}:{}-{}".format(keyvals[0], starts[i], stops[i]))
                i+=1
        return(nearest_splices)
    
    # setup the output data structures
    read_details={}
    splice_details={}
    splice_summary_numbers={"canonical_splices":0, "non_canonical_splices":0}
    if annotation_splices is not None and pwms is not None:
        splice_summary_numbers={"canonical_splices":0,
                                "non_canonical_splices":0,
                                "annotated_canonical_splices":0,
                                "annotated_canonical_undef_splices":0,
                                "annotated_canonical_U2_splices":0,
                                "annotated_canonical_U12_splices":0,
                                "annotated_non_canonical_splices":0,
                                "annotated_non_canonical_undef_splices":0,
                                "annotated_non_canonical_U2_splices":0,
                                "annotated_non_canonical_U12_splices":0,
                                "novel_canonical_splices":0,
                                "novel_canonical_undef_splices":0,
                                "novel_canonical_U2_splices":0,
                                "novel_canonical_U12_splices":0,
                                "novel_canonical_splices_with_nearby_annotated_canonical":0,
                                "novel_canonical_splices_with_nearby_annotated_non_canonical":0,
                                "novel_non_canonical_splices":0,
                                "novel_non_canonical_undef_splices":0,
                                "novel_non_canonical_U2_splices":0,
                                "novel_non_canonical_U12_splices":0,
                                "novel_non_canonical_splices_with_nearby_annotated_canonical":0,
                                "novel_non_canonical_splices_with_nearby_annotated_non_canonical":0,
                                "novel_non_canonical_splices_with_nearby_novel_canonical":0
                                }
    elif annotation_splices is not None and pwms is None:
        splice_summary_numbers={"canonical_splices":0,
                                "non_canonical_splices":0,
                                "annotated_canonical_splices":0,
                                "annotated_non_canonical_splices":0,
                                "novel_canonical_splices":0,
                                "novel_canonical_splices_with_nearby_annotated_canonical":0,
                                "novel_canonical_splices_with_nearby_annotated_non_canonical":0,
                                "novel_non_canonical_splices":0,
                                "novel_non_canonical_splices_with_nearby_annotated_canonical":0,
                                "novel_non_canonical_splices_with_nearby_annotated_non_canonical":0,
                                "novel_non_canonical_splices_with_nearby_novel_canonical":0,
                                }
    elif annotation_splices is None and pwms is not None:
        splice_summary_numbers={"canonical_splices":0,
                                "non_canonical_splices":0,
                                "canonical_undef_splices":0,
                                "canonical_U2_splices":0,
                                "canonical_U12_splices":0,
                                "non_canonical_undef_splices":0,
                                "non_canonical_U2_splices":0,
                                "non_canonical_U12_splices":0}
    
    if annotation_splices is not None:
        # split the annotation information by chromosome and position for efficient searching
        regex_match = re.compile("^(.+?):(.+?)-(.+?)$")
        split_annotation_splices = [regex_match.match(x).groups() for x in sorted(annotation_splices.keys())]
        annotation_splices_dict = {}
        for x,y,z in split_annotation_splices:
            try:
                annotation_splices_dict[x]["values"].append((int(y),int(z)))
            except KeyError:
                annotation_splices_dict[x]={"values":[]}
                annotation_splices_dict[x]["values"].append((int(y),int(z)))

        for key in annotation_splices_dict.keys():
            annotation_splices_dict[key]["values"].sort(key=lambda r: r[0])
            annotation_splices_dict[key]["starts"] = [r[0] for r in annotation_splices_dict[key]["values"]]
            annotation_splices_dict[key]["stops"] = [r[1] for r in annotation_splices_dict[key]["values"]]
    
    # Process the aligned reads looking for splices and classifying them
    nlogs=1
    counter=0
    ts = time.time()
    t0 = time.time()
    for read in read_iterator:
        
        current_read_pos = read.reference_start
        thisread_details={"canonical_splices":[],
                          "non_canonical_splices":[],
                          "is_reverse":False}
        
        if annotation_splices is not None:
            annot_details={"annotated_splices":[],
                           "novel_splices":[],
                           "transcripts_matching":{}}
            thisread_details.update(annot_details)
        
        if pwms is not None:
            classification_details={"unclassified_splices":[]}
            for key in pwms.keys():
                pwm = pwms[key]
                spliceend=None
                if "acceptor" in key:
                    spliceend = "acceptor"
                elif "donor" in key:
                    spliceend = "donor"
                splicetype = re.sub("_{}".format(spliceend),"",key)
                classification_details["{}_splices".format(splicetype)]=[]
            thisread_details.update(classification_details)
        
        # identify and process each splice in the read
        for tup in read.cigartuples:
            
            if tup[0]==3 and tup[1]>min_intron_length:
                
                # define the donor and acceptor splice regions in which we will search for alternative splices.
                donor_splice_site=None
                acceptor_splice_site=None
                strand = "+"
                if read.is_reverse:
                    strand = "-"
                    thisread_details["is_reverse"]=True
                    donor_splice_site = genome[read.reference_name][current_read_pos+tup[1]-2-splicepad:current_read_pos+tup[1]+splicepad].reverse_complement()
                    acceptor_splice_site = genome[read.reference_name][current_read_pos-splicepad:current_read_pos+2+splicepad].reverse_complement()
                else:
                    acceptor_splice_site = genome[read.reference_name][current_read_pos+tup[1]-2-splicepad:current_read_pos+tup[1]+splicepad]
                    donor_splice_site = genome[read.reference_name][current_read_pos-splicepad:current_read_pos+2+splicepad]
                
                # define the splice genomic coordinates and the terminal dinucleotides
                keystr = "{}:{}-{}".format(read.reference_name,
                                               current_read_pos,
                                               current_read_pos+tup[1])
                donor_splice_string = donor_splice_site.seq[splicepad:splicepad+2]
                acceptor_splice_string = acceptor_splice_site.seq[splicepad:splicepad+2]
                
                # if the splice has been seen before then we can just record that its been seen in a new read
                # otherwise we have to classify it.
                if keystr in splice_details:
                    splice_details[keystr]["reads"].append(read.query_name)
                else:
                    splice_details[keystr]={"reads":[read.query_name],
                                            "sites":(donor_splice_string, acceptor_splice_string)}
                    
                    # classify splice as cannonical or not.
                    if donor_splice_string!="GT" or acceptor_splice_string!="AG":
                        splice_details[keystr]["is_canonical"]=False
                    else:
                        splice_details[keystr]["is_canonical"]=True
                    
                    # classify the mapped site as U2/U12 - we only need to do this the first time this splice is seen
                    if pwms is not None:
                        classification, options = classifySplice(pwms, read.reference_name, current_read_pos,
                                                                 current_read_pos+tup[1], strand, genome,
                                                                 minthresh=pwmscorethreshold)
                        splice_details[keystr]["U2/U12_classification"]=classification
                        splice_details[keystr]["U2/U12_scores"]=options
                    
                    # classify if splice is annotated or not
                    if annotation_splices is not None:
                        if keystr in annotation_splices.keys():
                            splice_details[keystr]["is_annotated"]=True
                        else:
                            splice_details[keystr]["is_annotated"]=False
                            
                            # locate nearby annotated splice sites - the +1 here is so that we include the dinucleotide
                            # motif and then the pad region around it...
                            nearby_annotated_splices = locNearestAnnotated(keystr, annotation_splices_dict, splicepad+1)
                            
                            # if the splice is not annotated, search for nearby splice sites
                            nearby_canonical_splices = locNearestCanonicals(donor_splice_site.seq, acceptor_splice_site.seq,
                                                                            splicepad, read.reference_name, current_read_pos,
                                                                            current_read_pos+tup[1])
                            annot_alt_cannonical=[]
                            annot_alt_cannonical_classification=[]
                            annot_alt_cannonical_scores=[]
                            novel_alt_cannonical=[]
                            novel_alt_cannonical_classification=[]
                            novel_alt_cannonical_scores=[]
                            
                            for alt in nearby_canonical_splices:
                                annotated=False
                                
                                if alt in annotation_splices.keys():
                                    nearby_annotated_splices.remove(alt)
                                    annotated=True
                                
                                # classify the alternative splices as U2/U12
                                if pwms is not None:
                                    match = re.match("(.+):([0-9]+)-([0-9]+)", alt)
                                    classification, options = classifySplice(pwms, match.group(1), int(match.group(2)),
                                                                             int(match.group(3)), strand, genome,
                                                                             minthresh=pwmscorethreshold)
                                    if annotated:
                                        annot_alt_cannonical.append(alt)
                                        annot_alt_cannonical_classification.append(classification)
                                        annot_alt_cannonical_scores.append(options)
                                    else:
                                        novel_alt_cannonical.append(alt)
                                        novel_alt_cannonical_classification.append(classification)
                                        novel_alt_cannonical_scores.append(options)
                            
                            splice_details[keystr]["annotated_alt_canonical"]=annot_alt_cannonical
                            splice_details[keystr]["annotated_alt_non_canonical"]=nearby_annotated_splices
                            splice_details[keystr]["annotated_alt_canonical_U2/U12_classification"]=annot_alt_cannonical_classification
                            splice_details[keystr]["annotated_alt_canonical_U2/U12_scores"]=annot_alt_cannonical_scores
                            splice_details[keystr]["novel_alt_canonical"]=novel_alt_cannonical
                            splice_details[keystr]["novel_alt_canonical_U2/U12_classification"]=novel_alt_cannonical_classification
                            splice_details[keystr]["novel_alt_canonical_U2/U12_scores"]=novel_alt_cannonical_scores
                    
                    # build summary information
                    try:
                        
                        if splice_details[keystr]["is_canonical"]:
                            splice_summary_numbers["canonical_splices"]+=1
                        else:
                            splice_summary_numbers["non_canonical_splices"]+=1
                        
                        if "is_annotated" in splice_details[keystr].keys() and "U2/U12_classification" in splice_details[keystr].keys():
                            
                            if splice_details[keystr]["is_annotated"] and splice_details[keystr]["is_canonical"]:
                                splice_summary_numbers["annotated_canonical_splices"]+=1
                                if splice_details[keystr]["U2/U12_classification"] is None:
                                    splice_summary_numbers["annotated_canonical_undef_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U2"):
                                    splice_summary_numbers["annotated_canonical_U2_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U12"):
                                    splice_summary_numbers["annotated_canonical_U12_splices"]+=1
                            
                            elif splice_details[keystr]["is_annotated"] and not splice_details[keystr]["is_canonical"]:
                                splice_summary_numbers["annotated_non_canonical_splices"]+=1
                                if splice_details[keystr]["U2/U12_classification"] is None:
                                    splice_summary_numbers["annotated_non_canonical_undef_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U2"):
                                    splice_summary_numbers["annotated_non_canonical_U2_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U12"):
                                    splice_summary_numbers["annotated_non_canonical_U12_splices"]+=1
                            
                            elif not splice_details[keystr]["is_annotated"] and splice_details[keystr]["is_canonical"]:
                                splice_summary_numbers["novel_canonical_splices"]+=1
                                if len(splice_details[keystr]["annotated_alt_canonical"])>0:
                                    splice_summary_numbers["novel_canonical_splices_with_nearby_annotated_canonical"]+=1
                                if len(splice_details[keystr]["annotated_alt_non_canonical"])>0:
                                    splice_summary_numbers["novel_canonical_splices_with_nearby_annotated_non_canonical"]+=1
                                if splice_details[keystr]["U2/U12_classification"] is None:
                                    splice_summary_numbers["novel_canonical_undef_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U2"):
                                    splice_summary_numbers["novel_canonical_U2_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U12"):
                                    splice_summary_numbers["novel_canonical_U12_splices"]+=1
                            
                            elif not splice_details[keystr]["is_annotated"] and not splice_details[keystr]["is_canonical"]:
                                splice_summary_numbers["novel_non_canonical_splices"]+=1
                                if len(splice_details[keystr]["annotated_alt_canonical"])>0:
                                    splice_summary_numbers["novel_non_canonical_splices_with_nearby_annotated_canonical"]+=1
                                if len(splice_details[keystr]["annotated_alt_non_canonical"])>0:
                                    splice_summary_numbers["novel_non_canonical_splices_with_nearby_annotated_non_canonical"]+=1
                                if len(splice_details[keystr]["novel_alt_canonical"])>0:
                                    splice_summary_numbers["novel_non_canonical_splices_with_nearby_novel_canonical"]+=1
                                if splice_details[keystr]["U2/U12_classification"] is None:
                                    splice_summary_numbers["novel_non_canonical_undef_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U2"):
                                    splice_summary_numbers["novel_non_canonical_U2_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U12"):
                                    splice_summary_numbers["novel_non_canonical_U12_splices"]+=1
                        
                        elif "annotated" in splice_details[keystr].keys() and not "U2/U12_classification" in splice_details[keystr].keys():
                            
                            if splice_details[keystr]["is_annotated"] and splice_details[keystr]["is_canonical"]:
                                splice_summary_numbers["annotated_canonical_splices"]+=1
                            elif splice_details[keystr]["is_annotated"] and not splice_details[keystr]["is_canonical"]:
                                splice_summary_numbers["annotated_non_canonical_splices"]+=1
                            elif not splice_details[keystr]["is_annotated"] and splice_details[keystr]["is_canonical"]:
                                splice_summary_numbers["novel_canonical_splices"]+=1
                                if len(splice_details[keystr]["annotated_alt_canonical"])>0:
                                    splice_summary_numbers["novel_canonical_splices_with_nearby_annotated_canonical"]+=1
                                if len(splice_details[keystr]["annotated_alt_non_canonical"])>0:
                                    splice_summary_numbers["novel_canonical_splices_with_nearby_annotated_non_canonical"]+=1
                            elif not splice_details[keystr]["is_annotated"] and not splice_details[keystr]["is_canonical"]:
                                splice_summary_numbers["novel_non_canonical_splices"]+=1
                                if len(splice_details[keystr]["annotated_alt_canonical"])>0:
                                    splice_summary_numbers["novel_non_canonical_splices_with_nearby_annotated_canonical"]+=1
                                if len(splice_details[keystr]["annotated_alt_non_canonical"])>0:
                                    splice_summary_numbers["novel_non_canonical_splices_with_nearby_annotated_non_canonical"]+=1
                                if len(splice_details[keystr]["novel_alt_canonical"])>0:
                                    splice_summary_numbers["novel_non_canonical_splices_with_nearby_novel_canonical"]+=1
                        
                        elif "annotated" not in splice_details[keystr].keys() and "U2/U12_classification" in splice_details[keystr].keys():
                            
                            if splice_details[keystr]["is_canonical"]:
                                if splice_details[keystr]["U2/U12_classification"] is None:
                                    splice_summary_numbers["canonical_undef_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U2"):
                                    splice_summary_numbers["canonical_U2_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U12"):
                                    splice_summary_numbers["canonical_U12_splices"]+=1
                            
                            elif not splice_details[keystr]["is_canonical"]:
                                if splice_details[keystr]["U2/U12_classification"] is None:
                                    splice_summary_numbers["non_canonical_undef_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U2"):
                                    splice_summary_numbers["non_canonical_U2_splices"]+=1
                                elif splice_details[keystr]["U2/U12_classification"][0].startswith("U12"):
                                    splice_summary_numbers["non_canonical_U12_splices"]+=1

                    except:
                        print("splice_details")
                        for key in sorted(splice_details.keys()):
                            print(key, splice_details[key])
                        print("splice_summary_numbers")
                        for key in sorted(splice_summary_numbers.keys()):
                            print(key, splice_summary_numbers[key])
                        raise
                
                if splice_details[keystr]["is_canonical"] :
                    thisread_details["canonical_splices"].append(keystr)
                else:
                    thisread_details["non_canonical_splices"].append(keystr)
                
                if annotation_splices is not None:
                    if splice_details[keystr]["is_annotated"]:
                        thisread_details["annotated_splices"].append(keystr)
                        thisread_details["transcripts_matching"][keystr] = annotation_splices[keystr]["transcripts"]
                    else:
                        thisread_details["novel_splices"].append(keystr)
                
                if pwms is not None:
                    if splice_details[keystr]["U2/U12_classification"] is not None:
                        thisread_details["{}_splices".format(splice_details[keystr]["U2/U12_classification"][0])].append(keystr)
                    else:
                        thisread_details["unclassified_splices"].append(keystr)
                    
                #print(read.query_name, current_read_pos, current_read_pos+tup[1], donor_splice_site.seq, acceptor_splice_site.seq, read_is_all_canonical)
                current_read_pos +=  tup[1]
            elif tup[0]==0 or tup[0]==2:
                current_read_pos +=  tup[1]
        
        read_details[read.query_name] = thisread_details
        
        counter+=1
       
        if (counter % LOG_EVERY_N)==0:
            msg="processed {these} reads (dt = {sec:.2f}s) ...".format(these=(nlogs*LOG_EVERY_N), sec=time.time()-t0)
            if logger is not None:
                logger.info(msg)
            else:
                print(msg)
            nlogs+=1
            t0=time.time()
    
    msg = "Finished processed {these} reads (dt = {sec:.2f}s).".format(these=(nlogs*LOG_EVERY_N)+counter,
                                                                       sec=time.time()-t0)
    
    return(read_details, splice_details, splice_summary_numbers)

def getAnnotationIntrons(annot, genome, chr_synonym_dic={}, logger=None, LOG_EVERY_N=10000):
    
    """ return a dictionary with all the introns in a given annotation """
    
    if logger is not None:
        logger.info("parsing transcript exon structures....")
    
    annot.clear_all()
    annot.set_feature("exons")
    exons = annot.get_selection()
    tx_exons={}
    for exon in exons:
        if type(exon.desc["parent"]) is list:
            for parent in exon.desc["parent"]:
                try:
                    tx_exons[parent].append(exon)
                except KeyError:
                    tx_exons[parent]=[exon]
        else:
            try:
                tx_exons[exon.desc["parent"]].append(exon)
            except KeyError:
                tx_exons[exon.desc["parent"]]=[exon]
    
    splice_details={}
    nlogs=1
    counter=0
    for transcript in tx_exons.keys():
        these_introns = computeIntrons(tx_exons[transcript], logger=logger)
        for intron in these_introns:
            donor_splice_site=None
            acceptor_splice_site=None
            
            this_chrid = intron.chrid
            if intron.chrid not in genome.keys():
                if intron.chrid in chr_synonym_dic.keys() and chr_synonym_dic[intron.chrid] in genome.keys():
                    this_chrid = chr_synonym_dic[intron.chrid]
                else:
                    msg = "There is a mismatch between the annotation and " \
                          "genome chromosome IDs that is not accounted for " \
                          "in the provided chromosome synonyms list. Details: " \
                          "Annotation ID: {}, genome IDs: {}, synonyms: {} " \
                          "".format(intron.chrid, genome.keys(), chr_synonym_dic)
                    raise ValueError(msg)
            
            if intron.strand=="-":
                donor_splice_site = genome[this_chrid][intron.stop-2:intron.stop].reverse_complement()
                acceptor_splice_site = genome[this_chrid][intron.start-1:intron.start+1].reverse_complement()
            else:
                acceptor_splice_site = genome[this_chrid][intron.stop-2:intron.stop]
                donor_splice_site = genome[this_chrid][intron.start-1:intron.start+1]
            
            keystr = "{}:{}-{}".format(intron.chrid,
                                       intron.start-1,
                                       intron.stop)
            
            try:
                splice_details[keystr]["transcripts"].append(transcript)
            except:
                splice_details[keystr]={"transcripts":[transcript],
                                        "sites":(str(donor_splice_site.seq),
                                                 str(acceptor_splice_site.seq))}
                if donor_splice_site.seq!="GT" or acceptor_splice_site.seq!="AG":
                    splice_details[keystr]["is_canonical"]=False
                else:
                    splice_details[keystr]["is_canonical"]=True
        
        counter+=1
        
        if (counter % LOG_EVERY_N)==0:
            msg="processed {these} transcripts...".format(these=(nlogs*LOG_EVERY_N))
            if logger is not None:
                logger.info(msg)
            else:
                print(msg)
            nlogs+=1
    
    return(splice_details)

def writeBamFiles(thisbam, splice_read_details, prefix, splits, pwms=None, logger=None):
     
    """ write bam files separating out the reads into their categories """
    
    annot_opts = ["annotated","novel"]
    cann_opts = ["canonical","non_canonical"]
    
    pwm_opts=["unclassified"]
    if pwms is not None:
        for key in pwms.keys():
            spliceend=None
            if "acceptor" in key:
                spliceend = "acceptor"
                splicetype = re.sub("_{}".format(spliceend),"",key)
                pwm_opts.append(splicetype)
    
    bamfiles = {}
    if "annotated" in splits and "U2" in splits and "canonical" in splits:
        combos = list(itertools.product(annot_opts, pwm_opts, cann_opts))
        for combo in combos:
            pwmcombo = re.sub("_","",combo[1])
            bamfile = "{}{}_{}_{}_splices.bam".format(prefix, combo[0], pwmcombo, combo[2])
            bamfiles[combo]={"filename":bamfile,
                             "filehandle":pysam.AlignmentFile(bamfile, "wb", template=thisbam)}
    elif "annotated" in splits and "U2" in splits and "canonical" not in splits:
        combos = list(itertools.product(annot_opts, pwm_opts))
        for combo in combos:
            pwmcombo = re.sub("_","",combo[1])
            bamfile = "{}{}_{}_splices.bam".format(prefix, combo[0], pwmcombo)
            bamfiles[combo]={"filename":bamfile,
                             "filehandle":pysam.AlignmentFile(bamfile, "wb", template=thisbam)}
    elif "annotated" in splits and "U2" not in splits and "canonical" in splits:
        combos = list(itertools.product(annot_opts, cann_opts))
        for combo in combos:
            bamfile = "{}{}_{}_splices.bam".format(prefix, combo[0], combo[1])
            bamfiles[combo]={"filename":bamfile,
                             "filehandle":pysam.AlignmentFile(bamfile, "wb", template=thisbam)}
    elif "annotated" not in splits and "U2" in splits and "canonical" in splits:
        combos = list(itertools.product(pwm_opts, cann_opts))
        for combo in combos:
            pwmcombo = re.sub("_","",combo[0])
            bamfile = "{}{}_{}_splices.bam".format(prefix, pwmcombo, combo[1])
            bamfiles[combo]={"filename":bamfile,
                             "filehandle":pysam.AlignmentFile(bamfile, "wb", template=thisbam)}
    elif "annotated" in splits and "U2" not in splits and "canonical" not in splits:
        combos = annot_opts
        for combo in combos:
            bamfile = "{}{}_splices.bam".format(prefix, combo[0])
            bamfiles[combo]={"filename":bamfile,
                             "filehandle":pysam.AlignmentFile(bamfile, "wb", template=thisbam)}
    elif "annotated" not in splits and "U2" in splits and "canonical" not in splits:
        combos = pwm_opts
        for combo in combos:
            pwmcombo = re.sub("_","",combo[0])
            bamfile = "{}{}_splices.bam".format(prefix, pwmcombo)
            bamfiles[combo]={"filename":bamfile,
                             "filehandle":pysam.AlignmentFile(bamfile, "wb", template=thisbam)}
    elif "annotated" not in splits and "U2" not in splits and "canonical" in splits:
        combos = cann_opts
        for combo in combos:
            bamfile = "{}{}_splices.bam".format(prefix, combo[0])
            bamfiles[combo]={"filename":bamfile,
                             "filehandle":pysam.AlignmentFile(bamfile, "wb", template=thisbam)}
    
    read_summary_numbers={}
    
    for read in thisbam.fetch():
        if read.query_name in splice_read_details.keys():
            this_read = splice_read_details[read.query_name]
            for combo in combos:
                writetothiscombo=True
                for val in combo:
                    if len(this_read["{}_splices".format(val)])==0:
                        writetothiscombo=False
                if writetothiscombo:
                    bamfiles[combo]["filehandle"].write(read)
                    previouslywritten = combo
                    try:
                        read_summary_numbers[bamfiles[combo]["filename"]]+=1
                    except KeyError:
                        read_summary_numbers[bamfiles[combo]["filename"]]=1
    
    return(read_summary_numbers)

def calAltDists(mapsplice, altsplices):
    mapsplicematch = re.match(".+:([0-9]+)-([0-9]+)", mapsplice)
    dists=[]
    for val in altsplices:
        altsplicematch = re.match(".+:([0-9]+)-([0-9]+)", val)
        thisdist = abs(int(altsplicematch.group(1))-(int(mapsplicematch.group(1))))+abs(int(altsplicematch.group(2))-(int(mapsplicematch.group(2))))
        dists.append(thisdist)
    dists=numpy.array(dists)
    return(dists)

def mkAnnotPWMBoxPlots(splice_details, plotfile, pointsize=4, title="", legend_fontsize=10, legendloc=(0.2,0.7),
                       logger=None):
    
    """ plot the results of the analysis as a set of cool box plots..."""
    
    if logger is not None:
        logger.info("Plotting classified annotated splicing details boxplots...")
    
    annotated_splice_exprs_data={}
    novel_splice_exprs_data={}
    alt_splice_dists = {}
    alt_splice_annotated = {}
    
    for splice in splice_details.keys():
        this_splice = splice_details[splice]
        exprs = len(this_splice["reads"])
        splice_class = "None"
        if this_splice['U2/U12_classification'] is not None:
            splice_class = this_splice['U2/U12_classification'][0]
        if this_splice['is_annotated']:
            try:
                annotated_splice_exprs_data[splice_class].append(exprs)
            except KeyError:
                annotated_splice_exprs_data[splice_class]=[exprs]
        else:
            dist=0
            annot=False
            
            if "annotated_alt_canonical" in splice_details[splice].keys():
                if len(splice_details[splice]["annotated_alt_canonical"])>0:
                    dist = calAltDists(splice, splice_details[splice]["annotated_alt_canonical"]).min()
                    annot=True
            
            if not annot and "novel_alt_canonical" in splice_details[splice].keys():
                if len(splice_details[splice]["novel_alt_canonical"])>0:
                    dist = calAltDists(splice, splice_details[splice]["novel_alt_canonical"]).min()
            
            try:
                novel_splice_exprs_data[splice_class].append(exprs)
                alt_splice_dists[splice_class].append(dist)
                alt_splice_annotated[splice_class].append(annot)
            except KeyError:
                novel_splice_exprs_data[splice_class]=[exprs]
                alt_splice_dists[splice_class]=[dist]
                alt_splice_annotated[splice_class]=[annot]
    
    fig=plt.figure(figsize=(15,10))
    
    data = []
    labels = list(set(sorted(annotated_splice_exprs_data.keys())+sorted(novel_splice_exprs_data.keys())))
    for label in labels:
        if label in annotated_splice_exprs_data.keys():
            data.append(numpy.log10(numpy.array(annotated_splice_exprs_data[label])))
        else:
            data.append(numpy.array([]))
        if label in novel_splice_exprs_data.keys():
            data.append(numpy.log10(numpy.array(novel_splice_exprs_data[label])))
        else:
            data.append(numpy.array([]))
    
    plt.ylim((-0.1,5))
    axtickspos = (numpy.arange(len(labels))*2)+1.5
    boxpos=[]
    for val in axtickspos:
        boxpos.append(val-0.25)
        boxpos.append(val+0.25)
    
    plotaxs=[]
    plotays=[]
    plotuaxs=[]
    plotuays=[]
    plotann=[]
    plotzerodist=[]
    
    i=0
    while i<len(boxpos):
        ay = data[i]
        uay = data[i+1]
        ax = numpy.random.normal(boxpos[i], 0.05, len(ay))
        
        if len(uay)==0:
            xdata = []
        else:
            xdata = numpy.array(alt_splice_dists[labels[int((i+1)/2)]])
        for val in xdata:
            if val==0:
                plotzerodist.append(True)
            else:
                plotzerodist.append(False)
                
        scalewidth=0.4
        if len(xdata)>0:
            if xdata.max()==0:
                uax=xdata
            else:
                uax = boxpos[i+1]+((xdata/xdata.max())*scalewidth)-(scalewidth/2)
        else:
            uax=xdata
        
        for xval in ax:
            plotaxs.append(xval)
        for xval in uax:
            plotuaxs.append(xval)
        for yval in ay:
            plotays.append(yval)
        for yval in uay:
            plotuays.append(yval)
        
        if len(uay)>0:
            for val in alt_splice_annotated[labels[int((i+1)/2)]]:
                plotann.append(val)
        i+=2
    
    plotaxs=numpy.array(plotaxs)
    plotays=numpy.array(plotays)
    plotuaxs=numpy.array(plotuaxs)
    plotuays=numpy.array(plotuays)
    plotann=numpy.array(plotann)
    plotzerodist=numpy.array(plotzerodist)
    
    i=0
    j=0
    while i<len(labels):
        if labels[i] in annotated_splice_exprs_data.keys():
            plt.text(boxpos[j], 4.6, sum(numpy.array(annotated_splice_exprs_data[labels[i]])), size=legend_fontsize, ha="center", color="indianred")
            plt.text(boxpos[j], 4.8, len(numpy.array(annotated_splice_exprs_data[labels[i]])), size=legend_fontsize, ha="center", color="indianred")
        if labels[i] in novel_splice_exprs_data.keys():
            plt.text(boxpos[j+1], 4.6, sum(numpy.array(novel_splice_exprs_data[labels[i]])), size=legend_fontsize, ha="center", color="steelblue")
            plt.text(boxpos[j+1], 4.8, len(numpy.array(novel_splice_exprs_data[labels[i]])), size=legend_fontsize, ha="center", color="steelblue")
        i+=1
        j+=2

    bp = plt.boxplot(data, sym="", widths=0.2, positions=boxpos)
    acaregories = plt.scatter(plotaxs, plotays, s=pointsize, c="indianred", alpha=0.3, label="annotated splices")
    if len(plotann)>0:
        uacaregories_uanearby = plt.scatter(plotuaxs[numpy.invert(plotann)], plotuays[numpy.invert(plotann)], s=pointsize, c='y', alpha=1.0, label="novel splices with novel canonical splices nearby")
        uacaregories_anearby = plt.scatter(plotuaxs[plotann], plotuays[plotann], s=pointsize, c='darkcyan', alpha=0.3, label="novel splices with annotated canonical splices nearby")
    if len(plotzerodist)>0:
        uacaregories_nonearby = plt.scatter(plotuaxs[plotzerodist], plotuays[plotzerodist], s=pointsize, c='black', alpha=1.0, label="novel splices with no canonical splices nearby")
    
    lw=1.5
    for key in bp.keys():
        for box in bp[key]:
            box.set(linewidth=lw, color='0.4')
    
    plt.ylabel(r"$log_{10}(counts)$")
    ax = plt.gca()
    ax.set_xticks(axtickspos)
    xtixks = ax.set_xticklabels([re.sub("_"," ",x) for x in labels], rotation=20, ha="right", size=10)
    plt.title(title)
    plt.legend(loc=legendloc, fontsize=legend_fontsize)
    plt.tight_layout()
    
    plt.savefig(plotfile, dpi=300, format="svg")

def mkPWMBoxPlots(splice_details, plotfile, pointsize=4, title="", logger=None):
    
    """ plot the results of the analysis as a set of cool box plots..."""
    
    if logger is not None:
        logger.info("Plotting classified splicing details boxplots...")
    
    splice_exprs_data = {}
    
    for splice in splice_details.keys():
        this_splice = splice_details[splice]
        exprs = len(this_splice["reads"])
        splice_class = "None"
        if this_splice['U2/U12_classification'] is not None:
            splice_class = this_splice['U2/U12_classification'][0]
        try:
            splice_exprs_data[splice_class].append(exprs)
        except KeyError:
            splice_exprs_data[splice_class]=[exprs]
    
    fig=plt.figure(figsize=(15,10))
    
    data = []
    labels = list(sorted(splice_exprs_data.keys()))
    for label in labels:
        data.append(numpy.log10(numpy.array(splice_exprs_data[label])))
    
    plt.ylim((-0.1,5))
    boxpos = (numpy.arange(len(labels))*2)+1.5
    
    plotxs=[]
    plotys=[]
    
    i=0
    while i<len(boxpos):
        y = data[i]
        x = numpy.random.normal(boxpos[i], 0.05, len(y))
                
        for xval in x:
            plotxs.append(xval)
        for yval in y:
            plotys.append(yval)
        i+=1
    
    plotxs=numpy.array(plotxs)
    plotys=numpy.array(plotys)
    
    i=0
    j=0
    while i<len(labels):
        if labels[i] in splice_exprs_data.keys():
            plt.text(boxpos[j], 4.6, sum(numpy.array(splice_exprs_data[labels[i]])), size=10, ha="center", color="indianred")
            plt.text(boxpos[j], 4.8, len(numpy.array(splice_exprs_data[labels[i]])), size=10, ha="center", color="indianred")
        i+=1
        j+=1

    bp = plt.boxplot(data, sym="", widths=0.2, positions=boxpos)
    caregories = plt.scatter(plotxs, plotys, s=pointsize, c="indianred", alpha=0.3)
    
    lw=1.5
    for key in bp.keys():
        for box in bp[key]:
            box.set(linewidth=lw, color='0.4')
    
    plt.ylabel(r"$log_{10}(counts)$")
    ax = plt.gca()
    ax.set_xticks(boxpos)
    xtixks = ax.set_xticklabels([re.sub("_"," ",x) for x in labels], rotation=20, ha="right", size=10)
    plt.title(title)
    plt.tight_layout()
    
    plt.savefig(plotfile, dpi=300, format="svg")

def mkAnnotBoxPlots(splice_details, plotfile, pointsize=4, title="", legend_fontsize=8, logger=None):
    
    """ plot the results of the analysis as a set of cool box plots..."""
    
    if logger is not None:
        logger.info("Plotting annotated splicing details boxplots...")
    
    splice_exprs_data={"annotated_canonical":[],
                       "annotated_non_canonical":[],
                       "novel_canonical_with_nearby_annotated_canonical":[],
                       "novel_canonical_no_nearby_annotated_canonical":[],
                       "novel_non_canonical_with_nearby_annotated_canonical":[],
                       "novel_non_canonical_with_nearby_novel_canonical":[],
                       "novel_non_canonical_no_nearby_canonical":[]}
    
    alt_splice_distances = {"novel_canonical_with_nearby_annotated_canonical":[],
                            "novel_non_canonical_with_nearby_annotated_canonical":[],
                            "novel_non_canonical_with_nearby_novel_canonical":[]} 
    
    for splice in splice_details.keys():
        exprs = len(splice_details[splice]["reads"])
        
        if splice_details[splice]['is_annotated'] and splice_details[splice]['is_canonical']:
            splice_exprs_data["annotated_canonical"].append(exprs)
        
        elif splice_details[splice]['is_annotated'] and not splice_details[splice]['is_canonical']:
            splice_exprs_data["annotated_non_canonical"].append(exprs)
        
        elif splice_details[splice]['is_canonical'] and not splice_details[splice]['is_annotated']:
            if "annotated_alt_canonical" in splice_details[splice].keys():
                if len(splice_details[splice]["annotated_alt_canonical"])>0:
                    splice_exprs_data["novel_canonical_with_nearby_annotated_canonical"].append(exprs)
                    dists = calAltDists(splice, splice_details[splice]["annotated_alt_canonical"])
                    alt_splice_distances["novel_canonical_with_nearby_annotated_canonical"].append(dists.min())
                else:
                    splice_exprs_data["novel_canonical_no_nearby_annotated_canonical"].append(exprs)
            else:
                splice_exprs_data["novel_canonical_no_nearby_annotated_canonical"].append(exprs)
        
        elif not splice_details[splice]['is_canonical'] and not splice_details[splice]['is_annotated']:
            if "annotated_alt_canonical" in splice_details[splice].keys():
                if len(splice_details[splice]["annotated_alt_canonical"])>0:
                    splice_exprs_data["novel_non_canonical_with_nearby_annotated_canonical"].append(exprs)
                    dists = calAltDists(splice, splice_details[splice]["annotated_alt_canonical"])
                    alt_splice_distances["novel_non_canonical_with_nearby_annotated_canonical"].append(dists.min())
                elif "novel_alt_canonical" in splice_details[splice].keys():
                    if len(splice_details[splice]["novel_alt_canonical"])>0:
                        splice_exprs_data["novel_non_canonical_with_nearby_novel_canonical"].append(exprs)
                        dists = calAltDists(splice, splice_details[splice]["novel_alt_canonical"])
                        alt_splice_distances["novel_non_canonical_with_nearby_novel_canonical"].append(dists.min())
                    else:
                        splice_exprs_data["novel_non_canonical_no_nearby_canonical"].append(exprs)
            elif "novel_alt_canonical" in splice_details[splice].keys():
                if len(splice_details[splice]["novel_alt_canonical"])>0:
                    splice_exprs_data["novel_non_canonical_with_nearby_novel_canonical"].append(exprs)
                    dists = calAltDists(splice, splice_details[splice]["novel_alt_canonical"])
                    alt_splice_distances["novel_non_canonical_with_nearby_novel_canonical"].append(dists.min())
                else:
                    splice_exprs_data["novel_non_canonical_no_nearby_canonical"].append(exprs)
    
    fig=plt.figure(figsize=(15,10))
    
    data = []
    labels = list(splice_exprs_data.keys())
    for label in labels:
        data.append(numpy.log10(numpy.array(splice_exprs_data[label])))
    
    plt.ylim((-0.1,5))
    
    plotxs=[]
    plotys=[]
    cplotxs=[]
    cplotys=[]
    cplotcol=[]
    for i in numpy.arange(len(labels)):
        y = data[i]
        if labels[i] in alt_splice_distances.keys() and len(y)>0:
            #xvals = numpy.random.normal(i+1, 0.05, len(y))
            xdata = numpy.array(alt_splice_distances[labels[i]])
            scalewidth=0.4
            xvals = (i+1)+((xdata/xdata.max())*scalewidth)-(scalewidth/2)
            for xval in xvals:
                cplotxs.append(xval)
            for yval in y:
                cplotys.append(yval)
            for cval in alt_splice_distances[labels[i]]:
                cplotcol.append(cval)
        else:
            x = numpy.random.normal(i+1, 0.05, len(y))
            for xval in x:
                plotxs.append(xval)
            for yval in y:
                plotys.append(yval)
        plt.text(i+1, 4.8, len(numpy.array(splice_exprs_data[labels[i]])), size=10, ha="center")
        plt.text(i+1, 4.6, sum(numpy.array(splice_exprs_data[labels[i]])), size=10, ha="center")
    
    nearbys = plt.scatter(cplotxs, cplotys, s=pointsize, c=cplotcol, alpha=0.3, cmap="viridis_r")
    cbar = plt.colorbar(shrink=0.4, pad=0.02)
    cbar.ax.tick_params(labelsize=legend_fontsize)
    cbar.set_label('distance to nearest alternative (bp)', rotation=270, labelpad=13, fontsize=legend_fontsize)
    caregories = plt.scatter(plotxs, plotys, s=pointsize, c="indianred", alpha=0.3)
    bp = plt.boxplot(data, sym="")
    
    lw=1.5
    for key in bp.keys():
        for box in bp[key]:
            box.set(linewidth=lw, color='0.4')
    
        plt.ylabel(r"$log_{10}(counts)$")
    ax = plt.gca()
    xticks = ax.set_xticklabels([re.sub("_"," ",x) for x in labels], rotation=20, ha="right", size=10)
    plt.title(title)
    plt.tight_layout()
    
    plt.savefig(plotfile, dpi=300, format="svg")

def mkBoxPlots(splice_details, plotfile, pointsize=4, title="", logger=None):
    
    """ plot the results of the analysis as a set of cool box plots..."""
    
    if logger is not None:
        logger.info("Plotting basic splicing details boxplots...")
    
    splice_exprs_data = {"canonical":[],
                         "non_canonical":[]}
    
    for splice in splice_details.keys():
        this_splice = splice_details[splice]
        exprs = len(this_splice["reads"])
        splice_class = "None"
        try:
            if this_splice['is_canonical']:
                splice_exprs_data["canonical"].append(exprs)
            else:
                splice_exprs_data["non_canonical"].append(exprs)
        except KeyError:
            if this_splice['is_canonical']:
                splice_exprs_data["canonical"]=[exprs]
            else:
                splice_exprs_data["non_canonical"]=[exprs]
    
    fig=plt.figure(figsize=(15,10))
    
    data = []
    labels = list(sorted(splice_exprs_data.keys()))
    for label in labels:
        data.append(numpy.log10(numpy.array(splice_exprs_data[label])))
    
    plt.ylim((-0.1,5))
    boxpos = (numpy.arange(len(labels))*2)+1.5
    
    plotxs=[]
    plotys=[]
    
    i=0
    while i<len(boxpos):
        y = data[i]
        x = numpy.random.normal(boxpos[i], 0.05, len(y))
                
        for xval in x:
            plotxs.append(xval)
        for yval in y:
            plotys.append(yval)
        i+=1
    
    plotxs=numpy.array(plotxs)
    plotys=numpy.array(plotys)
    
    i=0
    j=0
    while i<len(labels):
        if labels[i] in splice_exprs_data.keys():
            plt.text(boxpos[j], 4.6, sum(numpy.array(splice_exprs_data[labels[i]])), size=10, ha="center", color="indianred")
            plt.text(boxpos[j], 4.8, len(numpy.array(splice_exprs_data[labels[i]])), size=10, ha="center", color="indianred")
        i+=1
        j+=1

    bp = plt.boxplot(data, sym="", widths=0.2, positions=boxpos)
    caregories = plt.scatter(plotxs, plotys, s=pointsize, c="indianred", alpha=0.3)
    
    lw=1.5
    for key in bp.keys():
        for box in bp[key]:
            box.set(linewidth=lw, color='0.4')
    
    plt.ylabel(r"$log_{10}(counts)$")
    ax = plt.gca()
    ax.set_xticks(boxpos)
    xtixks = ax.set_xticklabels([re.sub("_"," ",x) for x in labels], rotation=20, ha="right", size=10)
    plt.title(title)
    plt.tight_layout()
    
    plt.savefig(plotfile, dpi=300, format="svg")

def makePWM(pwmtext):
    
    nreg = '\s-?[1-9]+[0-9]*.?[0-9]*E-?\+?[0-9]+\s?'
    pwmdict={"A":[], "T":[], "G":[], "C":[]}
    alphorder = None
    for line in pwmtext:
        if alphorder is not None:
            vals = [float(x) for x in line.strip().split("\t")]
            i=0
            while i<len(alphorder):
                pwmdict[alphorder[i]].append(vals[i])
                i+=1

        alphmatch = re.match("(A|T|G|C)\t(A|T|G|C)\t(A|T|G|C)\t(A|T|G|C)\n", line)
        if alphmatch:
            alphorder = alphmatch.groups()
    
    pwm = motifs.matrix.PositionWeightMatrix(unambiguous_dna, pwmdict)
    return(pwm)

class Encoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (numpy.int_, numpy.intc, numpy.intp, numpy.int8,
                            numpy.int16, numpy.int32, numpy.int64, numpy.uint8,
                            numpy.uint16, numpy.uint32, numpy.uint64)):
            return int(obj)
        elif isinstance(obj, (numpy.float_, numpy.float16, numpy.float32, 
                              numpy.float64)):
            return float(obj)
        elif isinstance(obj,(numpy.ndarray,)): #### This is the fix
            return obj.tolist()
        elif isinstance(obj, Seq):
            return(str(obj))
        return json.JSONEncoder.default(self, obj)

if __name__ == '__main__':

    # parse command line options
    # Set standard parser
    parser, pos_args, kw_args = sp.standard_parser(__version__,
                                                   prog = __scriptname__, 
                                                   usage = __usage__,
                                                   description = __progdesc__,
                                                   epilog = __progepi__,
                                                   infile = False,
                                                   outfile = False,
                                                   tmpdir = False)
        
    parser, pos_args, kw_args = addScriptOptions(parser, pos_args, kw_args)
    
    args = parser.parse_args()
    
    # setup standard logging information
    script_logger = sl.standard_logger(__version__, sys.argv, args, pos_args, 
                                       kw_args, script_name=__scriptname__)
    
    splits = args.spliton.split(",")
    
    # ok, first lets check that the parameters specified are consistent.
    if args.annotation is None and args.splitreads and "annotated" in splits:
        script_logger.warn("Warning: cannot split on annotation because no annotation was provided.")
        splits.remove("annotated")
    
    if not args.pwm and "U2" in splits:
        script_logger.warn("Warning: cannot split on U2/U12 classification because no " \
                           "position weight matrices were provided.")
        splits.remove("U2")
    
    # load the genome
    script_logger.info("loading genome....")
    genome = {}
    for seq_record in SeqIO.parse(args.genomefile, "fasta"):
        genome[seq_record.name] = seq_record
    
    # load the annotation and enable novel splicing detection
    annot_details = None
    if args.annotation is not None:
        script_logger.info("Novel splicing detection enabled, loading annotation...")
        annot = annotation(args.annotation, filetype=args.input_format, stripChr=args.stripchr)
        
        # build a quick dictionary of chromosome synonyms to check for....
        chr_synonym_dic = {}
        for val in args.chr_synonyms.split(","):
            vals = val.split(":")
            if len(vals)==2:
                chr_synonym_dic[vals[0]]=vals[1]
                chr_synonym_dic[vals[1]]=vals[0]
        
        annot_details = getAnnotationIntrons(annot, genome, chr_synonym_dic,
                                             logger=script_logger)
    
        # load the annotation and enable novel splicing detection
    pwm_details = None
    if args.pwm:
        script_logger.info("U2/U12 calssification detection enabled, loading pwms...")
        PWMfiles = glob.glob("{}*.pwm".format(os.path.join(args.pwm_path,args.pwm_species)))
        pwm_details={}
        for PWMfile in PWMfiles:
            with open(PWMfile,"r") as fh:
                pwmkey = re.sub("{}_".format(args.pwm_species), "", os.path.splitext(os.path.basename(PWMfile))[0])
                pwm_details[pwmkey] = makePWM(fh.readlines())
    
    # classify intron splices    
    script_logger.info("loading bamfile and classifying data....")
    thisbam = pysam.AlignmentFile(args.bamfile, "rb")

    read_details, splice_details, splice_summary_numbers = charcterizeIntrons((read for read in thisbam.fetch()),
                                                                               genome,
                                                                               annotation_splices = annot_details,
                                                                               splicepad = args.altsplicepad,
                                                                               min_intron_length = args.minintronsize,
                                                                               pwms = pwm_details,
                                                                               pwmscorethreshold = args.pwm_thresh,
                                                                               logger=script_logger)
    
    # get output file prefix sorted
    prefix = args.prefix
    if prefix=="":
        prefix = os.path.splitext(args.bamfile)[0]
    
    plotfile = "{}spliceplot.svg".format(prefix)
    if args.pwm and args.annotation is not None:
        mkAnnotPWMBoxPlots(splice_details, plotfile, pointsize=4, title=os.path.basename(args.bamfile),
                           legend_fontsize=10, legendloc=(0.2,0.7), logger=script_logger)
    elif args.pwm and args.annotation is None:
        mkPWMBoxPlots(splice_details, plotfile, pointsize=4, title=os.path.basename(args.bamfile),
                      logger=script_logger)
    elif not args.pwm and args.annotation is not None:
        mkAnnotBoxPlots(splice_details, plotfile, pointsize=4, title=os.path.basename(args.bamfile),
                        legend_fontsize=8, logger=script_logger)
    elif not args.pwm and args.annotation is None:
        mkBoxPlots(splice_details, plotfile, pointsize=4, title=os.path.basename(args.bamfile)
                   , logger=script_logger)
    
    # log some details
    script_logger.info("Dataset {bamfile} contains {nsplices} detected splicing event " \
                       "(from {nreads} reads encompassing one or more splicing events)." \
                       "".format(bamfile = args.bamfile, nsplices = len(splice_details.keys()),
                                 nreads = len(read_details.keys())))
    
    strs=["The breakdown of the splicing events is:"]
    
    for key in sorted(splice_summary_numbers.keys()):
        thisstr = "\t{}:\t{} ({:.2%})".format(re.sub("_"," ",key),
                                              splice_summary_numbers[key],
                                              splice_summary_numbers[key]/len(splice_details.keys()))
        strs.append(thisstr)
    
    script_logger.info("\n".join(strs))
    
    strs=["Splices have been split into different alignment files:"]
    if args.splitreads:
        read_summary_numbers = writeBamFiles(thisbam, read_details, prefix, splits, pwms=pwm_details, logger=script_logger)
        for key in sorted(read_summary_numbers.keys()):
            thisstr = "\t{}:\t{} ({:.2%})".format(key, read_summary_numbers[key],
                                                  read_summary_numbers[key]/len(read_details.keys()))
            strs.append(thisstr)
    
    script_logger.info("\n".join(strs))
    
    script_logger.info("Graphical summary of the splicing details has been written to {}".format(plotfile))
     
    nc_stats_file = "{}splice_stats.json".format(prefix)
    script_logger.info("writing details of all detected splices to {}...".format(nc_stats_file))
    fh = open(nc_stats_file,"w")
    json.dump(splice_details, fh, sort_keys=True, indent=4, cls=Encoder)
    fh.close()
     
    nc_stats_file = "{}read_stats.json".format(prefix)
    script_logger.info("writing details of the splices in each read to {}...".format(nc_stats_file))
    fh = open(nc_stats_file,"w")
    json.dump(read_details, fh, sort_keys=True, indent=4, cls=Encoder)
    fh.close()
     
    script_logger.info("Finished. Have a nice day and don't forget to index the new bam files! ;)")
