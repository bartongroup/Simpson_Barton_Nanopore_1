'''
----------------------------------------
common.parsing_routines.gff_gtf_tools.py
----------------------------------------

This module contains routines for reading, parsing and summarizing data in 
gft and gff files. There are at least three versions of the gff format (so this
will be fun!) and then gtf (v2.2). This is a good page for information about the
versions: http://www.broadinstitute.org/igv/GFF

Note that becuase this is based on gffs and gtfs, its all 1-based.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.12
:created_on: 2013-11-27

==========
change log
==========
 - tweak to allow more flexable reading of go-terms in gtfs (mainly to facilitate the gtfs from Musa)
 - update to support writing annotation regions in bed6 format.
 - update to allow you to rebuild _ordering_array on the fly, if necessary,
   primarily to support indexing using a different 'name' id.
 - update annotation class to make better choices for the _ordering_array name
   column. It now allows you to specify a keyword from the annotation 
   description to use for the region name rather than just use the default. Also
   defaults to id, the name, then gene_id and finally parent if no keyword is
   specified.
 - update annotation class to use 'gene_id' rather than 'name' as it id value
 - added function for stripping leading strings from the file (typically Chr!) 
'''

import numpy, tempfile, time, warnings, copy, re
from parsing_routines.general_classes_and_functions import openfile, region
from parsing_routines.general_classes_and_functions import generic_set_region
from Bio import SeqIO

ver = 1.11

def stripChrStr (fname, stripstr="Chr", verbose=False):
    
    lc = 0
    lrep = 50000
    linc = 0
    retlist=[]
    fin = openfile(fname,"r")
    for line in fin:
        newline = re.sub("^%s" % stripstr,"",line)
        retlist.append(newline)
        lc+=1
        if lc==lrep:
            linc+=1
            lc=0
            if verbose:
                print( "processed %i lines" % int(linc*lrep))
    
    return(retlist)


def parse_bed6(fname, progress_report=False, forceStrand=None,
               generated_by="common.gff_gtf_tools.parse_gff_v3(v%s)" \
                            "" % str(ver),
               preserve_source=False, stripChr=False, splitChar="\t"):
    
    """ read a bed6 file and convert it into a set of 'regions' 
    
    forceStrand overrides any strand information in the bedfile. Acceptable values 
    are '-' or '+'."""
    
    valid_strands=["-","+"]
    
    if stripChr:
        filedata = stripChrStr(fname)
    else:
        filedata = openfile(fname)
    
    headerlines=[]
    type_region_dic={}
    type_region_dic["bedregions"]={}
    
    progress_counter=0
    progress_mult=10000
    progress_done=0
    progress_time=time.time()
    progress_inner_time=time.time()
    
    line_count=0
    region_count = 0
    for line in filedata:
        line = line.strip()
        if line!="":            
            linedata = line.split(splitChar)
            if len(linedata)==6 or (len(linedata)==5 and forceStrand is not None):
                region_name = str(line_count)
                
                # munge attributes
                attrib_dic={}
                
                # add some columns as new attributes
                this_generated_by = generated_by
                if linedata[4]!=".":
                    attrib_dic["score"]=float(linedata[4])
                else:
                    attrib_dic["score"]=None
                
                attrib_dic["id"]=str(linedata[3])
                attrib_dic["name"]=str(linedata[3])               
                
                region_name = attrib_dic["name"]
                
                chrid=None
                if linedata[0]!=".":
                    chrid=linedata[0]
                
                start=None
                if linedata[1]!=".":
                    start=int(linedata[1])

                stop=None
                if linedata[2]!=".":
                    stop=int(linedata[2])
                    
                if forceStrand is None:
                    strand=linedata[5]
                else:
                    if forceStrand in valid_strands:
                        strand = forceStrand
                    else:
                        msg = "Line %i has five fields and forceStrand is " \
                              "either None or not '-' or '+'." % line_count
                        raise IndexError(msg)
                
                this_region = region(region_name, chrid, start, stop, 
                                     strand=strand, desc=attrib_dic,
                                     generated_by=this_generated_by)
                
                # add to data dictionaries                                                
                type_region_dic["bedregions"][region_count]=this_region                    
                region_count+=1
            elif len(linedata)==5 and forceStrand is None:
                msg = "Line %i has five fields and forceStrand is None." \
                      "If this file is missing the strand information, " \
                      "set forceStrand to '+' or '-'. If it is missing " \
                      "other columns, get another file!" % line_count
                raise IndexError(msg)
            else:
                msg = "Line %i does not contain the required fields as per bed6 " \
                      "spec." % line_count
                raise IndexError(msg)
        
        line_count+=1
        progress_counter+=1
        if progress_report and progress_counter==progress_mult:
            progress_done+=1
            print("processed %i new lines in %.4fs, (%i total lines in %.2fs)" \
                  "" % (progress_mult, time.time()-progress_inner_time,
                        progress_mult*progress_done, time.time()-progress_time))
            progress_inner_time=time.time()
            progress_counter=0
        
    return(type_region_dic, headerlines)

def parse_gff_gtf(fname, skip_fasta=True, force_type=None, 
                  progress_report=False, fixID=True, fixName=True,
                  generated_by="common.gff_gtf_tools.parse_gff_v3(v%s)" \
                              "" % str(ver),
                  preserve_source=False, stripChr=False, splitChar="\t"):
    
    """ read a gff or gtf file and convert it into a set of 'regions' 
    
    skip_fasta controls what to do with fasta sequences stored at the end of 
    gff files, force_type controls whether to autodetect the type from the gff
    file or whether to assume a known filetype. fixID forces each regions to 
    include an id. If one doesn't already exist, one will be added for a region
    using either <name>-<int> (if the name tag exists) or, failing that, 
    <parent>-<feature>-<int> or, failing that, just <int>. Similarly,
    fixName forces each region to include a name. If one doesn't already exist
    the ID tag will be used for the name, otherwise <parent>-<feature>-<int> 
    or, failing that, just <int>.
    
    """
    
    if stripChr:
        filedata = stripChrStr(fname)
    else:
        filedata = openfile(fname)
    
    current_filetype=None
    fasta_temp=None
    fasta_isOpen=False
    headerlines=[]
    type_region_dic={}
    type_region_key_dic={}
    
    if force_type is not None:
        current_filetype=force_type
    
    progress_counter=0
    progress_mult=10000
    progress_done=0
    progress_time=time.time()
    progress_inner_time=time.time()
    
    fixcount=1
    line_count=0
    for line in filedata:
        line = line.strip()
        if line!="":
            
            if line.startswith("##gff-version 3"):
                current_filetype="gff3"

            if line.startswith("##gtf"):
                current_filetype="gtf"

            if line.startswith("##FASTA"):
                current_filetype="fasta"
            
            if line.startswith("#"):
                headerlines.append(line)
            elif current_filetype=="gff3" or current_filetype=="gtf":
                linedata = line.split(splitChar)
                if len(linedata)!=9:
                    msg = "Line %i does not contain 9 fields as per gff v3 " \
                          " and  GTF spec." % line_count 
                    raise IndexError(msg)
                else:
                    region_name = str(line_count)
                                        
                    # munge attributes
                    attrib_dic={}
                    if linedata[8]!=".":
                        other_attributes = linedata[8].split(";")
                        for entry in other_attributes:
                            if entry!="":
                                if current_filetype=="gff3":
                                    attrib_data = entry.strip().partition("=")
                                else:
                                    attrib_data = entry.strip().partition(' "')
                                if len(attrib_data)!=3:
                                    print("---",entry,"---")
                                    for x in attrib_data:
                                        print(x)
                                    msg="Error at line %i in the file: " \
                                        "Attribute %s does not have the " \
                                        "correct format: xx=xx,yy" \
                                        "" % (line_count, entry)
                                    raise IndexError(msg)
                                else:
                                    attrib_vals = attrib_data[2].strip('"').split(",")
                                    if len(attrib_vals)==1:
                                        attrib_vals = attrib_vals[0]
                                attrib_dic[attrib_data[0].lower()]=attrib_vals
                    
                    # add some columns as new attributes
                    this_generated_by = generated_by
                    if linedata[1]!=".":
                        attrib_dic["source"]=linedata[1]
                        # interestingly, these two lines here slow things down 
                        # quite a bit
                        if preserve_source: 
                            this_generated_by = linedata[1]
                    if linedata[2]!=".":
                        attrib_dic["type"]=linedata[2]
                    if linedata[5]!=".":
                        attrib_dic["score"]=float(linedata[5])
                    else:
                        attrib_dic["score"]=None
                    if linedata[7]!=".":
                        attrib_dic["phase"]=int(linedata[7])
                    else:
                        attrib_dic["phase"]=None
                    
                    # force ID and name is requested
                    addfixcount=False
                    if "id" not in attrib_dic and fixID:
                        if "name" in attrib_dic:
                            attrib_dic["id"]="%s" % attrib_dic["name"]
                        elif "parent" in attrib_dic:
                            attrib_dic["id"]="%s-%s%i" \
                                             "" % (attrib_dic["parent"],
                                                   attrib_dic["type"],
                                                   fixcount)
                            addfixcount=True
                        else:
                            attrib_dic["id"]=str(fixcount)
                            addfixcount=True
                    
                    if "name" not in attrib_dic and fixName:
                        if "id" in attrib_dic:
                            attrib_dic["name"]=attrib_dic["id"]
                        elif "parent" in attrib_dic:
                            attrib_dic["name"]="%s-%s%i" \
                                               "" % (attrib_dic["parent"],
                                                     attrib_dic["type"],
                                                     fixcount)
                            addfixcount=True
                        else:
                            attrib_dic["name"]=str(fixcount)
                            addfixcount=True
                    
                    if addfixcount:
                        fixcount+=1
                    
                    # make region objects
                    if "name" in attrib_dic:
                        region_name = attrib_dic["name"]
                    
                    chrid=None
                    if linedata[0]!=".":
                        chrid=linedata[0]
                    
                    start=None
                    if linedata[3]!=".":
                        start=int(linedata[3])

                    stop=None
                    if linedata[4]!=".":
                        stop=int(linedata[4])

                    strand=linedata[6]
                    
                    this_region = region(region_name, chrid, start, stop, 
                                         strand=strand, desc=attrib_dic,
                                         generated_by=this_generated_by)
                    
                    # add to data dictionaries
                    type_region_key = line_count
                    if linedata[2]==".":
                        type_region_dic[type_region_key] = this_region
                    else:
                        type_region_key = "%ss" % linedata[2]
                        if type_region_key not in type_region_dic.keys():
                            type_region_dic[type_region_key] = {}
                            type_region_key_dic[type_region_key] = 0
                        else:
                            type_region_key_dic[type_region_key]+=1
                                                    
                        type_region_dic[type_region_key][type_region_key_dic[type_region_key]]=this_region                    
                    
            elif current_filetype=="fasta":
                if not skip_fasta:
                    if not fasta_isOpen:
                        tmpdir = tempfile.mkdtemp()
                        fasta_temp = open("%s/thistempfile.fasta" % tmpdir,"w")
                        fasta_isOpen=True
                    fasta_temp.write("%s\n" % line)
                
        line_count+=1
        progress_counter+=1
        if progress_report and progress_counter==progress_mult:
            progress_done+=1
            print("processed %i new lines in %.4fs, (%i total lines in %.2fs)" \
                  "" % (progress_mult, time.time()-progress_inner_time,
                        progress_mult*progress_done, time.time()-progress_time))
            progress_inner_time=time.time()
            progress_counter=0
        
    if fasta_isOpen:
        fasta_temp.close()
        fasta_temp = open("%s/thistempfile.fasta" % tmpdir, "rU")
        for record in SeqIO.parse(fasta_temp, "fasta") :
            for entry in type_region_dic["chromosomes"].keys():
                if type_region_dic["chromosomes"][entry].name==record.id:
                    type_region_dic["chromosomes"][entry].sequence = record
                    break
    
    return(type_region_dic, headerlines)

class annotation:
    
    """A simple class for representing annotations from a gff or gtf file"""
    
    __sort_order=["chr","start","stop","strand"]
    __ordering_array_dtype=[("type","|U255"), ("id","int"), ("chr","|U255"),
                            ("start","int"),("stop","int"),("strand","|U10"),
                            ("name","|U255")]
    
    __known_types=["gff3","gtf","existing_annotation", "bed6"]
    
    def __init__(self, fname, skip_fasta=True, filetype=None, verbose=False,
                 forceGene_ID=True, forceTranscript_ID=True, generated_by=None,
                 existing_annotation_headerlines=None, preserve_source=False,
                 stripChr=False, featureNameKey=None, forceStrand=None):
        
        """ class constructor: this takes the filename of the gff or gtf file
        and reads the data from it into the class attributes. Filename should
        be a string and point to a file that exists. skip_fasta is passed to 
        the parsing options where appropriate and filetype controls whether to
        assume a particular filetype or to try and autodetect it from the file.
        
        If forceGene_ID or forceTranscript_id is specified then the gene_id
        or transcript_id attribute field will be checked for and, if not
        present, will be filled in by attempting to parse the
        gene<transcript<exon relationships you might expect to find in a gff
        file. This will slow things down a LOT if it has to be done to write a
        valid GTF file.
        
        The forceStrand option only applied to bed6 format and overrides any 
        strand information in the file. 
        
        """
        
        self._full_annotation=None
        self.file_header = None
        self.current_chr=None
        self._featurelist=[]
        self._chrid_dic={}
        self._ordering_array=None
        self.current_feature=None
        self._feature_index=None
        self.current_region=None
        self._region_index=None
        self.current_name=None
        self._name_index=None
        self.current_attributes=None
        self._features_not_found=None
        self._nupdated_geneids=0
        self._nupdated_transcriptids=0
        self._generated_by="common.gff_gtf_tools.annotation(v%s)" % str(ver)
        self.featureNameKey=featureNameKey
        
        if type(generated_by) is str:
            self._generated_by = generated_by
        
        # parse the annotation data
        if filetype in self.__known_types:
            if filetype=="gff3" or filetype=="gtf":
                annotation_data, headerlines = parse_gff_gtf(fname, 
                                                             skip_fasta=skip_fasta,
                                                             force_type=filetype,
                                                             progress_report=verbose,
                                                             generated_by=self._generated_by,
                                                             preserve_source=preserve_source,
                                                             stripChr=stripChr)
            elif filetype=="bed6":
                annotation_data, headerlines =  parse_bed6(fname, progress_report=verbose,
                                                           forceStrand=forceStrand,
                                                           generated_by=self._generated_by,
                                                           stripChr=stripChr)
            elif filetype=="existing_annotation":
                annotation_data = fname
                headerlines = existing_annotation_headerlines
                
            if len(annotation_data)==0:
                raise ValueError("Something went wrong with parsing the file " \
                                 "using filetype %s. Try specifying a " \
                                 "different one of the following with the " \
                                 "argument 'filetype='xxx': %s" \
                                 "" % (filetype, ", ".join(self.__known_types)))
                
        else:
            annotation_data, headerlines = parse_gff_gtf(fname, 
                                                        skip_fasta=skip_fasta,
                                                        progress_report=verbose)
            
            if len(annotation_data)==0:
                raise ValueError("Something went wrong with parsing the file " \
                                 "using autodetection for the filetype. Try " \
                                 "specifying one of the following with the " \
                                 "argument 'filetype='xxx': %s" \
                                 "" % ", ".join(self.__known_types))
        

        if forceGene_ID or forceTranscript_ID:
            # first parse all those regions without 'parent' attributes, give 
            # them gene_ids from their ids if they don't have a gene_id
            # already, and record them.
            gene_pass={}
            for feature_type in annotation_data.keys():
                for feature in annotation_data[feature_type].keys():
                    this_feature = annotation_data[feature_type][feature]
                    if "parent" not in this_feature.desc.keys():
                        if "gene_id" not in this_feature.desc.keys():
                            gene_pass[this_feature.desc["id"]]=this_feature.desc["id"]
                            if forceGene_ID:
                                this_feature.desc["gene_id"]=this_feature.desc["id"]
                                self._nupdated_geneids+=1
                        else:
                            gene_pass[this_feature.desc["id"]]=this_feature.desc["gene_id"]
            
            # then parse the regions for those with 'parent' attributes and,
            # for those with a parent attribute that matches in gene_pass,
            # give them the gene_id of the stored entry, and also assume that 
            # they are a transcript.
            transcript_pass_gene_id={}
            transcript_pass_transcript_id={}
            for feature_type in annotation_data.keys():
                for feature in annotation_data[feature_type].keys():
                    this_feature = annotation_data[feature_type][feature]
                    if "parent" in this_feature.desc.keys():
                        try:
                            parent_gene_id = gene_pass[this_feature.desc["parent"]]
                            if "gene_id" not in this_feature.desc.keys():
                                transcript_pass_gene_id[this_feature.desc["id"]]=this_feature.desc["parent"]
                                if forceGene_ID:
                                    this_feature.desc["gene_id"]=parent_gene_id
                                    self._nupdated_geneids+=1
                            else:
                                transcript_pass_gene_id[this_feature.desc["id"]]=this_feature.desc["gene_id"]
                            if "transcript_id" not in this_feature.desc.keys():
                                transcript_pass_transcript_id[this_feature.desc["id"]]=this_feature.desc["id"]
                                if forceTranscript_ID:
                                    this_feature.desc["transcript_id"]=this_feature.desc["id"]
                                    self._nupdated_transcriptids+=1
                            else:
                                transcript_pass_transcript_id[this_feature.desc["id"]]=this_feature.desc["transcript_id"]
                        except KeyError:
                            pass
                        except TypeError:
                            pass
                        
            # lastly, parse the remaining regions with 'parent' attributes and 
            # if the parent matches in transcript pass, give them the gene_id 
            # of the parent, otherwise give their parent as their gene_id
            for feature_type in annotation_data.keys():
                for feature in annotation_data[feature_type].keys():
                    this_feature = annotation_data[feature_type][feature]
                    if "parent" in this_feature.desc.keys():
                        if forceGene_ID and "gene_id" not in this_feature.desc.keys():
                            try:
                                this_feature.desc["gene_id"]=transcript_pass_gene_id[this_feature.desc["parent"]]
                            except KeyError:
                                this_feature.desc["gene_id"]=this_feature.desc["parent"]
                            except TypeError:
                                this_feature.desc["gene_id"]=",".join(this_feature.desc["parent"])
                            self._nupdated_geneids+=1
                        if forceTranscript_ID and "transcript_id" not in this_feature.desc.keys():
                            try:
                                this_feature.desc["transcript_id"]=transcript_pass_transcript_id[this_feature.desc["parent"]]
                            except KeyError:
                                this_feature.desc["gene_id"]=this_feature.desc["parent"]
                            except TypeError:
                                this_feature.desc["gene_id"]=",".join(this_feature.desc["parent"])
                            self._nupdated_transcriptids+=1
                      
        # store the final annotation
        self._full_annotation = annotation_data
        self.file_header = headerlines
        self._featurelist = self._full_annotation.keys()
        
        #=======================================================================
        # tuplelist=[]
        # for region_type in annotation_data:
        #     self._featurelist.append(region_type)
        #     setattr(self, str(region_type), annotation_data[region_type])
        #     
        #     try:            
        #         for each_region in annotation_data[region_type]:
        #             this_region = annotation_data[region_type][each_region]                    
        #             tuplelist.append(self.__region_to_tuple(this_region,
        #                                                     region_type,
        #                                                     each_region,
        #                                                     specifyNameID=self.featureNameKey))
        #     except TypeError:
        #             this_region = annotation_data[region_type]                    
        #             tuplelist.append(self.__region_to_tuple(this_region,
        #                                                     region_type,
        #                                                     region_type,
        #                                                     specifyNameID=self.featureNameKey))
        # 
        # iarray = numpy.zeros(len(tuplelist),
        #                      dtype=self.__ordering_array_dtype)
        # i=0
        # while i<len(tuplelist):
        #     iarray[i]=tuplelist[i]
        #     i+=1
        #=======================================================================
        
        self._ordering_array = self.__gen_ordering_array()
        
        try:
            for entry in self.chromosomes:
                self._chrid_dic[self.chromosomes[entry].chrid] = entry
        except:
            for entry in numpy.unique(self._ordering_array["chr"]):
                self._chrid_dic[str(entry)] = None        


    def __gen_ordering_array(self):
        
        """generates the main numpy array for working with the annotation """
        
        tuplelist=[]
        for region_type in self._full_annotation.keys():
            try:            
                for each_region in self._full_annotation[region_type]:
                    this_region = self._full_annotation[region_type][each_region]                    
                    tuplelist.append(self.__region_to_tuple(this_region,
                                                            region_type,
                                                            each_region,
                                                            specifyNameID=self.featureNameKey))
            except TypeError:
                    this_region = self._full_annotation[region_type]                    
                    tuplelist.append(self.__region_to_tuple(this_region,
                                                            region_type,
                                                            region_type,
                                                            specifyNameID=self.featureNameKey))
        
        iarray = numpy.zeros(len(tuplelist),
                             dtype=self.__ordering_array_dtype)
        i=0
        while i<len(tuplelist):
            iarray[i]=tuplelist[i]
            i+=1
        
        return(numpy.sort(iarray, order=self.__sort_order))
                
    def __region_to_tuple(self, thisregion, region_type, region_id,
                          specifyNameID=None):
        
        """ build a tuple appropriate for the ordering array from a region 
        
        This should contain 7 elements, see the __ordering_array_dtype"""
        
        thisid=None
        
        if specifyNameID is not None:
            if specifyNameID in thisregion.desc.keys():
                thisid = thisregion.desc[specifyNameID].lower()
        
        if thisid is None:
            thisid=""
            keyorder = ["id", "name", "gene_id", "parent"]
            for key in keyorder:
                if key in thisregion.desc.keys():
                    thisid = thisregion.desc[key].lower()
                    break
                
        this_tuple = (region_type, region_id, thisregion.chrid,
                      thisregion.start, thisregion.stop, 
                      thisregion.strand, thisid)
        
        return(this_tuple)
    
    def set_nameKey(self, newNameKey):
        
        """ sets a new keyword value for the 'name' column for ordering_array"""
        
        self.featureNameKey = newNameKey
        self._ordering_array = self.__gen_ordering_array()
    
    def set_feature(self, feature_str):
        
        """ set the feature type you're interested in """
        
        if type(feature_str) is not str:
            msg = "Feature type should be a string. Available " \
                  "features are: %s" % ", ".join(self._featurelist)
            raise TypeError(msg)
                
        if feature_str not in self._featurelist:
            msg = "Feature type %s not found in annotation. Available " \
                  "features are: %s" % (feature_str,
                                        ", ".join(self._featurelist))
            raise ValueError(msg)
        
        self.current_feature = feature_str
        self._feature_index = numpy.where(self._ordering_array["type"]==self.current_feature)[0]
    
    def clear_feature_selection(self):
        
        """ clears the variables defining the current feature selection """
        
        self.current_feature=None
        self._feature_index=None
    
    def set_region(self, region_rep, start=None, stop=None, strand=None):
        
        """ set the feature region you're interested in. 
        
        Input can either be a 'region' instance or a string that will be parsed 
        into a 'region' instance. See 
        parsing_routines.general_classes_and_functions.py for details on the 
        input for the 'region' constructor.
        """
        
        this_region = generic_set_region(region_rep, start, stop, strand)
                
        # check selected chromosome is in the wigData tracks
        if this_region.chrid not in self._chrid_dic.keys():
            msg = "Invalid region specification, chromosome id not found in " \
                  "annotation. Regions should be specified either as a " \
                  "string 'chrid:start-stop' or with three arguments as " \
                  "set_region('chrid', [start_int], [stop_int])."
            raise ValueError(msg)
                
        # sanitize stupid region start/stop values
        if this_region.start is None:
            this_region.start = 1
        if this_region.stop is None:

            # if there is a chromosome that matches get the max val from there, 
            # other wise get the max val from the max stop of an annotation 
            # on that chr.
            if self._chrid_dic[this_region.chrid] is not None:
                this_region.stop = self.chromosomes[self._chrid_dic[this_region.chrid]].stop
            else:            
                chr_id_ind = numpy.where(self._ordering_array["chr"]==this_region.chrid)[0]
                feature_ends = self._ordering_array["stop"][chr_id_ind]
                this_region.stop = feature_ends.max()+1 # the +1 is important!
        
        self.current_region = this_region
        
        chr_ind = numpy.where(self._ordering_array["chr"]==self.current_region.chrid)[0]
        
        # we want to include features where the start and/or stop occur partway 
        # through the feature. Hence why in indexes are stop>blah and start<blah 
        start_ind = numpy.where(self._ordering_array["stop"]>=self.current_region.start)[0]
        stop_ind = numpy.where(self._ordering_array["start"]<self.current_region.stop)[0]
        
        # if the strandedness is set in the region then use that as an index too
        strand_ind = numpy.arange(len(self._ordering_array["start"]))
        if self.current_region.strand is not None and self.current_region.strand!=".":
            strand_ind = numpy.where(self._ordering_array["strand"]==self.current_region.strand)[0]
        
        the1_ind = numpy.intersect1d(chr_ind, strand_ind)
        the1_ind = numpy.intersect1d(the1_ind, start_ind)
        the1_ind = numpy.intersect1d(the1_ind, stop_ind)
        self._region_index = the1_ind
        
    def clear_region_selection(self):
        
        """ clears the variables defining the current feature selection """
        
        self.current_region=None
        self._region_index=None

    def set_name_selection(self, name_str):
        
        """ set the name(s) of the feature(s) your are interested in """

        if type(name_str) is list or type(name_str) is numpy.ndarray:
            newlist = []
            for name in name_str:
                newlist.append(name.lower())
            # uniqueify
            self.current_name=list(set(newlist))
            mask = numpy.in1d(self._ordering_array["name"],
                              numpy.array(self.current_name))
            self._name_index = numpy.arange(len(self._ordering_array["name"]))[mask]
            if len(self._name_index)!=len(name_str):
                mask = numpy.in1d(numpy.array(self.current_name),
                                  self._ordering_array["name"],
                                  invert=True)
                self._features_not_found = numpy.array(self.current_name)[mask]
        elif type(name_str) is str: 
            self.current_name=name_str.lower()
            self._name_index = numpy.where(self._ordering_array["name"]==self.current_name)[0]
        else:
            msg = "Feature name should be a string, or list of strings"
            raise TypeError(msg)

    def clear_name_selection(self):
        
        """ clears the variables defining the current name selection """
        
        self.current_name=None
        self._name_index=None
        
    def set_attribute_filters(self, this_dict):
        
        """ filters regions for only those with a given attribute value"""
        
        if type(this_dict) is not dict:
            msg = "Attribute filters should be a dictionary of the form: " \
                  "{atribute_name:value, ... }."
            raise TypeError(msg)
        self.current_attributes=this_dict
        
    def clear_attribute_filter(self):
        
        """ clears the variables defining the current feature selection """
        
        self.current_attributes=None
    
    def clear_all(self):
        
        """ clears all the selectors and filter options """
        
        self.current_feature=None
        self._feature_index=None
        self.current_region=None
        self._region_index=None
        self.current_name=None
        self._name_index=None
        self.current_attributes=None
    
    def __build_final_index(self):
        
        """ return a combined index based on feature & region indexes """
        
        indexes = [numpy.arange(len(self._ordering_array["start"]))]
        if self._region_index is not None:
            indexes.append(self._region_index)
        if self._feature_index is not None:
            indexes.append(self._feature_index)
        if self._name_index is not None:
            indexes.append(self._name_index)
        
        if len(indexes)==1 and self.current_attributes is None:
            warnings.warn("No selections specified. Returning all regions!")
        
        final_ind = indexes[0]
        i=1
        while i<len(indexes):
            final_ind = numpy.intersect1d(final_ind, indexes[i])
            i+=1
        return(final_ind)
        
    def get_selection(self, match_case=False, return_indexes=False):
        
        """ return the selected annotations based on feature & region indexes """
        
        final_ind = self.__build_final_index()
        selection = self._ordering_array[final_ind]
        selected_regions=[]
        i=0
        for selected in selection:   
            try:
                this_region = self._full_annotation[selected[0]][selected[1]]
            except IndexError:
                this_region = self._full_annotation[selected[0]]
            
            include_this_region=True
            if self.current_attributes is not None:
                include_this_region=False
                for key in self.current_attributes:
                    if key in this_region.desc.keys():
                        if match_case and self.current_attributes[key]==this_region.desc[key]:
                            include_this_region=True
                        elif not match_case and self.current_attributes[key].lower()==this_region.desc[key].lower():
                            include_this_region=True
            
            if include_this_region:
                if return_indexes:
                    selected_regions.append(final_ind[i])
                else:
                    selected_regions.append(this_region)
            i+=1
        
        if len(selected_regions)==0:
            raise ValueError("No regions selected. Nothing to get!")
        
        return(selected_regions)
    
    def get_selection_names(self, index=None):

        """ return the names of the selected annotations 
        
        The keyword arg 'index' allows a specific index to be supplied. for the
        chrid position selection. This should be used mostly internally ;)
        """
        
        if index is None:
            final_ind = self.get_selection(return_indexes=True)
        else:
            final_ind = index
        
        names = self._ordering_array["name"][final_ind]
        return(names)

    def get_selection_attribute(self, attribute):

        """ return the keyword values of the selected annotations 
        
        The attribute here is a property tagged to a feature via the 
        description column in the annotation.
        """
        
        sel = self.get_selection()
        keywordvals = []
        for region in sel:
            if attribute in region.desc:
                keywordvals.append(region.desc[attribute])
            else:
                keywordvals.append("")
        
        return(numpy.array(keywordvals))

    def get_selection_types(self, index=None):

        """ return the types of the selected annotations 
        
        The keyword arg 'index' allows a specific index to be supplied. for the
        chrid position selection. This should be used mostly internally ;)
        """
        
        if index is None:
            final_ind = self.get_selection(return_indexes=True)
        else:
            final_ind = index
        
        types = self._ordering_array["type"][final_ind]
        return(types)

    def get_selection_chrids(self, index=None):

        """ return the chrids of the selected annotations 
        
        The keyword arg 'index' allows a specific index to be supplied. for the
        chrid position selection. This should be used mostly internally ;)
        """
        
        if index is None:
            final_ind = self.get_selection(return_indexes=True)
        else:
            final_ind = index
        
        chrids = self._ordering_array["chr"][final_ind]
        return(chrids)

    def get_selection_starts(self, index=None):

        """ return the start positions of the selected annotations 
        
        The keyword arg 'index' allows a specific index to be supplied. for the
        start position selection. This should be used mostly internally ;)
        """
        
        if index is None:
            final_ind = self.get_selection(return_indexes=True)
        else:
            final_ind = index
        
        starts = self._ordering_array["start"][final_ind]
        return(starts)

    def get_selection_stops(self, index=None):

        """ return the stop positions of the selected annotations

        The keyword arg 'index' allows a specific index to be supplied. for the
        stop position selection. This should be used mostly internally ;)
        """
        
        
        if index is None:
            final_ind = self.get_selection(return_indexes=True)
        else:
            final_ind = index

        stops = self._ordering_array["stop"][final_ind]
        return(stops)

    def get_selection_strands(self, index=None):

        """ return the strands of the selected annotations 
        
        The keyword arg 'index' allows a specific index to be supplied. for the
        chrid position selection. This should be used mostly internally ;)
        """
        
        if index is None:
            final_ind = self.get_selection(return_indexes=True)
        else:
            final_ind = index
        
        strands = self._ordering_array["strand"][final_ind]
        return(strands)
    
    def get_selection_lengths(self):

        """ return the lengths of the selected annotations """
        
        final_ind = self.get_selection(return_indexes=True)
        
        return(1+(self.get_selection_stops(final_ind)-self.get_selection_starts(final_ind)))
    
    def get_selection_sum_length(self):
        
        """ return the sum of the lengths of the selected annotations """
        
        return(self.get_selection_lengths().sum())

    def get_selection_mean_length(self):
        
        """ return the mean of the lengths of the selected annotations """
        
        return(self.get_selection_lengths().mean())
    
    def get_selection_median_length(self):
    
        """ return the median of the lengths of the selected annotations """
        
        return(numpy.median(self.get_selection_lengths()))

    def get_selection_stddev_length(self):
    
        """ return the unbiassed standard deviation of the lengths of the 
        selected annotations. 
        
        Here we correct numpys typical standard deviation routine to make it an
        unbiased estimator, according to GURLAND, J. & TRIPATHI, R. C. 1971. 
        "Simple Approximation for Unbiased Estimation of Standard Deviation." 
        American Statistician, 25, 30. Note that this correction will propogate 
        through to the standard error on the mean. 

        """
        
        return(self.get_selection_lengths().std(ddof=1)*
               (1.0+(1.0/(4.0*(len(self.get_selection_lengths())-1.0)))))

    def get_selection_stderr_length(self):

        """ return the unbiassed standard error on the mean of the lengths of 
        the selected annotations."""
                
        return(self.get_selection_stddev_length()/len(self.get_selection_lengths()))
        
    def get_selection_lengthstats(self):

        """ return the lengths of the selected annotations """

        return(self.get_selection_mean_lengths(),
               self.get_selection_median_lengths(),
               self.get_selection_stddev_lengths(),
               self.get_selection_stderr_lengths())
    
    def get_features(self):
        
        """ get the list of features found in the annotation """
        
        return(self._featurelist)
    
    def get_chromosomes(self):
        
        """ get the list of chromosomes in the whole annotation """
        
        return(sorted(self._chrid_dic.keys()))

    def get_gene_length_bp(self, feature="exons", inc_UTRs=True,
                           ret_type="dict"):
        
        """ return the bp length of all exons and UTRs in each gene 
        
        Note that region, name and attribute selections are preserved
        """

        allowed_ret_types=["dict","structured_array"]
        
        if ret_type not in allowed_ret_types:
            warnings.warn("Unknown return type requested. Returning a " \
                          "dictionary")
            ret_type = "dict"
        
        # change current feature to chosen one
        if self.current_feature!=feature:
            warnings.warn("Currently set feature type (%s) doesn't match the " \
                          "requested feature (%s). The set feature type will " \
                          "be changed to the requested type for this call, " \
                          "but will then be reset." % (self.current_feature,
                                                       feature))
        
        oldfeature = self.current_feature
        
        this_feature = None
        if feature=="exons" or feature=="CDS":
            this_feature=feature
        else:
            raise ValueError("please choose either 'exons' or 'CDS' as " \
                             "the feature type.")
        
        self.set_feature(this_feature)
        exons = self.get_selection()
        geneiddic={}
        for exon in exons:
            try:
                geneiddic[exon.desc["gene_id"]]+=exon.get_length()
            except KeyError:
                geneiddic[exon.desc["gene_id"]]=exon.get_length()
        
        if feature=="CDS" and inc_UTRs:
            for utytype in ["five_prime_UTRs","three_prime_UTRs"]:
                self.set_feature(utytype)
                utrs = self.get_selection()
                for utr in utrs:
                    try:
                        geneiddic[utr.desc["gene_id"]]+=utr.get_length()
                    except KeyError:
                        geneiddic[utr.desc["gene_id"]]=utr.get_length()

        # return old feature choice to memory
        self.current_feature = oldfeature
        
        if ret_type=="structured_array":
            newarr = numpy.zeros(len(geneiddic.keys()),
                                 dtype=[("Geneid","|S255"),("length (bp)",int)])
            
            i=0
            for geneid in geneiddic.keys():
                newarr[i] = (geneid,geneiddic[geneid])
                i+=1
            
            return(newarr)
        
        elif ret_type=="dict":
            return(geneiddic)
    
    def display(self, fmt="gff3", limit=None, printthem=True):
        
        regions_to_display = self.get_selection()
        
        if limit is None or len(regions_to_display)<limit:
            limit = len(regions_to_display)
        
        
        printstrs=[]
        i=0
        while i<limit:
            if fmt=="gff3":
                printstrs.append(regions_to_display[i].get_gff3line())
            elif fmt=="gtf":
                printstrs.append(regions_to_display[i].get_gtfline())
            i+=1
        
        if limit<len(regions_to_display):
            printstrs.append("%i (of %i) regions shown." \
                             "" % (limit,len(regions_to_display)))
        
        if printthem:
            print("\n".join(printstrs))
        else:
            return("\n".join(printstrs))
        
    def write(self, filename, fmt="gff3", regions_to_write=None, append=False,
              compress=True, verbose=False):
         
        """ write the existing annotation, inc. selections, to a file 
        
        The fmt option will determine the file format. Options are 'gff3', 
        'gtf' or bed6
        """
        
        if regions_to_write is None:
            regions_to_write = self.get_selection()
        
        if not append:
            fh = openfile(filename, "w", compress=compress)
        else:
            fh = openfile(filename, "a", compress=compress)
        
        if fmt=="gff3":
            fh.write("##gff-version 3\n")
        
        for region in regions_to_write:
            if fmt=="gff3":
                fh.write("%s\n" % region.get_gff3line())
            elif fmt=="gtf":
                fh.write("%s\n" % region.get_gtfline(verbose=verbose))
            elif fmt=="bed6":
                fh.write("%s\n" % region.get_bed6line())
            else:
                raise ValueError("Unknown format option: fmt=%s. Please " \
                                 "specify either 'gff3', 'gtf' or 'bed6'.")
        
        fh.close()
        
        return(filename)
     
    def __str__(self):
        return(self.display(limit=10))

    def __repr__(self):
        return(self.display(limit=10, printthem=False))
            
def mergeAnnotations(input_filenames, logger, filetype="gff3",
					 overwrite_strand_list=None, generated_by=None):
    """ merges two or more annotations, optionally overwriting strnad info """
    
    new_annotation=None
    annotations=[]
    i=0
    while i<len(input_filenames):
        filename = input_filenames[i]        
        logger.info("Reading annotation from %s..." % filename)
        annot = annotation(filename, filetype=filetype,
                           generated_by=filename)
        annot.clear_all()
        logger.info("Found %i annotations..." \
                           "" % len(annot.get_selection()))
        
        annot_full = annot._full_annotation
        # Overwrite strand information if specified.        
        if overwrite_strand_list is not None:
            overwrite_strand_option = overwrite_strand_list[i]
            if overwrite_strand_option != "/":
                annot_full = overwrite_strand(annot_full,
                                              overwrite_strand_option,
                                              logger)
        
        if new_annotation is None:
            new_annotation = annot_full
        else:
            annotations.append(annot_full)
        
        i+=1
    
    for annot in annotations:
        for key in annot.keys():
            if key not in new_annotation.keys():
                new_annotation[key]=annot[key]
            else:
                current_val = max(new_annotation[key].keys())+1
                for entry in annot[key].keys():
                    new_annotation[key][current_val]=annot[key][entry]
                    current_val+=1
    
    full_new_annot = annotation(new_annotation, filetype="existing_annotation",
								generated_by=generated_by)
    return(full_new_annot)
             
def overwrite_strand(full_annotation, newstrand, logger):
    
    """ takes a full annotation region and overwrites the strand information
    
    This overwrite *all* strand information. Use carefully!"""
    
    if type(full_annotation) is not dict:
        raise TypeError("Please supply a full annotation dictionary in order " \
                        "to overwrite the strand information. This can be " \
                        "obtained from the annotation class by accessing " \
                        "annotation._full_annotation.")
    
    
    if type(newstrand) is not str:
        raise TypeError("New strand information to be specified must be a " \
                        "string. Allowed options are: '.','+','-'")
    
    strand_choices = [".","+","-"]
    if newstrand not in strand_choices:
        raise ValueError("New strand information %s not recognized. Allowed " \
                         "options are: '.','+','-'" % newstrand)
    
    logger.info("overwritting strand info with %s ..." % newstrand)
    
    count = 0
    for key in full_annotation.keys():
        for region in full_annotation[key]:
            full_annotation[key][region].strand = newstrand
            count +=1
    
    logger.info("\tupdated %i regions." % count)
    return(full_annotation)
            
         
