'''
--------------------------------------------------------
common.parsing_routines.generic_classes_and_functions.py
--------------------------------------------------------

This module contains generic classes used in the other parsing routines. 

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.11
:created_on: 2013-11-26

##### Change Log #####
1.11 added bed6 option for regions.
1.10 added the functionality to autogenerate n leading columns from the json
1.9 fixed some problems with get_gff3line and get_gtfline when desc is None, 
    and added generated_by passthrough.
1.8 updated missing_value default to None for genfromtext.

'''

ver=1.10

import re, os, numpy, copy, warnings, json
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna, generic_protein

def openfile(fname, openMode="r", compress=False, readfile=True,
             skipComments=False, commentChar="#", debugLineLimit=0):
    
    """ open and optionally read a gzipped or uncompressed file 
    
    If the file is read, return the read data. If the file is just opened, 
    return the filehandle."""
    
    import gzip
    
    # try and cast string to a string if it isn't already!
    if type(fname) is not str:
        raise TypeError("The filename passed to openfile is not a " \
                        "string. Type is %s" % type(fname))
    
    # open the file
    if openMode.startswith("r"):

        # resolve relative paths to a full path
        full_path = os.path.abspath(fname)
        
        # check that the file exists
        if not os.path.isfile(full_path):
            msg = "The file specified (%s) does not exist.\nPlease specify a " \
                  "valid file." % full_path
            raise IOError(msg)
        
        try:
            filehandle = gzip.open(fname, openMode)
            filehandle.readline()
            filehandle.rewind()
        except IOError:
            filehandle = open(fname, openMode)
                
        filedata=None
        if readfile:
            
            # read all data for speed
            filedata = filehandle.readlines()
            filehandle.close()
            
            # slice for debug is enabled:
            if debugLineLimit>0:
                filedata = filedata[0:len(filedata)-debugLineLimit]
            
            # skip comments if indicated
            returndata = filedata
            if skipComments:
                returndata=[]
                for line in filedata:
                    if not line.startswith(commentChar):
                        returndata.append(line)
                        
            return(returndata)
        else:
            return(filehandle)
    
    elif compress:
        filehandle = gzip.open(fname, openMode)
        return(filehandle)
    else:
        filehandle = open(fname, openMode)
        return(filehandle)

class region:
    
    """A simple class for representing a a chromosome region"""
    
    def __init__(self, region_name, chrid, start, stop,
                 strand=None, seq=None, seqtype=None, desc=None,
                 generated_by=None):

        """ class contrustor: this requires name and chromosome id as strings, 
        start and stop and ints, and optionally strand, sequence and 
        description as strings"""
                
        self.name=None
        self.chrid=None
        self.start=None
        self.stop=None
        self.strand=None
        self.sequence=None
        self.desc=None
        self._generated_by_script="common.general_classes_and_functions(%s)" \
                                  "" % str(ver)

        if type(region_name) is str and region_name!="":
            self.name=region_name
        else:
            raise TypeError("Region name must be a non-empty string.")
        
        if type(chrid) is str and chrid!="":
            self.chrid=chrid
        else:
            raise TypeError("Chromosome ID must be a non-empty string.")
        
        if type(start) is int or type(start) is numpy.int64 or start is None:
            self.start=start
        else:
            raise TypeError("Region start must be an integer (or None).")
        
        if type(stop) is int or type(stop) is numpy.int64 or stop is None:
            self.stop=stop
        else:
            raise TypeError("Region stop must be an integer (or None).")
        
        if strand is not None and type(strand) is str:
            if strand=="+" or strand=="-" or strand==".":
                self.strand=strand
            elif strand=="fwd" or strand=="forward":
                self.strand="+"
            elif strand=="rev" or strand=="reverse":
                self.strand="-"
            elif strand=="nostrand":
                self.strand="."
            else:
                raise ValueError("Strand specification should be one of " \
                                 "+/-/./fwd/rev/forward/reverse/nostrand.")
        elif strand is not None and type(strand) is not str:
            raise TypeError("Region stand must be a string.")
        elif strand is None:
            self.strand="."

        if seq is not None:
            if type(seq) is Seq:
                self.sequence=seq
            elif type(seq) is str and seqtype is None:
                self.sequence=Seq(seq)
            elif type(seq) is str and type(seqtype) is str:
                if seqtype=="dna" or seqtype=="DNA":
                    self.sequence=Seq(seq, generic_dna)
                elif seqtype=="rna" or seqtype=="RNA":
                    self.sequence=Seq(seq, generic_rna)
                elif seqtype=="prot" or seqtype=="protein":
                    self.sequence=Seq(seq, generic_protein)
                else:
                    raise ValueError("type of sequence must be one of " \
                                     "dna/DNA/rna/RNA/prot/protein.") 
            elif type(seq) is str and type(seqtype) is not str and seqtype is not None:
                raise TypeError("Region sequence 'type' must be a string or " \
                                "None.")            
            elif type(seq) is not str and type(seq) is not Seq:
                raise TypeError("Region sequence must be a string or a Seq " \
                                "record from Biopython or None.")
        
        if desc is not None and (type(desc) is str or type(desc) is dict):
            self.desc=desc
        elif desc is not None:
            raise TypeError("Region desc must be a string or a dictionary.")
        
        if type(generated_by) is str:
            self._generated_by_script=generated_by
    
    def get_length(self):
        
        """ get the length (start-stop) of the feature """
        
        return(1+(self.stop-self.start))
    
    def get_positions_str(self):
        
        """ get the chr:start-stop of the feature """
        
        pos_str = "%s:%s-%s" % (str(self.chrid), str(self.start), str(self.stop))
        
        return(pos_str)
    
    def get_desc_str(self):
        
        """ get the desc of the feature as a csv string"""
        
        desc_str_list=[]
        for key in self.desc.keys():
            desc_str = u"%s:%s" % (key, self.desc[key])
            try:
                desc_str_list.append(desc_str.encode('ascii', 'ignore'))
            except:
                print key
                print self.desc[key]
                raise
        
        final_desc_str = ", ".join(desc_str_list)
        
        return(final_desc_str)
    
    def get_gff3line(self):
        
        """ gets the current region as a gff3 format line.
        
        The format is specified here http://www.sequenceontology.org/gff3.shtml
        """
                
        redundant_keys = ["score", "phase", "type"]
        replace_keys = [("source", "original_source"), ("id","ID")]
        
        this_attributes={}
        if self.desc is not None:
            this_attributes = copy.deepcopy(self.desc)
        
        this_type = "unknown"
        if "type" in this_attributes:
            if this_attributes["type"] is not None:
                this_type = str(this_attributes["type"])
        
        this_score = "."
        if "score" in this_attributes:
            if this_attributes["score"] is not None:
                this_score = str(this_attributes["score"])

        this_phase = "."
        if "phase" in this_attributes:
            if this_attributes["phase"] is not None:
                this_phase = str(this_attributes["phase"])
        
        if self.strand is None:
            this_strand = "."
        else:
            this_strand=self.strand
        
        for key in redundant_keys:
            if key in this_attributes.keys(): del this_attributes[key]
        for key in replace_keys:
            if key[0] in this_attributes.keys():
                this_attributes[key[1]]=this_attributes[key[0]]
                del this_attributes[key[0]]
        
        attributes_list = []
        for key in this_attributes:
            keyval = key.capitalize()
            if key=="ID":
                keyval="ID"
            
            attributes_list.append('%s=%s' % (keyval,
                                              str(this_attributes[key]))
                                   )
        
        vals = [self.chrid, self._generated_by_script, this_type, 
                str(self.start), str(self.stop), this_score,
                this_strand, this_phase, ';'.join(attributes_list)]
                
        gff3line = "\t".join(vals)
        
        return(gff3line)

    def get_gtfline(self, gene_id=None, transcript_id=None, verbose=False):
        
        """ gets the current region as a gtf2 format line.
        
        http://www.ensembl.org/info/website/upload/gff.html
        """
        
        redundant_keys = ["score", "phase", "type"]
        replace_keys = [("source", "original_source"), ("id","ID")]
        
        this_attributes={}
        if self.desc is not None:
            this_attributes = copy.deepcopy(self.desc)
        
        this_type = "unknown"
        if "type" in this_attributes:
            if this_attributes["type"] is not None:
                this_type = str(this_attributes["type"])
        
        this_score = "."
        if "score" in this_attributes:
            if this_attributes["score"] is not None:
                this_score = str(this_attributes["score"])

        this_phase = "."
        if "phase" in this_attributes:
            if this_attributes["phase"] is not None:
                this_phase = str(this_attributes["phase"])
        
        if self.strand is None:
            this_strand = "."
        else:
            this_strand=self.strand
        
        if gene_id is not None:
            this_attributes["gene_id"]=gene_id
        
        if transcript_id is not None:
            this_attributes["transcript_id"]=transcript_id
        
        if "gene_id" not in this_attributes.keys():
            if "ID" in this_attributes.keys():
                name = this_attributes["ID"]
            elif "name" in this_attributes.keys():
                name = this_attributes["name"]
            else:
                name = "\t".join([self.chrid, this_type, 
                                  str(self.start), str(self.stop)])
            if this_type is not "gene" and verbose:
                warnings.warn("Annotation does not contain a gene_id " \
                              "attribute for gene %s. Potentially writing " \
                              "invalid GTF" % name)

        if "transcript_id" not in this_attributes.keys():
            if "ID" in this_attributes.keys():
                name = this_attributes["ID"]
            elif "name" in this_attributes.keys():
                name = this_attributes["name"]
            else:
                name = "\t".join([self.chrid, this_type, 
                                  str(self.start), str(self.stop)])
            if verbose:
                warnings.warn("Annotation does not contain a transcript_id " \
                              "attribute for gene %s. Potentially writing " \
                              "invalid GTF" % name)

        for key in redundant_keys:
            if key in this_attributes.keys(): del this_attributes[key]
        for key in replace_keys:
            if key[0] in this_attributes.keys():
                this_attributes[key[1]]=this_attributes[key[0]]
                del this_attributes[key[0]]
        
        attributes_list = []
        for key in this_attributes:
            keyval = key.lower()
            if key=="ID":
                keyval="ID"
            
            attributes_list.append('%s "%s"' % (keyval,
                                                str(this_attributes[key]))
                                   )
        
        vals = [self.chrid, self._generated_by_script, this_type, 
                str(self.start), str(self.stop), this_score,
                this_strand, this_phase, '; '.join(attributes_list)]
                
        gtfline = "\t".join(vals)
        
        return(gtfline)

    def get_bed6line(self):
        
        """ gets the current region as a bed6 format line.
        
        The format is specified here http://https://genome.ucsc.edu/FAQ/FAQformat.html
        """
                
        this_attributes={}
        if self.desc is not None:
            this_attributes = copy.deepcopy(self.desc)
                
        this_score = "."
        if "score" in this_attributes:
            if this_attributes["score"] is not None:
                this_score = str(this_attributes["score"])
                
        if self.strand is None:
            this_strand = "."
        else:
            this_strand=self.strand
        
        this_name = "None"
        if self.name is None:
            if "id" in this_attributes:
                if this_attributes["id"] is not None:
                    this_name = str(this_attributes["name"])
            else:
                this_name="None"
        else:
            this_name = self.name
                
        vals = [self.chrid, str(self.start), str(self.stop), this_name, this_score,
                this_strand]
                
        bed6line = "\t".join(vals)
        
        return(bed6line)
    
    def __str__(self):
        return(self.get_gff3line())

    def __repr__(self):
        return(self.get_gff3line())

def generic_set_region(region_rep, start=None, stop=None, strand=None,
                       generated_by=None):
    
    ''' Returns the input region, or builds one from the input details
    
    This is used in several of the parsers including the annotation, and wigData
    classes
    '''

    general_msg = "Invalid region specification. Regions should be either a " \
                  "valid 'region' instance, or be specified as a string of " \
                  "format 'chrid' (which results in selecting the entire " \
                  "chromosome) or 'chrid:start-stop' or with three arguments " \
                  "as set_region('chrid', start_int, stop_int)."
        
    if region_rep.__class__.__name__ is 'region':
        return(region_rep)
    
    if type(region_rep) is not str:
        try:
            region_rep = str(region_rep)
        except:
            raise ValueError("Unable to cast region identifier to a string")
    
    if type(region_rep) is str and region_rep!="":
        if start is not None and stop is None:
            raise ValueError("Please specify a stop position")
        if stop is not None and start is None:
            raise ValueError("Please specify a start position")
        if (start is not None and type(start) is not int and 
            type(start) is not numpy.int64):
            raise TypeError("start position should be an int")
        if (stop is not None and type(stop) is not int and 
            type(stop) is not numpy.int64):
            raise TypeError("stop position should be an int")
        
        # parse string for region details if necessary
        if start is None and stop is None and strand is None:
            try:
                split_1 = re.split(":",region_rep.strip())
                chrid = split_1[0]
            except IndexError:                
                raise ValueError(general_msg)
            
            if len(split_1)>1:
                try:
                    split_2 = re.split("-",split_1[1])
                    start = int(split_2[0])
                    stop = int(split_2[1])
                except IndexError:
                    raise ValueError(general_msg)            
        else:
            if start is not None or stop is not None:
                if type(start) is not int or type(stop) is not int:
                    raise ValueError(general_msg)
            chrid = region_rep
        
        if start is not None and stop is not None:
            # sanitize stupid region start/stop values
            if start<1:
                msg="Attempt to set start position less-than 1."
                raise ValueError(msg)
            if stop<1:
                msg="Attempt to set stop position less-than 1."
                raise ValueError(msg)
            if start==stop:
                msg="Attempt to set zero-length region."
                raise ValueError(msg)
            if start>stop:
                msg="Attempt to set stop position < start position."
                raise ValueError(msg)
        return(region("set_region", chrid, start, stop, strand,
                      generated_by = generated_by))
    elif region_rep=="":
        raise ValueError(general_msg)
    else:
        if type(region_rep) is not str:
            msg="chromosome id must be a string, and isn't. Try again...."
        else:
            raise TypeError(general_msg)

def makeStructuredArray(data, dtype, delimiter="\t"):
    
    ''' makes a structured array from delimited data with a column mapping '''
    this_array = numpy.zeros(len(data), dtype=dtype)
    rec_count=0
    i=0
    while i<len(data):
        line = data[i]
        linedata = line.strip().split(delimiter)
        if len(linedata)>1:
            linedata[0] = linedata[0].replace('"','')
            this_array[i]=tuple(linedata)
            rec_count+=1            
        i+=1
    
    if rec_count<len(data):
        this_array=this_array[:rec_count]
    
    return(this_array)        

def getDataFromFormat(filename, dataformat, available_formats, 
                      skip_comments=True, logger=None, verbose=False, slen=100):
    
    """reads a file in a particular format and gets the file data """
    
    # type checks
    
    if type(filename) is not str:
        raise TypeError("The specified filename (%s) is not a string object." \
                        "" % filename)
        
    if type(dataformat) is not str:
        raise TypeError("The specified file format (%s) is not a string " \
                        "object." % dataformat)
    
    if type(available_formats) is not dict and type(available_formats) is not str:    
        raise TypeError("The available formats are not in a dictionary, nor " \
                        "are they a string filename to a json format config " \
                        "file. Please provide a formats dictionary, or the " \
                        "path to a valid json format config file with " \
                        "formats in." )

    if type(available_formats) is str:
        try:
            available_formats = json.load(openfile(available_formats,
                                                   openMode="r",
                                                   readfile=False))
        except:
            print("The available formats appears to be a string filename to " \
                  "a json format config file but I can't read/parse it " \
                  "correctly, please try again!")
            raise 
    
    if verbose and logger is not None:
        logger.info("reading format information...")
    
    # identify the format
    thisformat=None
    try:
        thisformat = available_formats[dataformat]
    except:
        raise ValueError("Format %s not found in the formats dictionary " \
                         "file specified by %s." \
                         "" % (dataformat,
                               available_formats))
    
    # add user defined string lengths
    if "string_field_length" in thisformat:
        slen = int(thisformat["string_field_length"])
        
    names=[]
    fieldformats=[]
    
    # auto generated columns
    if "gencols" in thisformat.keys():
        if verbose and logger is not None: 
            logger.info("attempting to auto-generate columns...")
        try:
            input_filehandle = openfile(filename, readfile=False)
            lno=1
            while lno<thisformat["gencols_names_line"]+1:
                headnames = input_filehandle.readline()
                lno+=1
            headnames = numpy.array(headnames.split(thisformat["separator"]),
                                    dtype="str")
            gencol_names = headnames[thisformat["gencols"]]
            for name in gencol_names:
                names.append(str(name))
                if thisformat["gencols_type"]=="str":
                    fieldformats.append("S%s" % slen)
                else:
                    fieldformats.append(thisformat["gencols_type"])
            input_filehandle.close()
        except KeyError:
            if verbose and logger is not None:
                logger.info("something went wrong with column auto-generation")
            raise
    
    # pre-defined_columns
    if verbose and logger is not None:
        logger.info("processing pre-defined columns...")
    for col in thisformat["cols"]:
        # make sure to convert to str to avoid unicode issues
        names.append(str(col["name"]))
        if col["type"]=="str":
            fieldformats.append("S%s" % slen)
        else:
            fieldformats.append(col["type"])
    
    # set dtype
    dtype = {'names':names, 'formats':fieldformats}
    
    if verbose and logger is not None:
        logger.info("reading data (skipping the first %i lines)..." \
                    "" % thisformat["start_skip_lines"])
    
    comments="#" # this is the numpy default for loadtxt
    if skip_comments and "comment_char" in thisformat.keys():
        comments = thisformat["comment_char"]
    
    skiplines = 0
    if "start_skip_lines" in thisformat.keys():
        skiplines = thisformat["start_skip_lines"]
    
    skipfoot = 0
    if "start_skip_lines" in thisformat.keys():
        skipfoot = thisformat["end_skip_lines"]
    
    missing_values = None
    if "missing_value_strings" in thisformat.keys():
        missing_values = thisformat["missing_value_strings"]
        
    input_filehandle = openfile(filename, readfile=False)
    if "separator" in thisformat.keys():
#        data = numpy.loadtxt(input_filehandle, dtype=dtype, comments=comments,
#                             skiprows=skiplines,
#                             delimiter=thisformat["separator"])
        data = numpy.genfromtxt(input_filehandle, dtype=dtype, comments=comments,
                             skip_header=skiplines, skip_footer=skipfoot, 
                             autostrip=True, delimiter=thisformat["separator"],
                             missing_values=missing_values)
    else:
        data = numpy.genfromtxt(input_filehandle, dtype=dtype, comments=comments,
                             skip_header=skiplines, skip_footer=skipfoot, 
                             autostrip=True, missing_values=missing_values)
    
    if "sort_field" in thisformat.keys():
        data = numpy.sort(data, order=[thisformat["sort_field"]])[::-1]
    
    if verbose and logger is not None:
        logger.info("trimming the last %i lines)..." \
                    "" % thisformat["end_skip_lines"])
    
    if int(thisformat["end_skip_lines"])!=0:
        data = data[0:(-1*int(thisformat["end_skip_lines"]))]
    
    return(data, thisformat)

def addStructuredArrayfield(a, descr, logger=None, verbose=False):
    
    """Return a new array that is like "a", but has additional fields.

    Arguments:
      a     -- a structured numpy array
      descr -- a numpy type description of the new fields

    The contents of "a" are copied over to the appropriate fields in
    the new array, whereas the new fields are uninitialized.  The
    arguments are not modified.
    
    from:http://stackoverflow.com/questions/1201817/adding-a-field-to-a-structured-numpy-array
    
    >>> sa = numpy.array([(1, 'Foo'), (2, 'Bar')], \
                         dtype=[('id', int), ('name', 'S3')])
    >>> sa.dtype.descr == numpy.dtype([('id', int), ('name', 'S3')])
    True
    >>> sb = add_field(sa, [('score', float)])
    >>> sb.dtype.descr == numpy.dtype([('id', int), ('name', 'S3'), \
                                       ('score', float)])
    True
    >>> numpy.all(sa['id'] == sb['id'])
    True
    >>> numpy.all(sa['name'] == sb['name'])
    True
    """
    
    if a.dtype.fields is None:
        raise ValueError, "`A' must be a structured numpy array"
    
    if verbose and logger is not None:
        logger.info("adding new field to the structured data array")
    
    newdtype=[]
    for val in a.dtype.descr:
        newtup = (str(val[0]),val[1])
        newdtype.append(newtup)
    
    newdtype = newdtype+descr
    
    b = numpy.empty(a.shape, dtype=newdtype)

    for name in a.dtype.names:
        b[name] = a[name]
    
    return b

def computeIntrons(exon_regions, logger=None, generated_by=None, verbose=False):
    
    """ takes an set of exon regions and computes the corresponding introns
    
    Basically, place each exon in order and create the intron spaces, labelling
    the appropriately. All generic labels common to the exon set will be applied
    to the introns. In addition, the feature type will be intron and the introns
    will be named appropriately"""
    
    if verbose:
        logger.info("Identifying exon common description items...")
        
    exon_common_desc_keys=None
    for exon in exon_regions:
        for desckey in exon.desc.keys():
            if exon_common_desc_keys is None:
                exon_common_desc_keys=exon.desc
            elif desckey in exon_common_desc_keys:
                if exon_common_desc_keys[desckey]!=exon.desc[desckey]:
                    del exon_common_desc_keys[desckey]
    
    if verbose:
        logger.info("Computing intron regions...")
    
    introns = []
    i=0
    while i<len(exon_regions)-1:
        this_region_name = "%s.intron%i" % (exon_regions[i].desc["parent"], i+1)
        this_region_desc = copy.deepcopy(exon_common_desc_keys)
        this_region_desc["id"]=this_region_name
        this_region_desc["name"]=this_region_name
        this_region_desc["type"]="intron"
        this_region_desc["intron_index"]=i+1
        
        this_intron = region(this_region_name, exon_regions[i].chrid,
                             exon_regions[i].stop+1, exon_regions[i+1].start-1,
                             strand=exon_regions[i].strand,
                             desc=this_region_desc, generated_by=generated_by)
        introns.append(this_intron)
        i+=1
    
    return(introns)
    
    
