'''
------------------------------------
common.parsing_routines.wig_tools.py
------------------------------------

This module contains routines for reading, parsing and summarizing data in 
wig files. Python 3 conversion.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.6
:created_on: 2013-11-22
'''

from parsing_routines.general_classes_and_functions import openfile
from parsing_routines.general_classes_and_functions import generic_set_region
import re, numpy, os, tempfile, subprocess, warnings, sys, copy
ver=1.6

def convert_BigWig2Wig(fname, path_to_binary=None):
    
    """ uses the UCSC bigwigToWig to do the conversion, ready for reading """
        
    tmp_path = tempfile.mkdtemp()
    success=subprocess.check_call(["BigWig2Wig", fname,
                                   "%s/out.wig" % tmp_path])
    
    if success==0:
        return("%s/out.wig" % tmp_path)

def convert_Wig2BigWig(fname, chromSizes, outfile, verbose=False):
    
    """ uses the UCSC wigToBigWig to do the conversion """
    
    # try and cast string to a string if it isn't already!
    if type(fname) is not str:
        raise TypeError("The filename passed to convert_Wig2BigWig is not a " \
                        "string. Value: %s, type: %s" \
                        "" % (str(fname), type(fname)))
    
    # resolve relative paths to a full path
    full_path = os.path.abspath(fname)
    
    # check that the file exists
    if not os.path.isfile(full_path):
        msg = "The file specified (%s) does not exist.\nPlease specify a " \
              "valid file." % full_path
        raise IOError(msg)
    
    cmd = ["wigToBigWig", full_path, chromSizes, outfile]
    if verbose:
        print(cmd)

    success=subprocess.check_call(cmd)
    
    if success==0:
        return(outfile)

class wigData:
    
    """A simple class for representing the data in a wig file"""
    
    def __init__(self, fname=None, isBigWig=False, bigWigToWig_binary=None,
                 wigToBigWig_binary=None, process_mode=1):
        
        """ class constructor: this takes the filename of the wig file and
        reads the data from it into the class attributes. filename should
        be a string and point to a file that exists """
        
        # set some defaults
        self.__filename=None
        self.__file_full_path=None
        self.tracks={}
        self.region=None
        self._current_data=None
        self.converted_from_bigwig=False
        self.converted_from_bedgraph=False
        self._generated_by_script="common.parsing_routines.wigtools(%s)" \
                              "" % str(ver)
        
        # try and cast string to a string if it isn't already!
        if type(fname) is str:
                
            # resolve relative paths to a full path
            full_path = os.path.abspath(fname)
            
            # check that the file exists
            if not os.path.isfile(full_path):
                msg = "The file specified (%s) does not exist.\nPlease specify a " \
                      "valid file." % full_path
                raise IOError(msg)
            
            #store the file details
            self.__filename=fname
            self.__file_full_path=full_path
                        
            # read the file data
            if isBigWig:
                self.converted_from_bigwig=True
                fname = convert_BigWig2Wig(fname, self.bigWigToWig_binary)
                filedata = openfile(fname)
                if len(filedata[1].split("\t"))!=1:
                    warnings.warn("Looks like the original file that made your " \
                                  "BigWig was not a wig file! bigWig2Wig has " \
                                  "returned a file thats no in wig format. It's " \
                                  "probably in bedgraph format, so I'll try " \
                                  "converting the file from bedgraph to wig. If " \
                                  "its not in bedgraph format this will fail!")
                    try:
                        self.converted_from_bedgraph=True
                        fname = convert_bedgraph2Wig(fname)
                        filedata = openfile(fname)
                        filedata = numpy.array(filedata)
                    except:
                        print("OK, the attempt to convert from bedgraph format " \
                              "failed. I give up. Get a better file!")
                        raise
                else:
                    filedata = self.trimBigWigMultiLines(filedata)
    
            else:
                filedata = openfile(fname)
                filedata = numpy.array(filedata)
                    
            # this does a regular expression match across a large numpy array fast!
            r=re.compile("^[A-Za-z].+$")
            vmatch = numpy.vectorize(lambda x:bool(r.match(x)))
            sel_line_indexes = numpy.where(vmatch(filedata))[0]
            track_indexes = self.parseWigStructure(filedata, 
                                                   index=sel_line_indexes)        
            
            first_chrid=None
            for index in track_indexes:
                newTrack = wigTrack(filedata[index[0]:index[1]])
                self.tracks[newTrack.chr] = newTrack
                if first_chrid is None:
                    first_chrid = newTrack.chr
            
            # set default region to first chromosome
            self.set_region(first_chrid)
            self._current_data=self.tracks[self.region.chrid]
        else:    
            warnings.warn("The filename passed to wigData is " \
                          "not string. Creating empty container.")
    
    def trimBigWigMultiLines(self, datalines, return_keep=False):
        
        """ the bigWigToWig conversion introduces multiple track lines for a 
        chromosome with each track being ~1024 entries long. This is pissing 
        annoying! and needs to be removed in order to get one track per 
        chromosome. There might at some point be a valid reason for having more
        than one track for a chromosome however (although TBH I can't think of 
        one), so we only remove multiple entries from bigWig conversions, 
        rather than from everything."""
        
        # this does a regular expression match across a large numpy array fast!
        datalines = numpy.array(datalines)            
        r=re.compile("^[A-Za-z].+$")
        vmatch = numpy.vectorize(lambda x:bool(r.match(x)))
        mask = vmatch(datalines)
        mask_index = numpy.where(mask)[0]
        
        
        tracklines=datalines[mask_index]
        print(tracklines[0:10])
        iter_len=len(tracklines)
        i=0
        keep_inds=[]
        exclude=None
        while i<iter_len:
            # this does a regular expression match across a large numpy array fast!
            matchstr=re.match("^(.+ chrom=.+?)( .*|\n)$",tracklines[i]).group(1)
            r=re.compile("^%s .*$" % matchstr)
            vmatch=numpy.vectorize(lambda x:bool(r.match(x)))
            sel_line_indexes=numpy.where(vmatch(tracklines))[0]
            keep_inds.append(sel_line_indexes[0])
            if exclude is None:
                exclude=sel_line_indexes[1:]
            else:
                exclude=numpy.hstack((exclude,sel_line_indexes[1:]))
            
            i=sel_line_indexes[-1]+1
        
        newmask = numpy.ones(len(datalines), dtype='bool')
        newmask[mask_index[exclude]]=False        
        
        return(datalines[newmask])
    
    def parseWigStructure(self, datalines, index=None):
        
        """ parse the structure of the wig data ready to make tracks from """
        
        maxind = len(datalines)
        if index is not None:
            datalines = datalines[index]
        
        track_indexes=[]
        new_start_found=False
        track_start=0
        track_end=0
                
        for line in datalines:
            line = line.strip()
            if line!="":
                if not new_start_found:
                    if (line.startswith("fixedStep") or 
                        line.startswith("variableStep")):
                        new_start_found = True
                else:
                    if not re.match("^[0-9]",line):
                        if (line.startswith("fixedStep") or 
                            line.startswith("variableStep")):
                            new_start_found = True
                        else:
                            new_start_found=False
                        
                        if index is None:    
                            track_indexes.append([track_start, track_end])
                        else:
                            track_indexes.append([index[track_start], index[track_end]])
                        
                        track_start = track_end
                        
            track_end+=1
        
        # add in final end-point
        if index is None:
            track_indexes.append([track_start, track_end])
        else:
            track_indexes.append([index[track_start], maxind])

        return(track_indexes)
    
    def set_region(self, region_rep, start=None, stop=None):
        
        ''' sets the current region within the wigData object '''
        
        this_region = generic_set_region(region_rep, start, stop,
                                         generated_by=self._generated_by_script)
        
        # check selected chromosome is in the wigData tracks
        if this_region.chrid not in self.tracks.keys():
            msg = "Invalid region specification, chromosome id not found. " \
                  "Regions should be specified either as a string of format " \
                  "'chrid:start-stop' or with three arguments as " \
                  "set_region('chrid', [start_int], [stop_int])."
            raise ValueError(msg)

        # sanitize stupid region start/stop values
        if this_region.start is None:
            this_region.start=1
        if this_region.stop is None:
            this_region.stop=int(self.tracks[this_region.chrid].position.max()+1)
        
        self.region = this_region
        self.tracks[self.region.chrid].set_region(self.region)
        self._current_data = self.tracks[self.region.chrid]
    
    def clear_region_selection(self):
        
        """ clears the variables defining the current feature selection """
        
        self.region=None
        self._current_data=None
       
    def display(self, fmt="gff3", limit=None, printthem=True):
        
        printstrs=["Wiggle track data with %i chromosome tracks" \
                   "" % len(self.tracks.keys()),
                   "Original File: %s" % self.__file_full_path]
        if fmt=="gff3":
            printstrs.append("Currently set region: %s" \
                             "" % self.region.get_gff3line())
        elif fmt=="gtf":
            printstrs.append("Currently set region: %s" \
                             "" % self.region.get_gtfline())
                
        if printthem:
            print("\n".join(printstrs))
        else:
            return("\n".join(printstrs))

    def get_region_length(self):
        
        ''' return length, in base-pars of the selected region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self.region.stop-self.region.start)
        
    def get_region_ndatapoint(self):
        
        ''' return number of datapoints in the selected region '''
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_ndatapoint())
    
    def get_region_data(self, expand=False):
        
        ''' return the position and data information within the set region '''
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_data())

    def set_region_data(self, newdata):
        
        ''' sets the data information within the set region '''
        
        self.tracks[self.region.chrid].set_region_data(newdata)
            
    def get_region_min(self):
        
        ''' return the minimum value and the positions with this value 
        within the set region '''
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_min())

    def get_region_max(self):
        
        ''' return the maximum value and the positions with this value 
        within the set region '''
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_max())

    def get_region_mean(self):
        
        ''' return the mean value for the set region '''
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_mean())

    def get_region_median(self):
        
        ''' return the median value for the set region '''
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_median())

    def get_region_stddev(self):
        
        ''' return the std deviation for the set region '''
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_stddev())
    
    def get_region_stderr(self):
        
        ''' return the std error on the mean for the set region '''
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_stderr())
    
    def get_region_stats(self):
        
        """ return mean, median, stddev ans stderr for the region """
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_stats())
    
    def get_region_fracbases(self):
        
        """ return the fraction of bases in the region that have a signal """
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_fracbases())
    
    def get_region_mean_per_base(self):
        
        ''' return the maximum value and the positions with this value 
        within the set region '''
        
        if self._current_data is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self._current_data.region_mean_per_base())

    def writeWig(self, filename, isBigWig=False, fulldata=False, name=None,
                 desc=None, rgbcolstr=None, chromSizes=None,
                 multiple_definitionlines = True, verbose=False):
        
        ''' writes all, or a subset, of the wig data to a new output wigfile 
        
        chromSizes is a dictionary of key, value pairs comprising chromosome 
        names and lengths.'''
        
        outfile = filename
        if isBigWig:
            if verbose:
                print("isBigWig=True; defining temp dir and temp output wig " \
                      "file.")
            tmp_path = tempfile.mkdtemp()
            outfile = "%s/tmpwig.wig" % tmp_path
        
        tracks = self.tracks.keys()
        if not fulldata:    
            tracks = [self.region.chrid]

        if desc is not None:
            desc = "%s. (%s)" % (desc, self._generated_by_script)
        else:
            desc = self._generated_by_script
        
        if verbose:
            print("writing track data...")
        
        skipdefinitionline=False
        firstwrite=True
        for track in tracks:
            if verbose:
                print("writing data for chromosome %s..." % track)
            thistrack = self.tracks[track]
            if firstwrite:
                writemode="w"
            else:
                writemode="a"
            
            if isBigWig or not multiple_definitionlines:
                skipdefinitionline = True
            
            thistrack.writeWig(outfile, fulldata=fulldata, name=name,
                               desc=desc, rgbcolstr=None, writemode=writemode,
                               skipdefinitionline=skipdefinitionline)
            
            firstwrite=False
            
        if isBigWig:
            
            print("isBigWig=True; writing chromSizes file and converting wig " \
                  "to BigWig...")
            
            if chromSizes is None:
                chromSizes={}
                for track in self.tracks.keys():
                    chromSizes[track]=self.tracks[track].position.max()
            
            sizefile = "%s/chromSizes.txt" % tmp_path
            fh = open(sizefile, "w")
            for key in chromSizes:
                fh.write("%s\t%i\n" % (key, chromSizes[key]))
            fh.close()
            
            convert_Wig2BigWig(outfile, sizefile, filename, verbose=False)

class wigTrack:
    
    """A simple class for representing a chromosomes data from a wig file"""
    
    def __init__(self, trackdata=None):
        
        """ class constructor: this takes the filename of the wig file and
        reads the data from it into the class attributes. filename should
        be a string and point to a file that exists """
        
        self.wigtype=None
        self.chr=None
        self.start=None
        self.step=None
        self.span=1
        self.position=None
        self.data=None
        self.region=None
        self._region_index=None
        
        # process the header, stripping out and saving headerlines, reading the 
        # formatting information.
        if trackdata is not None:
            headerlines=[]
            start_found=False
            data_start=0
            for line in trackdata:
                line = line.strip()
                data_start+=1
                if not start_found:
                    if (line.startswith("fixedStep") or 
                        line.startswith("variableStep")):
                        start_found = True
                        headline = line.split()
                        if headline[0]=="fixedStep" or headline[0]=="variableStep":
                            self.wigtype = headline[0]
                        for info in headline[1:]:
                            match = re.match("chrom=(.+)", info)
                            if match:
                                self.chr = match.group(1)
                            match = re.match("start=([0-9]+)", info)
                            if match:
                                self.start = int(match.group(1))                    
                            match = re.match("step=([0-9]+)", info)
                            if match:
                                self.step = int(match.group(1))
                            match = re.match("span=([0-9]+)", info)
                            if match:
                                self.span = int(match.group(1))
                        break;
                    else:
                        headerlines.append(line)
    
            self.header = headerlines
            
            # check wig file formating before saving the data
            if self.wigtype is None:
                raise ValueError("No type specification found in wigfile. Please " \
                                 "specify fixedStep or variableStep.")
    
            if self.chr is None:
                raise ValueError("No chromosome specification found in wigfile. " \
                                 "Use 'chrom=xxxx' to specify a chromosome name.")
            
            if self.wigtype is "fixedStep" and self.start is None:
                raise ValueError("No start specified for fixedStep wigfile. " \
                                 "Use 'start=x' to specify a start position.")
    
            if self.wigtype is "fixedStep" and self.step is None:
                raise ValueError("No step-size specified for fixedStep wigfile. " \
                                 "Use 'step=x' to specify a start position.")
            
            if self.wigtype is "fixedStep" and self.step==1 and self.span>1:
                raise ValueError("No span-size greater than step-size for " \
                                 "fixedStep wigfile.")
    
            # save the data, but leave the data and position info blank for chromosomes
            # without data.
            if len(trackdata)>data_start: 
                if self.wigtype=="fixedStep":
                    in_data = numpy.array(trackdata[data_start:], dtype="float")
                    if self.step==1:
                        self.position=numpy.arange(len(in_data))+self.start
                        self.data=in_data
                    elif self.span==1:
                        self.position=(numpy.arange(len(in_data))*self.step)+self.start
                        self.data = in_data
                    else:
                        dstack = in_data
                        i=1
                        while i<self.span:
                            dstack = numpy.vstack((dstack,in_data))
                            i+=1
                        self.data = dstack.flatten('F') # flatten the 2d array fortran style
                    
                        pstack = numpy.arange(self.span)
                        i=1
                        while i<len(in_data):
                            pstack = numpy.vstack((pstack, numpy.arange(self.span)))
                            i+=1
                        steps = (numpy.arange(len(in_data))*self.step)+self.start
                        pos2d = pstack+numpy.transpose(numpy.atleast_2d(steps))
                        self.position = pos2d.flatten() # flatten the 2d array fortran style
                else:
                    in_data = numpy.fromstring(" ".join(trackdata[data_start:]), sep=" ")
                    in_data = in_data.reshape(len(in_data)/2,2)
                    ipos = in_data[:,0]
                    idata = in_data[:,1]
                    if self.span>1:
                        i=1
                        while i<self.span:
                            ipos = numpy.vstack((ipos, in_data[:,0]+i))
                            idata = numpy.vstack((idata, in_data[:,1]))
                            i+=1
                
                    self.position = ipos.flatten('F').astype(int)
                    self.data = idata.flatten('F')
            
                self.set_region(self.chr, 1, int(self.position.max()+1))
        else:
            warnings.warn("The data passed to wigTrack is " \
                          "empty. Creating empty container.")
    
    def set_region(self, region_rep, start=None, stop=None):
                
        ''' sets the current region within the wigTrack object '''

        if start is None and stop is not None:
            start=1
        elif start is not None and stop is None:
            stop=int(self.position.max()+1)
                
        self.region = generic_set_region(region_rep, start, stop)
    
        if self.region.start is None and self.region.stop is None:
            self.region.start=1
            self.region.stop=int(self.position.max()+1)

        ind1 = numpy.where(self.position>=self.region.start)
        ind2 = numpy.where(self.position<self.region.stop)
        self._region_index = numpy.intersect1d(ind1[0],ind2[0])
    
    def region_length(self):
        
        ''' return length, in base-pars of the selected region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self.region.stop-self.region.start)
        
    def region_ndatapoint(self):
        
        ''' return number of datapoints in the selected region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(len(self._region_index))
    
    def region_data(self, expand=False):
        
        ''' return the position and data information within the set region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
        
        if not expand:
            return((self.position[self._region_index],
                    self.data[self._region_index]))
        else:
            e_pos = numpy.arange(self.region.stop)+1
            e_data = numpy.zeros(self.region.stop, dtype='float')
            e_data[self.position[self._region_index]] = self.data
            ind1 = numpy.where(e_pos>=self.region.start)
            ind2 = numpy.where(e_pos<self.region.stop)
            e_index = numpy.intersect1d(ind1[0],ind2[0])
            return((e_pos[e_index],e_data[e_index]))
    
    def set_region_data(self, newdata):
        
        ''' sets the data information within the set region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
        
        if type(newdata) is not numpy.ndarray:
            try:
                newdata = numpy.array(newdata)
            except:
                msg = "Cannot cast data to numpy array. New data not set."
                raise TypeError(msg)
        
        if len(newdata)!=len(self._region_index):
            msg = "New data array is not the same length (%i) as the current " \
                  "region selected (%i bp). Not setting data."
            raise ValueError(msg)
        
        self.data[self._region_index]=newdata    
    
    def region_min(self):
        
        ''' return the minimum value and the positions with this value 
        within the set region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
        
        minval = self.data[self._region_index].min()
        val_index = numpy.where(self.data[self._region_index]==minval)[0]
        val_positions = self.position[self._region_index][val_index]
                                        
        return((minval,val_positions))

    def region_max(self):
        
        ''' return the maximum value and the positions with this value 
        within the set region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
        
        maxval = self.data[self._region_index].max()
        val_index = numpy.where(self.data[self._region_index]==maxval)[0]
        val_positions = self.position[self._region_index][val_index]
                                        
        return((maxval,val_positions))

    def region_mean(self):
        
        ''' return the mean value for the set region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
                                                
        return(self.data[self._region_index].mean())

    def region_median(self):
        
        ''' return the median value for the set region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
                                                
        return(numpy.median(self.data[self._region_index]))

    def region_stddev(self):
        
        ''' return the std deviation for the set region 
        
        Here we correct numpys typical standard deviation routine to make it an
        unbiased estimator, according to GURLAND, J. & TRIPATHI, R. C. 1971. 
        "Simple Approximation for Unbiased Estimation of Standard Deviation." 
        American Statistician, 25, 30. Note that this correction will propogate 
        through to the standard error on the mean. 
        '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
        
        return(self.data[self._region_index].std(ddof=1)*
               (1.0+(1.0/(4.0*(self.region_ndatapoint()-1.0)))))
    
    def region_stderr(self):
        
        ''' return the std error on the mean for the set region '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
                                                
        return(self.region_stddev()/self.region_ndatapoint())
    
    def region_stats(self):
        
        """ return mean, median, stddev ans stderr for the region """
        
        return({"mean":self.region_mean(),
                "median":self.region_median(),
                "standard_deviation":self.region_stddev(),
                "standard error": self.region_stderr()})
    
    def region_fracbases(self):
        
        """ return the fraction of bases in the region that have a signal """
        
        return(float(len(self._region_index))/self.region_length())
    
    def region_mean_per_base(self):
        
        ''' return the sum of the data divided by its length '''
        
        if self.region is None:
            msg = "No region set"
            raise ValueError(msg)
                                                
        return(self.data[self._region_index].sum()/self.region_length())
        
    def writeWig(self, filename, fulldata=False, name=None, desc=None,
                 rgbcolstr=None, writemode="a", skipdefinitionline=False):
        
        ''' writes all, or a subset, of the tract data to the filehandle 
        
        Optionally can include a name, desc and colour for the track, specify
        whether to write all the data or just the data for the region selected,
        and whether to append to a file, or overwrite.'''
        if self.data is not None and self.position is not None:
            fh=None
            try:
                fh = open(filename, writemode)
            except:
                msg = "cannot open output wig file for writing: %s" % filename
                raise IOError(msg)
        
            if not skipdefinitionline:
                definitionline_opts = ["track", "type=wiggle_0"]
                if name is not None:
                    definitionline_opts.append('name="%s"' % name)            
                if desc is not None:
                    definitionline_opts.append('description="%s"' % desc)
                if rgbcolstr is not None:
                    definitionline_opts.append("color=%s" % rgbcolstr)
                fh.write("%s\n" % " ".join(definitionline_opts))
        
            regstart = self.position[0]
            wdata = self.data
            if not fulldata:
                regstart = self.region.start
                wdata = self.region_data()[1]
        
            declarationline_opts = ["%s" % self.wigtype, "chrom=%s" % self.chr]
            if self.wigtype=="fixedStep":
                declarationline_opts.append("start=%i" % regstart)
                declarationline_opts.append("step=%i" % self.step)
        
            fh.write("%s\n" % " ".join(declarationline_opts))
        
            match = re.match("^0.(0+)[1-9][0-9]*$", str(wdata.min()))
            ndp=1
            if match:
                ndp = len(match.group(1))+1
        
            if self.wigtype=="fixedStep":
                for value in wdata:
                    fstring = "%s.%if\n" % ("%", ndp)
                    fh.write(fstring % value)
            else:
                i=0
                while i<len(wdata):
                    fstring = "%s.%if\n" % ("%i %", ndp)
                    fh.write(fstring % (self.position[i], wdata[i]))
                    i+=1
        
            fh.close()

def mathWigs(wig1, wig2, operation="sum", verbose=False):
    
    """performs a maths operation on two wig files, like summing them """
    
    # find overlapping chrs and unique chrs in each set of wig data. maths will only be
    # done on data that exists on matching chrs for both datasets
    intersect_Chr = set.intersection(set(wig1.tracks.keys()), set(wig2.tracks.keys()))
    wig1only_Chr = set.difference(set(wig1.tracks.keys()), set(wig2.tracks.keys()))
    wig2only_Chr = set.difference(set(wig2.tracks.keys()), set(wig1.tracks.keys()))
    
    newwig = copy.deepcopy(wig1)
    
    # add unique data tracks to the output
    for val in wig1only_Chr:
        newwig.tracks[val]=wig1.tracks[val]
    
    for val in wig2only_Chr:
        newwig.tracks[val]=wig2.tracks[val]
    
    for val in intersect_Chr:
        if wig1.tracks[val].position is None:
            newwig.tracks[val]=wig2.tracks[val]
        elif wig2.tracks[val].position is None:
            newwig.tracks[val]=wig1.tracks[val]
        else:
            union_positions = numpy.union1d(wig1.tracks[val].position, wig2.tracks[val].position)
            new_data = numpy.zeros(len(union_positions), dtype=float)
            wig1_index = numpy.in1d(union_positions, wig1.tracks[val].position)
            wig2_index = numpy.in1d(union_positions, wig2.tracks[val].position)
        
            new_data[wig1_index] = wig1.tracks[val].data
            if operation=="sum":
                new_data[wig2_index] = new_data[wig2_index] + wig2.tracks[val].data
            elif operation=="diff":
                new_data[wig2_index] = new_data[wig2_index] - wig2.tracks[val].data
            elif operation=="mean":
                new_data[wig2_index] = numpy.mean(numpy.vstack((new_data[wig2_index], wig2.tracks[val].data)),axis=0)        
        
            newwig.tracks[val].position = union_positions
            newwig.tracks[val].data = new_data
    
    newwig._wigData__file_full_path=""
    newwig._wigData__filename = "%s of %s & %s" % (operation, wig1._wigData__filename, wig2._wigData__filename)
    newwig.converted_from_bedgraph=False
    newwig.converted_from_bigwig=False
        
    return(newwig)
    
def scaleWig(infile, scaleFactor, outfile, isBigWig=False,
             bigWigToWig=None, wigToBigWig=None, verbose=False):
    
    """ reads, scales and output bigwig data"""
    
    thiswig = wigData(infile, isBigWig=isBigWig, wigToBigWig_binary=wigToBigWig,
                      bigWigToWig_binary=bigWigToWig)
    
    for key in thiswig.tracks.keys():
        thiswig.set_region(key)
        data = thiswig.get_region_data()
        scaleddata = data[1]*scaleFactor
        thiswig.set_region_data(scaleddata)
    
    return(thiswig)
    
    
    
    