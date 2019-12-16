'''
-------------------
custom_callables.py
-------------------
A module containing common-use subroutines for options parsing in scripts. 
These are designed to be used with the `argparse package 
<http://docs.python.org/2.7/library/argparse.html>`_ that should be 
used for all argument parsing from this point on!

Basically these are a bunch of custom callables that takes a single string 
argument and return a converted value for use with the 'type' option.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.0
:created_on: 2013-04-05

'''

__version__ = "1.0"

import argparse, os, warnings

def input_path(string):
    
    ''' Checks to see if the string is a valid input path
    
    Returns the full path of the string if true, otherwise raise a custom 
    exception.
    '''
    
    # try and cast string to a string if it isn't already!
    if type(string) is not str:
        raise TypeError("The value passed to input_path is not a string. " \
                        "value: %s, type: %s" % (str(string), type(string))
                        )
    
    # resolve relative paths to a full path
    path = os.path.abspath(string)
    
    if os.path.exists(path):
        return(path)
    else:
        msg = "The path specified (%s) does not exist.\nPlease specify a " \
              "valid path. See -h|--help for help with running this script." \
              "" % path
        raise argparse.ArgumentTypeError(msg)

def input_file(string):
    
    ''' Checks to see if the string is a valid input file
    
    Returns the full path of the string if true, otherwise raise a custom 
    exception.
    '''
    
    # try and cast string to a string if it isn't already!
    if type(string) is not str:
        raise TypeError("The value passed to input_file is not a string. " \
                        "value: %s, type: %s" % (str(string), type(string))
                        )
    
    # resolve relative paths to a full path
    path = os.path.abspath(string)
    
    if os.path.isfile(path):
        return(path)
    else:
        msg = "The file specified (%s) does not exist.\nPlease specify a " \
              "valid file. See -h|--help for help with running this script." \
              "" % path
        raise argparse.ArgumentTypeError(msg)

def output_path(string):
    
    ''' Checks to see if the string is a valid output path; makes it if not'''
    
    # try and cast string to a string if it isn't already!
    if type(string) is not str:
        raise TypeError("The value passed to output_path is not a string. " \
                        "value: %s, type: %s" % (str(string), type(string))
                        )
    
    # resolve relative paths to a full path
    path = os.path.abspath(string)
    
    if not os.path.exists(path):
        msg = "The path specified (%s), does not exist... creating it." % path
        warnings.warn(msg, UserWarning)
        os.makedirs(path)
    
    return(path)

def output_file(string):
    
    ''' Checks to see if the the dir of the string is a valid output path
    
    If not, it makes it and then returns the absolute path for the file'''
    
    # try and cast string to a string if it isn't already!
    if type(string) is not str:
        raise TypeError("The value passed to output_file is not a string. " \
                        "value: %s, type: %s" % (str(string), type(string))
                        )
    
    # resolve relative paths to a full path
    path = os.path.dirname(os.path.abspath(string))
    
    if not os.path.exists(path):
        msg = "The path (%s) for the specified file (%s), does not exist... " \
              "creating it." % (path, os.path.basename(string))
        warnings.warn(msg, UserWarning)
        os.makedirs(path)
    
    filepath = os.path.join(path, os.path.basename(string))
    return(filepath)