'''
-------------------------------
script_logging.standard_logging
-------------------------------

This module contains a set of commands that initialize 
:py:class:`logging.Logger` objects with standard sets of pre-defined 
options. The idea here is that I have a standard basic logger with set suite
of information logged as a header to the logfile. In addition then, I can also
have a 'cluster logger' with a set of pre-defined cluster oriented logging 
options that can be added to the standard logger for scripts that need
cluster support, etc etc

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.3
:created_on: 2013-04-09

'''

__version__ = "1.3"

import logging, os, sys, re

def get_logfile_handler(script_args):

    ''' Set up a logfile hander with a standard set of formatting options.
    
    The standard logger setup here is designed to work with a set of command
    line options that has been parsed with the `argparse package 
    <http://docs.python.org/2.7/library/argparse.html>`_ passed to the routine.
    The parser should have *verbose* and *log* attributes in the arguments. 
    Here we setup the file hander and the log message formats based on the 
    verbose and log parameters. The resulting handler should be used in the 
    main script to set up the root logger.
    '''
    
    # check for the correct attributes
    try:
        verbose = script_args.verbose
    except AttributeError:
        raise AttributeError("The args passed to standard_logger do not " \
                             "include the required 'verbose' attribute.")
    
    try:
        logfile =  script_args.log
    except AttributeError:
        raise AttributeError("The args passed to standard_logger do not " \
                             "include the required 'log' attribute.")
    
    # open the logfile and set the log format 
    format_str = "\n %(asctime)s : %(message)s"
    if verbose:
        format_str = "\n %(asctime)s : %(name)s : %(levelname)s " \
                     "\n\t%(message)s"
        
    # setup the logfile details
    logfile_handler = logging.FileHandler(logfile, mode="w")
    logfile_handler.setFormatter(logging.Formatter(format_str,
                                                   datefmt='%H:%M:%S'
                                                   )
                                 )

    return(logfile_handler)

def standard_logger(script_version, cmd_line, script_args, pos_args, kw_args,
                    logbase=None, script_name=None, modules=None):
    
    ''' Set up a logger with a standard set of options and an automatic header.
    
    The standard info logged here required that a root logger has already been 
    instantized and is designed to work with a set of command line options that 
    has been parsed with the `argparse package 
    <http://docs.python.org/2.7/library/argparse.html>`_ passed to the routine.
    The raw command line and master script version is also passed to the logger
    and recorded to that 1) its clear what command line was used and what 
    version of the script was run and 2) its clear whether the command line was 
    parsed correctly.
    
    '''

    # check for the correct attributes
    try:
        verbose = script_args.verbose
    except AttributeError:
        raise AttributeError("The args passed to standard_logger do not " \
                             "include the required 'verbose' attribute.")
    
    try:
        logfile =  script_args.log
    except AttributeError:
        raise AttributeError("The args passed to standard_logger do not " \
                             "include the required 'log' attribute.")
    
    std_logger = logging.getLogger('__main__')
    std_logger.setLevel(logging.INFO)
    if script_args.verbose:
        std_logger.setLevel(logging.DEBUG)
    logfile_handler = get_logfile_handler(script_args)
    std_logger.addHandler(logfile_handler)
        
    # write the title and commandline to the log
    if script_name is None:
        script_name=os.path.basename(cmd_line[0])
    
    log_title = "%s v%s" % (script_name, script_version)
    title_lines = "=" * len(log_title)    
    
    loggerstr = "\n" \
                " %s\n" \
                " %s\n" \
                " %s\n" \
                "\n" \
                " Command-line: %s\n" \
                "\n" \
                " Script Arguments: \n" \
                "" % (title_lines, log_title, title_lines, " ".join(cmd_line))
    
    # get maximum argument label length so that the log formatting looks good!
    arg_dict = script_args.__dict__
    maxlen=0
    for arg in arg_dict.keys():
        if len(arg)>maxlen:
            maxlen = len(arg)
    
    # pad it a little
    maxlen=maxlen+3
    
    # report positional arguments first since they should be required! 
    for arg in pos_args:
        argstr = '\t%s:  %s' % (arg[0].ljust(maxlen), arg_dict[arg[0]])
        if arg_dict[arg[0]] == arg[1]:
            argstr = argstr + " (default)"
        loggerstr = loggerstr + argstr + "\n"
    
    # report keyword arguments (denote them with --)
    for arg in kw_args:
        argstr = '\t%s:  %s' % (("--"+arg[1]).ljust(maxlen), arg_dict[arg[0]])
        if arg_dict[arg[0]] == arg[2]:
            argstr = argstr + " (default)"
        loggerstr = loggerstr + argstr + "\n\n"
    

    if verbose:
        # munge module information
        module_strings=["Module imports:\n"]
                
        if modules is None:
            modules = sys.modules
                
        for module in modules.keys():
            try:
                module_strings.append("%s: %s" % (module,
                                                 str(modules[module].version)))
            except AttributeError:
                try:                    
                    module_strings.append("%s: %s" % (module,
                                                     ".".join(str(modules[module].version_info))))
                except AttributeError:
                    module_strings.append("%s" % module)
        
        loggerstr = loggerstr + "\n".join(module_strings) + "\n\n"
        
        # print relevent environment variables
        env_vars=["Environment variables:\n"]
        for val in os.environ:
            varstr = "%s:" % val
            try:
                evar = re.split(":",os.environ[val])
                evarstr = ""
                for line in evar:
                    if line!="":
                        evarstr = "%s %s" % (evarstr, line.strip())
            except KeyError:
                varstr = "%s not set" % varstr
            env_vars.append("%s\t%s" % (varstr,evarstr))
        
        loggerstr = loggerstr + "\n".join(env_vars) + "\n"
    
    #write to log
    std_logger.info(loggerstr)
    
    return(std_logger)
    
     
