"""
=============================
Read pooling pipeline
=============================

:Author: Matt Jackson
:Release: $Id$
:Date: |today|
:Tags: Python


Overview
========

Pool reads from multiple fasta/q files into one with optional dereplication of the resulting file.
  
Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 


Input
-----

Reads
+++++

Reads are imported by placing files are linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   <filename>.<suffix>

The ``suffix`` determines the file type.  The
following suffixes/file types are possible:

fasta,fastq

fasta.1,fastq.1

fasta.gz,fastq.qz

fasta.1.gz,fastq.1.gz

NOTE:Different formats are considered, but all filles must be the same within one pooling run

Where unnumbered files relate to single end reads or inter-leaved paired-end files.
Interleaved files should be detected automatically.
Files containg 1 in the suffix will be assumed paired end and should be in the same 
directory as the read 2 files, which share the same file name except for the read number.

.. note::

   Quality scores need to be of the same scale for all input files. Thus it might be
   difficult to mix different formats.

Optional inputs
+++++++++++++++

Requirements
------------

On top of the default CGAT setup and the MetaAssemblyKit pipeline, the pipeline requires the following
software to be in the path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|  BBMap             | 38.0.0            | Dereplicate pooled reads                       |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The main output is a fasta/q file(s) containing all the input reads and/or a dereplicated copy. 

Glossary
========

.. glossary::

Code
====

"""

#load modules
from ruffus import *
import os
import sys
import math

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
       ["%s.ini" % __file__[:-len(".py")],
       "../pipeline.ini",
       "pipeline.ini" ] )
PARAMS = P.PARAMS

#add Pipeline functions
import PipelineMetaAssemblyKit
import PipelinePoolReads

#get all files within the directory to process
SEQUENCEFILES = ("*.fasta", "*.fasta.gz", "*.fasta.1.gz", "*.fasta.1",
                 "*.fna", "*.fna.gz", "*.fna.1.gz", "*.fna.1",
                 "*.fa", "*.fa.gz", "*.fa.1.gz", "*.fa.1", 
                 "*.fastq", "*.fastq.gz", "*.fastq.1.gz","*.fastq.1")

#using log files to track progression as ruffus cant parse inputs within merge

####################################################
# Merge the inputs into one pooled file
####################################################
@follows(mkdir("pooled.dir"))
@merge(SEQUENCEFILES,"pooled.dir/pool.log")
def poolReads(infiles,outfile):
    statementlist = []
    #get file type from first file
    ftype = PipelineMetaAssemblyKit.SequencingData(infiles[0])
    #generate output filename
    outname = "pooled.dir/"+PARAMS["General_output_prefix"]+"."+ftype.fileformat
    #pool the reads
    statementlist.append(PipelinePoolReads.poolReads(ftype,infiles,outname))
    #create the log to ensure jpb isn't rerun
    statementlist.append('echo "Pooled {} files to {}" >> pooled.dir/pool.log'.format(len(infiles),outname))
    statement = " && ".join(statementlist)
    P.run()

###########################################################
# Add numbers to to read names to make unique if required
###########################################################
@active_if(PARAMS["General_unique_names"]=="true")
@follows(poolReads)
@merge(SEQUENCEFILES,"pooled.dir/rename.log")
def renameNumbers(infiles,outfile):
    statementlist = []
    #get the pooled output file location from the file name and params
    pooledfile = PipelinePoolReads.pooledName(infiles,PARAMS)
    #process the sequencing data object from the pooled output
    statementlist.append(PipelinePoolReads.renameNums(pooledfile,os.getcwd()+"/pooled.dir/"))
    statementlist.append('echo "Renamed reads within {}" >> pooled.dir/rename.log'.format(pooledfile.filename))
    statement = " && ".join(statementlist)
    P.run()
       

###########################################################
# Remove duplicate reads if required
###########################################################
@active_if(PARAMS["General_dereplicate"]=="true")
@follows(renameNumbers)
@merge(SEQUENCEFILES,"pooled.dir/dereplicate.log")
def dereplicatePooled(infiles,outfile):
    statementlist = []
    #get the pooled (and potentially renamed) output file location from the filename and params
    pooledfile = PipelinePoolReads.pooledName(infiles,PARAMS)
    #generate the dereplication statement for CD-HIT
    statementlist.append(PipelinePoolReads.derepReads(pooledfile,os.getcwd()+"/pooled.dir/",PARAMS))
    #remove the original before derpelication (based on params) 
    if PARAMS['General_keep_with_reps'] == "false":
        statementlist.append("rm {}".format(os.getcwd()+"/pooled.dir/"+pooledfile.filename))
        if pooledfile.paired == True:
            statementlist.append("rm {}".format(os.getcwd()+"/pooled.dir/"+pooledfile.pairedname))
    #make log
    statementlist.append('echo "Dereplicated reads within {}" >> pooled.dir/dereplicate.log'.format(pooledfile.filename))
    statement = " && ".join(statementlist)
    job_memory = "{}G".format(math.ceil(os.path.getsize(os.getcwd()+"/pooled.dir/"+pooledfile.filename)/1e9))
    job_threads = 2
    P.run()
    
@follows(dereplicatePooled)
def full():
    pass
    

if __name__ == "__main__":
    if sys.argv[1] == "plot":
        pipeline_printout_graph("test.pdf", "pdf", [full], no_key_legend=True,
                                size=(4, 4),
                                user_colour_scheme = {"colour_scheme_index": 1})
    else:
        sys.exit(P.main(sys.argv))



        
