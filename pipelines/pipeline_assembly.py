"""
=============================
Metagenome assembly pipeline
=============================

:Author: Nick Ilott & Matt Jackson
:Tags: Python

The assembly pipeline takes reads from meta'omic sequencing 
in fasta/fastq files and assembles them into contigs.

Overview
========

The pipeline assembles all the reads within a given fasta/fastq file using a range of possible 
assembly tools designed for metagenomic assembly. It can handle both single and 
paired end reads and provides summary statistics relating to the assembly.

Assembly stategy
----------------

While there exist many tools for assembling reads from single genomes,
only recently has software been specifically developed (or extended)
to allow for the assembly of metagenomes. The major factor affecting
the ability to assemble a metagenome is the diversity of the sample.
Accurate assembly of long contigs is difficult in the presence of many
species at differential abundances. This is in constrast to single
genome assembly where reads are (or should be) uniformly sampled
across the genome.

This pieline provides a range of available software for the assembly of metagenomes.

Considerations
--------------

Metagenomics is a young and rapidly developing field. There is, as
yet, no gold standard for assembly. It is likely that the presence of
multiple, highly related species will lead to the assembly of
chimeric contigs i.e. contigs derived from more than one species. It
is generally considered that longer K-mer lengths used in the
construction of the de-bruijn graph (for de-bruijn graph assemblers)
will result in fewer chimeras. Nevertheless, longer k-mers may also
result in more, short contigs being produced as a result of a
neccessity for a greater overlap between nodes in the graph. The
length of k-mer chosen is also dependent on the length of reads that
you are trying to assemble - longer reads means you can use longer
k-mers. Which k-mer to use in the assembly process is therefore
dependent on the data used and the expected complexity of the
sample. We make no effort here to advise on k-mer length.
Metatranscriptomic data does not fit the assumptions made by
many of the metagenomic assembly algorithms and the transcriptomic
assembly tools assume RNA from a single source organism. So whilst
this pipeline can be used for these data, the results must be considered
with caution and within this context.

  
Usage
=====

See cgat.org for general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file. 


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

Where unnumbered files relate to single end reads or inter-leaved paired-end files.
Interleaved files should be detected automatically.
Files containg 1 in the suffix will be assumed paired end and should be in the same 
directory as the read 2 files, which share the same file name except for the read number.
This pipeline has been largely tested/developed using paired-end, non-interleaved, fastq files.

Optional inputs
+++++++++++++++

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path (individual assemblers can be skipped if not required):

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Tested Version*   |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|megahit             |1.1.3              |Metagenomic Assembler                           |
+--------------------+-------------------+------------------------------------------------+
|spades              |3.11.1             |Metagenomic Assembler                           |
+--------------------+-------------------+------------------------------------------------+
|IDBA-UD             |1.1.1              |Metagenomic Assembler                           |
+--------------------+-------------------+------------------------------------------------+
|seqtk               |1.0.r82-dirty      |Tools for reformating sequencing data           |
+--------------------+-------------------+------------------------------------------------+
|BBMap               |38.0.0             |Short read aligner with tools for summary stats |
+--------------------+-------------------+------------------------------------------------+

NOTE: Spades and IDBA-UD required paired reads, and IDBA-UD reads must be <128bp (with default compile)

Pipeline output
===============

The main output is the genome assembly (contigs) - output as a fasta formatted file.
Additional outputs include sumary statistics relating to contig size etc.

Running build_report will generate a summary of the contigs for each sample and each assembly method used.

Code
====

"""

#load modules
from ruffus import *
import os
import sys
import math

from pipeline_assembly import PipelineAssembly


###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
# load options from the config file
import cgatcore.pipeline as P
P.get_parameters(
       ["%s/pipeline.yml" % __file__[:-len(".py")],
       "../pipeline.yml",
       "pipeline.yml" ] )
PARAMS = P.PARAMS


#get all files within the directory to process
SEQUENCEFILES = ("*.fasta", "*.fasta.gz", "*.fasta.1.gz", "*.fasta.1",
                 "*.fna", "*.fna.gz", "*.fna.1.gz", "*.fna.1",
                 "*.fa", "*.fa.gz", "*.fa.1.gz", "*.fa.1", 
                 "*.fastq", "*.fastq.gz", "*.fastq.1.gz","*.fastq.1")

SEQUENCEFILES_REGEX = regex(
    r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)")


###################################################
# Check the input file, count number of reads
###################################################
@follows(mkdir("filesummaries.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"filesummaries.dir/\1.seqsummary")
def checkFile(infile, outfile):
    seqdat=PipelineAssembly.SequencingData(infile)
    outf=open(outfile,'w')
    outf.write("name\t{}\nformat\t{}\ncompressed\t{}\npaired\t{}\ninterleaved\t{}\n".format(
        seqdat.filename,seqdat.fileformat,seqdat.compressed,seqdat.paired,seqdat.interleaved))
    seqdat.readCount()
    outf.write("read_count\t{}\n".format(seqdat.readcount))
    outf.close()

##################################################
#Run Selected Assemblers
##################################################
    
#get the list of assemblers to run on the data 
ASSEMBLERS = P.as_list(PARAMS.get("Assembler_assemblers", ""))
    
###################################################
# Run Megahit
###################################################
@active_if("megahit" in ASSEMBLERS)
@follows(checkFile)
@follows(mkdir("megahit_out.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"megahit_out.dir/\1/\1.contigs.fa")
def runMegahit(infile, outfile):
    job_memory = str(PARAMS["Megahit_clus_memory"])+"G"
    job_threads = int(PARAMS["Megahit_clus_threads"])
    seqdat=PipelineAssembly.SequencingData(infile)
    assembler = PipelineAssembly.Megahit(seqdat,"megahit_out.dir",PARAMS)
    statement = assembler.build()
    P.run(statement)

###################################################
# Run Metaspades
###################################################
@active_if("metaspades" in ASSEMBLERS)
@follows(checkFile)
@follows(mkdir("metaspades_out.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"metaspades_out.dir/\1/contigs.fasta")
def runMetaspades(infile,outfile):
    job_memory = str(PARAMS["Metaspades_memory"])+"G"
    job_threads = int(PARAMS["Metaspades_threads"])
    seqdat = PipelineAssembly.SequencingData(infile)
    if seqdat.paired == True:
        assembler = PipelineAssembly.Metaspades(seqdat,"metaspades_out.dir",PARAMS)
        statement = assembler.build()
        P.run(statement)
    else:
        print("cannot run metaspades on file {} as it requires paired-end data".format(seqdat.filename))


###################################################
# Run IDBA-UD
###################################################
#requires additional step to interleave non-interleaved fasta/fastq files
@active_if("idba_ud" in ASSEMBLERS)
@follows(checkFile)
@follows(mkdir("idbaud_out.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"idbaud_out.dir/\1/\1.interleaved.fa")
def idbaudInterleave(infile,outfile):
    seqdat = PipelineAssembly.SequencingData(infile)
    if os.path.exists(os.getcwd()+"/idbaud_out.dir/{}".format(seqdat.cleanname)) == False:
        os.mkdir(os.getcwd()+"/idbaud_out.dir/{}".format(seqdat.cleanname))
    statement = PipelineAssembly.IdbaudInterleave(seqdat,os.getcwd(),outfile)
    if statement != None:
        P.run(statement)

@active_if("idba_ud" in ASSEMBLERS)
@follows(idbaudInterleave)
@transform(idbaudInterleave,regex(r"(\S+)/(\S+).(interleaved.fa)"),r"\1/contig.fa")
def runIdbaud(infile,outfile):
    job_memory = str(PARAMS["IDBAUD_clus_memory"])+"G"
    job_threads = int(PARAMS["IDBAUD_clus_threads"])
    seqdat = PipelineAssembly.SequencingData(infile)
    assembler = PipelineAssembly.Idbaud(seqdat,"idbaud_out.dir",PARAMS)
    statement = assembler.build()
    P.run(statement)

###################################################
# Symlink data together from seperate per file folders
###################################################
@follows(runMegahit,runMetaspades,runIdbaud)
@follows(mkdir("contigs.dir"))
@transform([runMegahit,runMetaspades,runIdbaud],regex(r'(\S+)_out.dir/(\S+)/(\S+)\..*$'),r'contigs.dir/\1/\2.contigs.fasta')
def collateContigfiles(infile,outfile):
    dirs=os.getcwd()+"/"
    outdir="/".join(outfile.split("/")[0:-1])
    print(outdir)
    statement = "mkdir -p {} && ln -s {} {}".format(dirs+outdir,dirs+infile,dirs+outfile)
    P.run(statement)

###################################################
# Generate summary stats for contigs
###################################################
@follows(collateContigfiles)
@transform([runMegahit,runMetaspades,runIdbaud],regex(r'(\S+)/(\S+)/(\S+)\..*$'),r'\1/\2/\3.summary')
def summariseContigs(infile,outfile):
    #summarise each contigs file 
    statement = PipelineAssembly.SummariseContigs(infile,outfile)
    P.run(statement)

@follows(summariseContigs)
@merge(summariseContigs,'contigs.dir/Contigs.Summary')
def mergeSummaries(infiles,summaryfile):
    #file to store all the stats combined
    print(mergeSummaries)
    combstats = os.getcwd()+"/"+summaryfile
    statementlist = []
    in0 = os.getcwd()+"/"+infiles[0] 
    statementlist.append("touch {}".format(combstats))
    statementlist.append("head -1 {} >>{}".format(in0,combstats))
    statementlist.append("sed  -i '1s/^/{}\\t{}\\t /' {}".format("file","assembler",combstats))
    #extract filenames and assembler names to add to summary text file
    for infile in infiles:
        indir=os.getcwd()+"/"+infile
        insplit=infile.split("/")
        filen=insplit[1]
        assem=insplit[0].split("_")[0]
        #just append the last line and add filename and assembler name
        statementlist.append("tail -1 {} >> {}".format(indir,combstats))
        statementlist.append("sed -i '$s/^/{}\\t{}\\t /' {}".format(filen,assem,combstats))
    statement = " && ".join(statementlist)
    P.run(statement)
    
###################################################
# Clear-up unused directories
###################################################
@follows(mergeSummaries)
def cleanupDirs():
    statementlist=[]
    if "megahit" not in ASSEMBLERS:
        statementlist.append("rm -r "+os.path.join(os.getcwd(),"megahit_out.dir/"))
    if "metaspades" not in ASSEMBLERS:
        statementlist.append("rm -r "+os.path.join(os.getcwd(),"metaspades_out.dir/"))
    if "idba_ud" not in ASSEMBLERS:
        statementlist.append("rm -r "+os.path.join(os.getcwd(),'idbaud_out.dir/'))
    if statementlist != []:
        statement=" && ".join(statementlist)
        P.run(statement)


@follows(cleanupDirs)
def full():
    pass

@follows(mergeSummaries)
@follows(mkdir("report.dir"))
def build_report():
   scriptloc =  "/".join(os.path.dirname(sys.argv[0]).split("/")[0:-1])+"/scripts/assembly_report.Rmd"
   statement = 'R -e "rmarkdown::render(\'{}\',output_file=\'{}/report.dir/assembly_report.html\')" --args {}/contigs.dir/Contigs.Summary'.format(scriptloc,os.getcwd(),os.getcwd())
   P.run(statement)
   
    
if __name__ == "__main__":
    if sys.argv[1] == "plot":
        pipeline_printout_graph("test.pdf", "pdf", [full], no_key_legend=True,
                                size=(4, 4),
                                user_colour_scheme = {"colour_scheme_index": 1})
    else:
        sys.exit(P.main(sys.argv))
