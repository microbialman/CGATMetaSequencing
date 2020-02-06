"""
=============================
Metagenome filtering pipeline
=============================

:Author: Matt Jackson

Filters NGS data from metagenomics experiments to remove one or both of rRNA and host contaminant reads.
Optionally retaining filtered reads.

Overview
========

The pipeline uses the SortMeRNA package to filter rRNAs using a given database.
Host contaminant reads are removed by mapping to a given genome and subsequently filtering
mapped reads.

  
Usage
=====

See cgat.org for  general information how to use CGAT pipelines.

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

Requirements
------------

On top of the default CGAT setup and the MetaSequencing pipeline, the filtering requires the following
software to be in the path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Tested Version*   |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|SortMeRNA           |2.1b               | Filter rRNA sequences                          |
+--------------------+-------------------+------------------------------------------------+
|seqtk               |1.0.r82-dirty      |Tools for reformating sequencing data           |
+--------------------+-------------------+------------------------------------------------+
|bowtie2             |2.3.0              |Mapper to align to host reference genome        |
+--------------------+-------------------+------------------------------------------------+
|samtools            |1.3.1              |Convert and handle bam and sam files            |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |2.25.0             |Convert bed to fastx                            |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The main output is fasta/fastq files filtered to remove the chosen rRNA/host reads.
Optionally, the filtered reads can also be retained.

Running build_report will generate a summary of reads filtered per sample.

Code
====

"""

#load modules
from ruffus import *
import os
import sys

from pipeline_filter import PipelineFilter
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

####################################################
# SortMeRNA steps
####################################################
#first make an index for the rRNA reference files
@active_if(PARAMS["General_rrna_filter"] == "true" and PARAMS["SortMeRNA_rna_index"] == "false")
@follows(mkdir("ref_index.dir"))
@transform(PARAMS["SortMeRNA_rna_refs"].split(","),
          regex(r"^(.+)/([^/]+)$"),
           r"ref_index.dir/\2-db.stats"
)
def makeSortMeRNAIndices(infile,outfile):
    #command to generate index files
    statement = "indexdb_rna --ref {},{}".format(infile,outfile.strip(".stats"))
    P.run(statement)

#run SortMeRNA
@follows(makeSortMeRNAIndices)
@follows(mkdir("rrna_filter_out.dir"))
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,
          r"rrna_filter_out.dir/\1/other_\1.\2")
def runSortMeRNA(infile,outfile):
    seqdat = PipelineAssembly.SequencingData(infile)
    if PARAMS["General_rrna_filter"] == "true":
        sortmerna = PipelineFilter.SortMeRNA(seqdat,outfile,PARAMS)
        if PARAMS["SortMeRNA_memory"] != "false":
            job_memory = str(PARAMS["SortMeRNA_memory"])+"G"
        else:
            job_memory = "1G"
        job_threads = int(PARAMS["SortMeRNA_threads"])
        statement = sortmerna.build()
    else:
        #if skipping rRNA filtering symlink files and make appropriate directory
        statementlist= []
        statementlist.append('mkdir -p rrna_filter_out.dir/{}'.format(seqdat.cleanname))
        statementlist.append('ln -s {} {}'.format(os.getcwd()+"/"+infile,outfile))
        if seqdat.paired == True and seqdat.interleaved == False:
            statementlist.append('ln -s {} rrna_filter_out.dir/{}/other_{}'.format(os.getcwd()+"/"+seqdat.pairedname,seqdat.cleanname,seqdat.pairedname))
        statement = " && ".join(statementlist)
    P.run(statement)

###################################################
# Genome Alignment step
###################################################
#align reads to chosen reference genome using bowtie2
#convert the output from sam to bam
@active_if(PARAMS["General_host_filter"] == "true")
@follows(runSortMeRNA)
@follows(mkdir("genome_filter_out.dir"))
@transform(runSortMeRNA,regex(r'rrna_filter_out.dir/(\S+)/other_(\S[^.]+.\S+)$'),
           r"genome_filter_out.dir/\1/\2.mapped.bam")
def mapBowtie2(infile,outfile):
    job_threads = int(PARAMS["Bowtie_threads"])
    job_memory = str(PARAMS["Bowtie_memory"])+"G"
    seqdat = PipelineAssembly.SequencingData(infile)
    bowtie = PipelineFilter.Bowtie2(seqdat,outfile,PARAMS,PARAMS["Bowtie_genome_db"])
    statementlist = []
    #remove all comments from read names in files (trimming can add comments making non-matching readnames in pairs)
    #can skip this if crashed in previous run using the skip_file_prep parameter
    if PARAMS["Bowtie_skip_file_prep"] != "true":
        statementlist.append(bowtie.cleanNames())
    #directory for output
    statementlist.append("mkdir -p {}".format(os.path.dirname(outfile)))
    #call to bowtie
    statementlist.append(bowtie.build())
    #convert sam to bam
    statementlist.append("samtools view -bS {} > {}".format(outfile.replace(".bam",".sam"),outfile))
    #remove the sam file
    statementlist.append("rm {}".format(outfile.replace(".bam",".sam")))
    statement = " && ".join(statementlist)
    P.run(statement)

#filter the reads from the mapping output and convert to fasta/q file(s)
@active_if(PARAMS["General_host_filter"] == "true")
@follows(mapBowtie2)
@transform(mapBowtie2,regex(r"genome_filter_out.dir/(\S+)/(\S+).mapped.bam"),
           r"genome_filter_out.dir/\1/hostfiltered_\2")
def filterMapping(infile,outfile):
    #use the original sequencing file to pull pairedness, file format and compression
    seqdat = PipelineAssembly.SequencingData(os.path.basename(infile.strip(".mapped.bam")))
    filterer = PipelineFilter.FilterFromBam(infile,outfile,seqdat,PARAMS)
    statementlist = []
    statementlist.append(filterer.build())
    statement = " && ".join(statementlist)
    P.run(statement)

#symlink the appropriate outputfiles to a final clean directory
@follows(filterMapping)
@follows(mkdir("filtered_reads.dir"))
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"filtered_reads.dir/filtered-\1.\2")
def cleanUp(infile,outfile):
    seqdat = PipelineAssembly.SequencingData(infile)
    statement = PipelineFilter.CleanUp(seqdat,outfile,PARAMS)
    P.run(statement)
    
    
@follows(cleanUp)
def full():
    pass
    
#report building
#summarise number of reads in each file at each stage
@follows(mkdir("report.dir/per_file_summaries/"))
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"report.dir/per_file_summaries/\1.filtersummary.txt")
def summariseCounts(infile,outfile):
    statement = PipelineFilter.CountReads(infile,outfile,PARAMS)
    P.run(statement)
    
#combine these into one file
@follows(summariseCounts)
@merge(summariseCounts,"report.dir/combined.filtersummary.txt")
def mergeSummaries(infiles, outfile):
    combinedcounts = open(outfile,"w",1)
    combinedcounts.write("File\tInput\tPost_rRNA_Filtering\tPost_Genome_Filtering\n")
    for i in infiles:
        sumfile = open(i,"rU").readlines()
        combinedcounts.write(sumfile[1])
    combinedcounts.close()
    
#generate an HTML report from the summary files
@follows(mergeSummaries)
def build_report():
    scriptloc = "/".join(os.path.dirname(sys.argv[0]).split("/")[0:-1])+"/scripts/filter_report.Rmd"
    statement = 'R -e "rmarkdown::render(\'{}\',output_file=\'{}/report.dir/filter_report.html\')" --args {}/report.dir/combined.filtersummary.txt'.format(scriptloc,os.getcwd(),os.getcwd())
    P.run(statement)
    

if __name__ == "__main__":
    if sys.argv[1] == "plot":
        pipeline_printout_graph("test.pdf", "pdf", [full], no_key_legend=True,
                                size=(4, 4),
                                user_colour_scheme = {"colour_scheme_index": 1})
    else:
        sys.exit(P.main(sys.argv))



        
