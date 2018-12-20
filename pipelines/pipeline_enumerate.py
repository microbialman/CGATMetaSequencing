"""
=============================
Metagenome enumeration pipeline
=============================

:Author: Matt Jackson
:Tags: Python

Takes filtered reads from metagenomic data and enumerates annotations per sample from pre-defined ORFs within the data.

Overview
========

This pipeline maps the filtered reads from pipeline_filter to the contigs produced by pipeline_assembly.
The resulting mapping is then used to count to reads mapping to the annotated ORFs produced by pipeline_annotated.
Count tables are then generated at the level of each annotation. 

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

The pre-filtered reads to be aligned are imported by placing files are linking to files in the
:term:`working directory`. Directories for the contigs and ORF annotations must be provided in the pipeline.yml file.
The script assumes filenames match across these directories and have the default suffixes of the other pipelines.

The default file format assumes the following convention, with one file/paired file per sample to be enumerated:

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

On top of the default CGAT setup and the MetaSequencing pipeline, the pipeline requires the following
software to be in the path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Tested Version*   |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|bowtie2             | 2.3.0             | Mapping of reads to contigs                    |
+--------------------+-------------------+------------------------------------------------+
|featureCounts       |                   | Counting of reads per ORF annotation           |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The main outputs are tables of counts of functional and taxonomic annotations for each sample.

Glossary
========

.. glossary::

Code
====

"""
#load modules
from ruffus import *
import os, re
import sys, glob
import subprocess


###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
# load options from the config file
import cgatcore.pipeline as P
P.get_parameters(
       ["%s.yml" % __file__[:-len(".py")],
       "../pipeline.yml",
       "pipeline.yml" ] )
PARAMS = P.PARAMS

import PipelineAssembly
import PipelineEnumerate
import PipelineFilter

#get all files within the directory to process
SEQUENCEFILES = ("*.fasta", "*.fasta.gz", "*.fasta.1.gz", "*.fasta.1",
                 "*.fna", "*.fna.gz", "*.fna.1.gz", "*.fna.1",
                 "*.fa", "*.fa.gz", "*.fa.1.gz", "*.fa.1", 
                 "*.fastq", "*.fastq.gz", "*.fastq.1.gz","*.fastq.1")

SEQUENCEFILES_REGEX = regex(
    r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)")

#############################################################
# Make a bowtie indexed database for each sample's contigs
#############################################################
@follows(mkdir("contig_databases.dir"))
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"contig_databases.dir/\1.contigs.bowtie.1.bt2l")
def makeBowtieDbs(infile,outfile):
    filename=re.match(r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)",infile).group(1)
    filemap=PipelineEnumerate.enumerateMapper(filename,PARAMS)
    #call to bowtie2-build
    job_memory = str(PARAMS["BowtieDB_memory"])+"G"
    job_threads = int(PARAMS["BowtieDB_threads"])
    statement = PipelineEnumerate.buildBowtieDB(filemap.contigpath,outfile.replace(".1.bt2l",""),PARAMS)
    P.run(statement)
    

#######################################################
# Map each samples reads against the contig database
#######################################################
@follows(makeBowtieDbs)
@follows(mkdir("sample_mappings.dir"))
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"sample_mappings.dir/\1/\1.mapped.bam")
def mapSamples(infile,outfile):
    filename=re.match(r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)",infile).group(1)
    filemap=PipelineEnumerate.enumerateMapper(filename,PARAMS)
    #get the mapping DB
    bowtiedb="contig_databases.dir/{}.contigs.bowtie".format(filemap.samplename)
    job_threads = int(PARAMS["Bowtie_threads"])
    job_memory = str(PARAMS["Bowtie_memory"])+"G"
    seqdat = PipelineAssembly.SequencingData(infile)
    bowtie = PipelineFilter.Bowtie2(seqdat,outfile,PARAMS,bowtiedb)
    #need to reset the working directory in the bowtie function as it is running on files in one directory
    bowtie.indir = ""
    statementlist = []
    #remove all comments from read names in files (trimming can add comments making non-matching pairs)
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
'''
#################################################################
# Chunk the GTF
#################################################################
@follows(mapSamples)
@follows(mkdir("chunk_gtf.dir"))
@split(PARAMS["General_gtf_file"],"chunk_gtf.dir/*.gtf")
def chunkGTF(infile,outfiles):
    statement = 'split -n {} {} {} && for f in chunk_gtf.dir/*; do mv "$f" "$f.gtf"; done'.format(PARAMS["featureCounts_chunk"],infile,"chunk_gtf.dir/chunk_")
    P.run(statement)

#################################################################
# Count reads mapping to each ORF using featureCounts
#################################################################
@follows(chunkGTF)
@follows(mkdir("chunk_counts.dir"))
@transform(chunkGTF,regex(r"chunk_gtf.dir/(\S+).gtf"),r"chunk_counts.dir/\1.tsv")
def countOrfs(infile,outfile):
    feat = "gene_id"
    samplemappings = [x for x in glob.iglob('sample_mappings.dir/**/*.bam', recursive=True)]
    #get pairdness from first sample
    paired = True
    num_paired = subprocess.check_output(["samtools","view","-c","-f 1","{}".format(os.getcwd()+"/"+samplemappings[0])]).decode(sys.stdout.encoding)
    if int(num_paired.strip("\n")) == 0:
        paired = False
    #generate counts per orf across all samples
    job_threads = int(PARAMS["featureCounts_threads"])
    job_memory = str(PARAMS["featureCounts_memory"])+"G"
    statement = PipelineEnumerate.countFeatures(feat,infile,paired,outfile,samplemappings,PARAMS)
    P.run(statement)

################################################################
# Combine chunk counts
################################################################
@follows(countOrfs)
@follows(mkdir("feature_counts.dir"))
@merge(countOrfs,"feature_counts.dir/gene_id_counts.txt")
def mergeCounts(infiles,outfile):
    statement="head -1 {} >> {} && awk 'FNR>1' {} >> {}".format(infiles[0],outfile," ".join(infiles),outfile)
    P.run(statement)

################################################################
# Generate formatted ORF counts
################################################################
@follows(countOrfs)
@follows(mkdir("formatted_counts.dir"))
@split(mergeCounts,["formatted_counts.dir/raw_counts/gene_id_raw.tsv","formatted_counts.dir/tpm_counts/gene_id_tpm.tsv"])
def tpmOrfs(infile,outfiles):
    job_memory = str(PARAMS["normalise_memory"])+"G"
    job_threads = int(PARAMS["normalise_threads"])
    statement = "Rscript {}scripts/tpmCounts.R {} {} {}".format(os.path.dirname(__file__).rstrip("pipelines"),infile[0],outfiles[0],outfiles[1])
    P.run(statement)

#########################################################
# Do counting for other features from gene_id counts
#########################################################
@follows(tpmOrfs)
@originate(["formatted_counts.dir/raw_counts/{}_raw.tsv".format(x) for x in PARAMS["General_feature_list"].split(",")])
def countRawFeatures(outfile):
    feat = re.search("formatted_counts.dir/raw_counts/(\S+)_raw.tsv",outfile).group(1)
    #generate counts for other features from ORF counts
    job_threads = int(PARAMS["featureCounts_threads_otherfeats"])
    job_memory = str(PARAMS["featureCounts_memory_otherfeats"])+"G"
    statement = "python {}scripts/countFeat.py --orfinput {} --feature {} --gtf {} --output {}".format(os.path.dirname(__file__).rstrip("pipelines"),
                                                                                                                      "formatted_counts.dir/raw_counts/gene_id_raw.tsv",
                                                                                                                      feat,
                                                                                                                      PARAMS["General_gtf_file_full"],
                                                                                                                      outfile)
    P.run(statement)


#######################################
# Move out log files
######################################
@follows(countRawFeatures)
@follows(countTpmFeatures)
@follows(mkdir("count_logs.dir"))
@originate(["count_logs.dir/gtflog_{}_raw.tsv".format(x) for x in PARAMS["General_feature_list"].split(",")])
def moveLogs():
    statement="mv -f formatted_counts.dir/*/gtflog* count_logs.dir/"
    P.run(statement)

@follows(moveLogs)
def full():
    pass


@follows(mkdir("report.dir"))
def build_report():
    job_memory = str(PARAMS["report_memory"])+"G"
    scriptloc = "/".join(os.path.dirname(sys.argv[0]).split("/")[0:-1])+"/scripts/enumeration_report.Rmd"
    statement = 'R -e "rmarkdown::render(\'{}\',output_file=\'{}/report.dir/enumeration_report.html\')" --args {} {} {}'.format(scriptloc,os.getcwd(),PARAMS["report_readcounts"],os.getcwd()+"/formatted_counts.dir/raw_counts/gene_id_raw.tsv",os.getcwd()+"/formatted_counts.dir/raw_counts")
    P.run(statement)
'''
if __name__ == "__main__":
    if sys.argv[1] == "plot":
        pipeline_printout_graph("test.pdf", "pdf", [full], no_key_legend=True,
                                size=(4, 4),
                                user_colour_scheme = {"colour_scheme_index": 1})
    else:
        sys.exit(P.main(sys.argv))
