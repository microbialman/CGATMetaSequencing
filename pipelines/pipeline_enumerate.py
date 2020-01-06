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
from ruffus.combinatorics import *
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
       ["%s/pipeline.yml" % __file__[:-len(".py")],
       "../pipeline.yml",
       "pipeline.yml" ] )
PARAMS = P.PARAMS

FEATURES = P.as_list(PARAMS.get("General_feature_list"))
FEATUREPAIRS = P.as_list(PARAMS.get("General_feature_pairs"))
FEATUREPAIRS = ["{}_BY_{}".format(x.split(":")[0],x.split(":")[1]) for x in FEATUREPAIRS]
ALLFEATURES = FEATURES+FEATUREPAIRS

from pipeline_assembly import PipelineAssembly
from pipeline_enumerate import PipelineEnumerate
from pipeline_filter import PipelineFilter

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
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"contig_databases.dir/\1.contigs.bowtie.log")
def makeBowtieDbs(infile,outfile):
    filename=re.match(r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)",infile).group(1)
    filemap=PipelineEnumerate.enumerateMapper(filename,PARAMS)
    #call to bowtie2-build
    job_memory = str(PARAMS["BowtieDB_memory"])+"G"
    job_threads = int(PARAMS["BowtieDB_threads"])
    statement = PipelineEnumerate.buildBowtieDB(filemap.contigpath,outfile.replace(".log",""),PARAMS)
    statement += ' && echo "Made file {}." > {}'.format(outfile.replace(".log",""),outfile)
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
    #only skip if there was a failure in a previous run at the bowtie step
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

######################################################
# Compress the bowtie databases to save space
######################################################
@follows(mapSamples)
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"contig_databases.dir/\1.contigs.bowtie.1.bt2l.gz")
def gzipDatabases(infile,outfile):
    filename=re.match(r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)",infile).group(1)
    statement = "gzip contig_databases.dir/{}.*.bt2l".format(filename)
    P.run(statement)
    
########################################################################
# Count reads mapping to each ORF using featureCounts and GTF files
########################################################################
@follows(mkdir("orf_counts.dir"))
@follows(gzipDatabases)
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"orf_counts.dir/\1.tsv.gz")
def countOrfs(infile,outfile):
    feat = "gene_id"
    outfile = outfile.replace(".gz","")
    filename=re.match(r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)",infile).group(1)
    filemap=PipelineEnumerate.enumerateMapper(filename,PARAMS)
    mapping="sample_mappings.dir/{}/{}.mapped.bam".format(filemap.samplename,filemap.samplename)
    #get pairedness
    paired = True
    num_paired = subprocess.check_output(["samtools","view","-c","-f 1","{}".format(mapping)]).decode(sys.stdout.encoding)
    if int(num_paired.strip("\n")) == 0:
        paired = False
    #generate counts per orf across all samples
    job_threads = int(PARAMS["featureCounts_threads"])
    job_memory = str(PARAMS["featureCounts_memory"])+"G"
    scratchgtf = PARAMS["featureCounts_tmp"]+filename+".gtf"
    statementlist = []
    statementlist.append("zcat {} > {}".format(filemap.shortgtfpath,scratchgtf))
    statementlist.append(PipelineEnumerate.countFeatures(feat,scratchgtf,paired,outfile,mapping,PARAMS))
    statementlist.append("gzip {}".format(outfile))
    statementlist.append("rm -rf {}".format(scratchgtf))
    statement = " && ".join(statementlist)
    P.run(statement)

################################################
# TPM normalise if enabled
################################################
@follows(countOrfs)
@active_if(PARAMS["General_tpm"]=="true")
@transform(countOrfs,regex(r"orf_counts.dir/(\S+).tsv.gz"),r"orf_counts.dir/\1.tpm.tsv.gz")
def makeTPM(infile,outfile):
    statement = "python {}scripts/makeTPM.py --orf_counts {} --output {}".format(
        os.path.dirname(__file__).rstrip("pipelines"),infile,outfile.replace(".gz",""))
    statement += " && gzip {}".format(outfile.replace(".gz",""))
    P.run(statement)

############################################################
# Counting taxa and function features from gene_id counts
############################################################
@follows(makeTPM)
@follows(mkdir("annotation_counts.dir"))
@follows(mkdir("annotation_counts.dir/logs"))
@mkdir(FEATURES,regex(r"(\S+)"),r"annotation_counts.dir/\1")
@transform(countOrfs,regex(r"orf_counts.dir/(\S+).tsv"),r"annotation_counts.dir/logs/\1.log")
def countFeatures(infile,outfile):
    filename=re.search("orf_counts.dir/(\S+).tsv",infile).group(1)
    filemap=PipelineEnumerate.enumerateMapper(filename,PARAMS)
    #generate counts for other features from ORF counts and full GTF
    job_threads = int(PARAMS["featureCounts_threads_otherfeats"])
    job_memory = str(PARAMS["featureCounts_memory_otherfeats"])+"G"
    statement = "python {}scripts/countFeat.py --orf_counts {} --features {} --multimethod {} --gtf {} --outdir annotation_counts.dir/ --logfile {}".format(
        os.path.dirname(__file__).rstrip("pipelines"),infile,",".join(FEATURES),PARAMS["General_multimethod"],filemap.gtfpath,outfile)
    #also count the tpm counts if enabled
    if PARAMS["General_tpm"]=="true":
        statement += " && python {}scripts/countFeat.py --orf_counts {} --features {} --multimethod {} --gtf {} --outdir annotation_counts.dir/ --logfile {}".format(os.path.dirname(__file__).rstrip("pipelines"),infile.replace(".tsv",".tpm.tsv"),",".join(FEATURES),PARAMS["General_multimethod"],filemap.gtfpath,outfile.replace(".log",".tpm.log"))
    P.run(statement)


############################################################
# Counting paired features from gene_id counts
############################################################
@follows(countFeatures)
@mkdir(FEATUREPAIRS,regex(r"(\S+)"),r"annotation_counts.dir/\1")
@transform(countOrfs,regex(r"orf_counts.dir/(\S+).tsv"),r"annotation_counts.dir/logs/\1_paired.log")
def countPairedFeatures(infile,outfile):
    filename=re.search("orf_counts.dir/(\S+).tsv",infile).group(1)
    filemap=PipelineEnumerate.enumerateMapper(filename,PARAMS)
    #generate counts for other features from ORF counts and full GTF
    job_threads = int(PARAMS["featureCounts_threads_otherfeats"])
    job_memory = str(PARAMS["featureCounts_memory_otherfeats"])+"G"
    statement = "python {}scripts/countFeatPairs.py --orf_counts {} --feature_pairs {} --multimethod {} --gtf {} --outdir annotation_counts.dir/ --logfile {}".format(
        os.path.dirname(__file__).rstrip("pipelines"),infile,",".join(FEATUREPAIRS),PARAMS["General_multimethod"],filemap.gtfpath,outfile)
    #count the tpm counts if enabled
    if PARAMS["General_tpm"]=="true":
        statement += " && python {}scripts/countFeatPairs.py --orf_counts {} --feature_pairs {} --multimethod {} --gtf {} --outdir annotation_counts.dir/ --logfile {}".format(
            os.path.dirname(__file__).rstrip("pipelines"),infile.replace(".tsv",".tpm.tsv"),",".join(FEATUREPAIRS),PARAMS["General_multimethod"],filemap.gtfpath,outfile.replace(".log",".tpm.log"))
    P.run(statement)
    
    
###########################################################
# Pool the counts from all the samples into combined files
###########################################################
@follows(countPairedFeatures)
@follows(mkdir("combined_counts.dir"))
@originate(["combined_counts.dir/{}_combined_counts.tsv".format(x) for x in ALLFEATURES])
def combineCounts(outfile):
    feat = re.search("combined_counts.dir/(\S+)_combined_counts\.tsv",outfile).group(1)
    statement = "python {}scripts/combineCounts.py --feature {} --countdir {} --outfile {} --tpm false".format(os.path.dirname(__file__).rstrip("pipelines"),feat,"annotation_counts.dir",outfile)
    if PARAMS["General_tpm"]=="true":    
        statement += " && python {}scripts/combineCounts.py --feature {} --countdir {} --outfile {} --tpm true".format(os.path.dirname(__file__).rstrip("pipelines"),feat,"annotation_counts.dir",outfile.replace("tsv","tpm.tsv"))
    job_threads = int(PARAMS["combineCounts_threads"])
    job_memory = str(PARAMS["combineCounts_memory"])+"G"
    P.run(statement)    

@follows(combineCounts)
def full():
    pass

#Make report (uses raw counts not tpms)
@follows(mkdir("report.dir"))
def build_report():
    job_memory = str(PARAMS["report_memory"])+"G"
    job_threads = int(PARAMS["report_threads"])
    scriptloc = "/".join(os.path.dirname(sys.argv[0]).split("/")[0:-1])+"/scripts/enumeration_report.Rmd"
    statement = 'R -e "rmarkdown::render(\'{}\',output_file=\'{}/report.dir/enumeration_report.html\')" --args {} {}'.format(
        scriptloc,os.getcwd(),",".join(FEATURES),os.getcwd()+"/combined_counts.dir")
    P.run(statement)

if __name__ == "__main__":
    if sys.argv[1] == "plot":
        pipeline_printout_graph("test.pdf", "pdf", [full], no_key_legend=True,
                                size=(4, 4),
                                user_colour_scheme = {"colour_scheme_index": 1})
    else:
        sys.exit(P.main(sys.argv))
