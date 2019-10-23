"""
=============================
Metagenome annotation pipeline
=============================

:Author: Matt Jackson
:Tags: Python

Takes contigs from assembled metagenomic data, identfies ORFs and carrys out functional and taxonomic annotation.

Overview
========

This pipeline uses Prodigal to detect ORFs in metagenomic assemblies.
Eggnog mapper is then used to assign functional annotations to each ORF.

Taxonomy can then be assigned either at the ORF or contig level.
ORF level assignment uses Diamond to align to the NCBI reference (using amino-acid sequence) and then MEGAN lowest-common ancestor approach on alignments.
Contig level assignment uses Kraken to assign taxonomy (using nucleotide sequence) and then assigns each ORF in the contig the contigs taxonomy.

Trade-off is that contig level has the potential to assign more accurate taxonomy but will be more highly influenced by chimeric contigs.
ORF level taxonomic assignments will be less influenced by chimeric contigs but might be less accurate in the case of horizontal gene transfer events.
ORF level assignments based on aa sequence will also suffer from reduced resolution comapred to DNA alignmnets and their shorter length relative to contigs.

  
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

Requirements
------------

On top of the default CGAT setup and the MetaSequencing pipeline, annotation requires the following
software to be in the path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Tested Version*   |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
| Prodigal           | 2.6.3             | Detect ORFs                                    |
+--------------------+-------------------+------------------------------------------------+
| eggnog-mapper      | 1.0.3             | Functional annotation of ORFs                  |  
+--------------------+-------------------+------------------------------------------------+
| Diamond            | 0.9.22            | Alignment of ORFs vs NCBI database             |
+--------------------+-------------------+------------------------------------------------+
| MEGAN CE           | 6.11.6            | Taxonomic assignment from Diamond alignments   |
+--------------------+-------------------+------------------------------------------------+
| Kraken2            | 2.0.8-beta        | Taxonomic assignment from contigs              |
+--------------------+-------------------+------------------------------------------------+
| Taxonkit           | 0.3.0             | Translation of NCBI taxids to taxonomic names  |
+--------------------+-------------------+------------------------------------------------+

Also requires suitable reference databases indexed with the appropriate software versions.

Pipeline output
===============

The main output is a set of functionally and taxonomically annotated contigs.

Glossary
========

.. glossary::

Code
====

"""

#load modules
from ruffus import *
import os
import sys, re
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

from pipeline_assembly import PipelineAssembly
from pipeline_annotate import PipelineAnnotate

#get all files within the directory to process
SEQUENCEFILES = ("*.fasta", "*.fasta.gz", "*.fasta.1.gz", "*.fasta.1",
                 "*.fna", "*.fna.gz", "*.fna.1.gz", "*.fna.1",
                 "*.fa", "*.fa.gz", "*.fa.1.gz", "*.fa.1", 
                 "*.fastq", "*.fastq.gz", "*.fastq.1.gz","*.fastq.1")

SEQUENCEFILES_REGEX = regex(
    r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)")


####################################################
# Find ORFs using prodigal
####################################################
@follows(mkdir("orfs.dir"))
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"orfs.dir/\1.orf_peptides.gz")
def detectOrfs(infile,outfile):
    statementlist=[]
    #set job memory and threads
    job_memory = str(PARAMS["Prodigal_memory"])+"G"
    job_threads = int(PARAMS["Prodigal_threads"])
    #command to generate index files
    seqdat = PipelineAssembly.SequencingData(infile)
    #generate outfile without gz
    outfile=outfile.replace(".gz","")
    #ensure input is FASTA
    if seqdat.paired == True:
        print("Cannot detect ORFs from paired-end reads.")
    else:
        if seqdat.fileformat == "fastq":
            statementlist.append("reformat.sh in={} out={}".format(infile,"orfs.dir/"+seqdat.cleanname+".fa"))
            infile = "orfs.dir/"+seqdat.cleanname+".fa"
        #generate the call to prodigal
        statementlist.append(PipelineAnnotate.runProdigal(infile,outfile,PARAMS))
        #remove the temp FASTA if created
        if seqdat.fileformat == "fastq":
            statementlist.append("rm {}".format("orfs.dir/"+seqdat.cleanname+".fa"))
        #compress the outputs
        statementlist.append("gzip {}".format(outfile))
        statementlist.append("gzip {}".format(outfile.replace("peptides","positions")))
        statement = " && ".join(statementlist)
        P.run(statement)

######################################################################################################################
# Find seed orthologs of ORFs using eggnog-mapper with chunk splitting
######################################################################################################################
@follows(detectOrfs)
@follows(mkdir("functional_annotations.dir"))
@follows(mkdir("functional_annotations.dir/emapper_chunks"))
@active_if(PARAMS["Eggnogmapper_skip"]=="false")
@subdivide(detectOrfs,regex(r"orfs.dir/(\S+).orf_peptides.gz"),r"functional_annotations.dir/emapper_chunks/\1.*.chunk.log",r"functional_annotations.dir/emapper_chunks/\1")
def splitFasta(infile,outfiles,outfileroot):
    statement = "python {}/fastaToChunks.py --input {} --output_prefix {} --chunk_size {}".format(os.path.dirname(os.path.abspath(__file__)).replace("pipelines","scripts"),infile,os.getcwd()+"/"+outfileroot,PARAMS["Eggnogmapper_chunksize"])
    P.run(statement)

@follows(splitFasta)
@active_if(PARAMS["Eggnogmapper_skip"]=="false")
@transform(splitFasta,regex(r"functional_annotations.dir/emapper_chunks/(\S+).chunk.log"),r"functional_annotations.dir/emapper_chunks/\1.chunk.emapper.seed_orthologs.log")
def functionalAnnotSeed(infile,outfile):
    job_memory = str(PARAMS["Eggnogmapper_memory"])+"G"
    job_threads = int(PARAMS["Eggnogmapper_threads"])
    #generate call to eggnog-mapper
    #requires older version of diamond to use the eggnog mapper databases
    statement = "module load bio/diamond/0.8.22 && "
    statement += PipelineAnnotate.runEggmapSeed(infile.replace(".log",""),infile.replace(".log",""),PARAMS)
    statement += ' && echo "Made file {}." > {}'.format(outfile.replace(".log",""),outfile)
    P.run(statement)

###############################################################
# Functional annotation of the seeds
###############################################################
@follows(functionalAnnotSeed)
@active_if(PARAMS["Eggnogmapper_skip"]=="false")
@transform(functionalAnnotSeed,regex(r"functional_annotations.dir/emapper_chunks/(\S+).emapper.seed_orthologs.log"),r"functional_annotations.dir/emapper_chunks/\1.emapper.annotations.log")
def functionalAnnotChunks(infile,outfile):
    job_memory = str(PARAMS["Eggnogmapper_memory_annot"])+"G"
    job_threads = int(str(PARAMS["Eggnogmapper_threads_annot"]))
    statement=[]
    if PARAMS["Eggnogmapper_scratch"] == "true":
        #copy the db into fast local SSD
        statement.append("cp {}eggnog.db $SCRATCH_DIR/eggnog.db".format(PARAMS["Eggnogmapper_eggdata"]))
        datadir="$SCRATCH_DIR"
    else:
        datadir=PARAMS["Eggnogmapper_eggdata"]
    #get annotation from seeds
    statement.append(PipelineAnnotate.runEggmapAnnot(infile.replace(".log",""),outfile.replace(".emapper.annotations.log",""),PARAMS,datadir))
    statement.append('echo "Made file {}." > {}'.format(outfile.replace(".log",""),outfile))
    statement = " && ".join(statement)
    #run the annotation step
    P.run(statement)

###################################################
# Merge the functional annotations
##################################################
@follows(functionalAnnotChunks)
@active_if(PARAMS["Eggnogmapper_skip"]=="false")
@collate(functionalAnnotChunks,regex(r"functional_annotations\.dir/emapper_chunks/(\S+)\.[0-9]+\.chunk\.emapper\.annotations.log"),r"functional_annotations.dir/\1.functional.annotations.gz")
def functionalAnnot(infiles,outfile):
    infiles=[x.replace(".log","") for x in infiles]
    filename=re.search("(\S+).contigs.(\S+).chunk.emapper.annotations",infiles[0]).group(1)
    statementlist = ["cat {} >> {}".format(" ".join(infiles),outfile.replace(".gz",""))]
    statementlist.append("gzip {}".format(outfile.replace(".gz","")))
    #clean up intermediate files
    statementlist.append('find {}.* ! -name "*.log" -delete'.format(os.path.join(os.getcwd(),filename)))
    statement = " && ".join(statementlist)
    P.run(statement)

##################################################
# Taxonomic alignment of ORFs using DIAMOND
#################################################
@follows(functionalAnnot)
@follows(mkdir("taxonomic_annotations.dir"))
@active_if(PARAMS["General_tax_method"]=="diamond")
@transform(detectOrfs,regex(r"orfs.dir/(\S+).orf_peptides"),r"taxonomic_annotations.dir/\1.daa")
def diamondAlignORF(infile,outfile):
    #set memory and threads
    job_memory = str(PARAMS["Diamond_memory"])+"G"
    job_threads = int(PARAMS["Diamond_threads"])
    #generate call to diamond
    statement = PipelineAnnotate.runDiamond(infile,outfile,PARAMS)
    P.run(statement)

########################################################################################################
# Get taxanomic annotation (and optionally kegg functions) from DIAMOND alignment using MEGAN blast2lca
########################################################################################################
@follows(diamondAlignORF)
@active_if(PARAMS["General_tax_method"]=="diamond")
@transform(diamondAlignORF,regex(r"taxonomic_annotations.dir/(\S+).daa"),r"taxonomic_annotations.dir/\1.taxonomic.annotations")
def meganAnnot(infile,outfile):
    job_memory = str(PARAMS["Blast2lca_memory"])+"G"
    job_threads = int(PARAMS["Blast2lca_threads"])
    #generate call to blast2lca
    statement = PipelineAnnotate.runBlast2Lca(infile,outfile,PARAMS)
    P.run(statement)

###############################################################################################
# Get taxonomic assignment using Kraken2 and translate NCBI taxids to mpa style names
###############################################################################################
@follows(functionalAnnot)
@follows(mkdir("taxonomic_annotations.dir"))
@active_if(PARAMS["General_tax_method"]=="kraken")
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"taxonomic_annotations.dir/\1.kraken_translated.gz")
def krakenAlignContig(infile,outfile):
    #set memory and threads
    job_memory = str(PARAMS["Kraken_memory"])+"G"
    job_threads = int(PARAMS["Kraken_threads"])
    #generate call to diamond
    statementlist = [PipelineAnnotate.runKraken(infile,outfile.replace(".gz",""),PARAMS,"{}scripts/translateKraken2.py".format(os.path.dirname(__file__).rstrip("pipelines")))]
    #compress outputs
    statementlist.append("gzip {}".format(outfile.replace("translated.gz","out")))
    statementlist.append("gzip {}".format(outfile.replace(".gz","")))
    statement= " && ".join(statementlist)
    P.run(statement)

######################################################################
# Convert Kraken2 contig level output to match ORF level MEGAN output
######################################################################
@follows(krakenAlignContig)
@active_if(PARAMS["General_tax_method"]=="kraken")
@transform(krakenAlignContig,regex(r"taxonomic_annotations.dir/(\S+).kraken_translated.gz"),r"taxonomic_annotations.dir/\1.taxonomic.annotations.gz")
def krakenFormatAnnot(infile,outfile):
    #get the orf names for the sample
    sampleName=re.search("taxonomic_annotations.dir/(\S+).kraken_translated.gz",infile).group(1)
    orffile="orfs.dir/{}.orf_peptides.gz".format(sampleName)
    statement="python {}scripts/krakenFormat.py --orfs {} --contig-taxonomy {} --orf-taxonomy-output {}".format(
        os.path.dirname(__file__).rstrip("pipelines"),orffile,infile,outfile.replace(".gz",""))
    statement += " && gzip {}".format(outfile.replace(".gz",""))
    P.run(statement)


################################################
# Generate GTF summarising ORF annotations
################################################
@follows(meganAnnot)
@follows(krakenFormatAnnot)
@follows(mkdir("combined_annotations.dir"))
@transform(detectOrfs,regex(r"orfs.dir/(\S+).orf_peptides.gz"),r"combined_annotations.dir/\1.orf_annotations.gtf.gz")
def mergeAnnotations(infile,outfile):
    sampleName=re.search("orfs.dir/(\S+).orf_peptides.gz",infile).group(1)
    func="functional_annotations.dir/{}.functional.annotations.gz".format(sampleName)
    tax="taxonomic_annotations.dir/{}.taxonomic.annotations.gz".format(sampleName)
    gzless=outfile.replace(".gz","")
    job_memory = str(PARAMS["Merge_memory"])+"G"
    statementlist = ["python {}scripts/makeGtf.py --orfs {} --functions {} --taxonomy {} --output {} --output-short {}".format(os.path.dirname(__file__).rstrip("pipelines"), infile, func, tax, gzless, gzless.replace(".gtf","_short.gtf"))]
    statementlist.append("gzip {}".format(gzless.replace(".gtf","_short.gtf")))
    statementlist.append("gzip {}".format(gzless))
    statement = " && ".join(statementlist)
    P.run(statement)

@follows(mergeAnnotations)
def full():
    pass


###################################################
# Summarise annotations
###################################################
@follows(mkdir("report.dir"))
@originate("report.dir/gtf_summary.tsv")
def summariseGTFs(outfile):
    job_memory = str(PARAMS["Merge_memory"])+"G"
    statement = "python {}scripts/annotationSummaryTable.py --gtf-dir combined_annotations.dir --annot-output {} --orf-output {}".format(
        os.path.dirname(__file__).rstrip("pipelines"),outfile,"report.dir/orf_summary.tsv")
    P.run(statement)

@follows(summariseGTFs)
def build_report():
    job_memory = str(PARAMS["Merge_memory"])+"G"
    scriptloc = "/".join(os.path.dirname(sys.argv[0]).split("/")[0:-1])+"/scripts/annotation_report.Rmd"
    statement = 'R -e "rmarkdown::render(\'{}\',output_file=\'{}/report.dir/annotation_report.html\')" --args {}/report.dir/gtf_summary.tsv {}/report.dir/orf_summary.tsv'.format(scriptloc,os.getcwd(),os.getcwd(),os.getcwd())
    P.run(statement)


    
if __name__ == "__main__":
    if sys.argv[1] == "plot":
        pipeline_printout_graph("test.pdf", "pdf", [full], no_key_legend=True,
                                size=(4, 4),
                                user_colour_scheme = {"colour_scheme_index": 1})
    else:
        sys.exit(P.main(sys.argv))



        
