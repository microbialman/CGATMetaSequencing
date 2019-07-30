'''
classes and utility functions for pipeline_filter.py

'''

import os
import sys
sys.path.append('..')

from pipeline_assembly import PipelineAssembly

'''
class to build call to SortMeRNA
'''

class SortMeRNA:    
    def __init__(self,seqdat,outfile,params):
        
        self.seqdat = seqdat
        self.outfile = outfile
        self.outdir = os.getcwd()+"/"+os.path.dirname(outfile)
        self.indir = os.getcwd()+"/"
        self.params = params
        self.filelocation = ""
        self.statementlist = ["mkdir -p {}".format(self.outdir)]
        self.rminter = False
        self.checkInterleave()
        self.buildStatement()
        self.deInterleave()
        
    #add statements to generate interleaved files if necessary
    def checkInterleave(self):
        #if single end or already interleaved just use the original file location
        if self.seqdat.paired == False or self.seqdat.interleaved == True:
            self.filelocation = os.getcwd()+"/"+self.seqdat.filepath
        #else generate commands to make temp interleaved reads
        else:
            self.statementlist.append("mkdir -p {}/interleaved".format(self.outdir))
            self.filelocation = self.outdir+"/interleaved/{}".format(self.seqdat.cleanname+"."+self.seqdat.fileformat)
            self.statementlist.append("seqtk mergepe {} {} >{}".format(self.indir+self.seqdat.filepath,self.indir+self.seqdat.pairedname,self.filelocation))

    #make the main call to sortmerna
    def buildStatement(self):
        sortlist = ["sortmerna"]
        sortlist.append(self.refList())
        if self.params["SortMeRNA_paired"] == "out":
            sortlist.append("--paired_out")
        else:
            sortlist.append("--paired_in")
        sortlist.append("--reads {}".format(self.filelocation))
        sortlist.append("--aligned {}".format(self.outdir+"/aligned_"+self.seqdat.cleanname))
        sortlist.append("--other {}".format(self.outdir+"/other_"+self.seqdat.cleanname))
        if self.params["SortMeRNA_fastx"] == "true":
            sortlist.append("--fastx")
        if self.params["SortMeRNA_sam"] == "true":
            sortlist.append("--sam")
        if self.params["SortMeRNA_sq"] == "true":
            sortlist.append("--SQ")
        if self.params["SortMeRNA_blast"] != "false":
            sortlist.append("--blast {}".format(self.params["SortMeRNA_blast"]))
        if self.params["SortMeRNA_log"] == "true":
            sortlist.append("--log")
        if self.params["SortMeRNA_num_alignments"] != "false":
            sortlist.append("--num_alignments {}".format(self.params["SortMeRNA_num_alignments"]))
        if self.params["SortMeRNA_best"] != "false":
            sortlist.append("--best {}".format(self.params["SortMeRNA_best"]))
            sortlist.append("--min_lis {}".format(self.params["SortMeRNA_min_lis"]))
        if self.params["SortMeRNA_print_all_reads"] == "true":
            sortlist.append("--print_all_reads")
        sortlist.append("--match {}".format(self.params["SortMeRNA_match"]))
        sortlist.append("--mismatch {}".format(self.params["SortMeRNA_mismatch"]))
        sortlist.append("--gap_open {}".format(self.params["SortMeRNA_gap_open"]))
        sortlist.append("--gap_ext {}".format(self.params["SortMeRNA_gap_ext"]))
        sortlist.append("-N {}".format(self.params["SortMeRNA_n"]))
        if self.params["SortMeRNA_f"] == "true":
            sortlist.append("-F")
        if self.params["SortMeRNA_r"] == "true":
            sortlist.append("-R")
        sortlist.append("-a {}".format(self.params["SortMeRNA_threads"]))
        sortlist.append("-e {}".format(self.params["SortMeRNA_e"]))
        if self.params["SortMeRNA_memory"] != "false":
            sortlist.append("-m {}".format(str(int(self.params["SortMeRNA_memory"])*900*int(self.params["SortMeRNA_threads"]))))
        if self.params["SortMeRNA_v"] == "true":
            sortlist.append("-v")
        self.statementlist.append(" ".join(sortlist))

    #if input was interleaved remove interleaved temp file (if necessary)
    def deInterleave(self):
        if self.seqdat.paired == False:
            pass
        elif self.seqdat.paired == True:
            if self.seqdat.interleaved == False:
                #remove temp interleave file if one was made
                self.statementlist.append("rm -r {}".format(self.outdir+"/interleaved/"))
                #undo interleaving
                otherf = self.outdir+"/other_"+self.seqdat.cleanname+"."+self.seqdat.fileformat
                pair1 = self.outdir+"/other_"+self.seqdat.filename
                pair2 = self.outdir+"/other_"+self.seqdat.pairedname
                self.statementlist.append("seqtk seq -l0 -1 {} > {}".format(otherf,pair1.strip(".gz")))
                self.statementlist.append("seqtk seq -l0 -2 {} > {}".format(otherf,pair2.strip(".gz")))
                self.statementlist.append("rm {}".format(otherf))
                #compress the outputs if inputs were
                self.statementlist.append("gzip {}".format(pair1.strip(".gz")))
                self.statementlist.append("gzip {}".format(pair2.strip(".gz")))
                    
                
    #abstract out making reference command as it is long 
    def refList(self):
        reffastas = self.params["SortMeRNA_rna_refs"].split(",")
        if self.params["SortMeRNA_rna_index"] != "false":
            refindex = self.params["SortMeRNA_rna_index"].split(",")
        else:
            refindex = reffastas
            refindex = [os.path.basename(i) for i in refindex]
            refindex = [os.getcwd()+"/ref_index.dir/{}-db".format(i) for i in refindex]
        paired = [",".join([x[0],x[1]]) for x in zip(reffastas,refindex)]
        return("--ref {}".format(":".join(paired)))
            

    def build(self):
        return(" && ".join(self.statementlist))
            

'''
Class to call bowtie2 mapping command
'''
class Bowtie2:
    def __init__(self,seqdat,outfile,params,db):
        
        self.seqdat = seqdat
        self.db = db
        self.outfile = outfile
        self.indir = os.path.dirname(self.seqdat.filepath)+"/"
        self.outdir = os.path.dirname(outfile)
        self.params = params

    #function to strip comments from fastq and fasta names that interfere with bowtie paired mode
    #skipped if prep_files is false (used to skip if crash during first bowtie run)
    def cleanNames(self):
        statementlist = []
        wd = self.indir
        statementlist.append("zcat -f {} | awk '{{print $1}}' > {}temp1{}".format(wd+self.seqdat.filename,wd,self.seqdat.filename))
        statementlist.append("rm {} && mv {}temp1{} {}".format(wd+self.seqdat.filename,wd,self.seqdat.filename,wd+self.seqdat.filename.rstrip(".gz")))
        if self.seqdat.paired == True and self.seqdat.interleaved == False:
            statementlist.append("zcat -f {} | awk '{{print $1}}' > {}temp2{}".format(wd+self.seqdat.pairedname,wd,self.seqdat.pairedname))
            statementlist.append("rm {} && mv {}temp2{} {}".format(wd+self.seqdat.pairedname,wd,self.seqdat.pairedname,wd+self.seqdat.pairedname.rstrip(".gz")))
        if self.seqdat.compressed == True:
            statementlist.append("gzip {}".format(wd+self.seqdat.filename.strip(".gz")))
            if self.seqdat.paired == True and self.seqdat.interleaved == False:
                statementlist.append("gzip {}".format(wd+self.seqdat.pairedname.strip(".gz")))
        return(" && ".join(statementlist))

    #main call to bowtie implements most arguments    
    def build(self):
        statementlist = ["bowtie2"]
        statementlist.append("-x {}".format(self.db))
        if self.seqdat.paired == False:
            statementlist.append("-U {}".format(self.seqdat.filepath))
        elif self.seqdat.interleaved == True:
            statementlist.append("--interleaved {}".format(self.seqdat.filepath))
        else:
            statementlist.append("-1 {} -2 {}".format(self.seqdat.filepath,self.indir+self.seqdat.pairedname))
        statementlist.append("-S {}".format(self.outfile.replace("bam","sam")))
        if self.seqdat.fileformat == "fasta":
            statementlist.append("-f")
        else:
            statementlist.append("-q")
        if self.params["Bowtie_phred_type"] == "64":
            statementlist.append("--phred64")
        statementlist.append("--{}".format(self.params["Bowtie_mode"]))
        if self.params["Bowtie_preset"] != "false":
            statementlist.append("--{}".format(self.params["Bowtie_preset"]))
        statementlist.append("-N {}".format(self.params["Bowtie_n"]))
        statementlist.append("-L {}".format(self.params["Bowtie_l"]))
        statementlist.append("-i {}".format(self.params["Bowtie_i"]))
        statementlist.append("--n-ceil {}".format(self.params["Bowtie_n_ceil"]))
        statementlist.append("--dpad {}".format(self.params["Bowtie_dpad"]))
        statementlist.append("--gbar {}".format(self.params["Bowtie_gbar"]))
        if self.params["Bowtie_ignore_quals"] == "true":
            statementlist.append("--ignore-quals")
        if self.params["Bowtie_ignore_quals"] == "true":
            statementlist.append["--ignore-quals"]
        if self.params["Bowtie_nofw"] == "true":
            statementlist.append["--nofw"]
        if self.params["Bowtie_norc"] == "true":
            statementlist.append["--norc"]
        if self.params["Bowtie_no_1mm_upfront"] == "true":
            statementlist.append["--no-1mm-upfront"]
        statementlist.append("--ma {}".format(self.params["Bowtie_ma"]))
        statementlist.append("--mp {}".format(self.params["Bowtie_mp"]))
        statementlist.append("--np {}".format(self.params["Bowtie_np"]))
        statementlist.append("--rdg {}".format(self.params["Bowtie_rdg"]))
        statementlist.append("--rfg {}".format(self.params["Bowtie_rfg"]))
        statementlist.append("--score-min {}".format(self.params["Bowtie_score_min"]))
        if self.params["Bowtie_reporting"] != "default":
            statementlist.append("--{}".format(self.params["Bowtie_reporting"]))
        statementlist.append("-D {}".format(self.params["Bowtie_d"]))
        statementlist.append("-R {}".format(self.params["Bowtie_r"]))
        statementlist.append("--minins {}".format(self.params["Bowtie_minins"]))
        statementlist.append("--maxins {}".format(self.params["Bowtie_maxins"]))
        if self.params["Bowtie_no_mixed"] == "true":
            statementlist.append["--no-mixed"]
        if self.params["Bowtie_no_discordant"] == "true":
            statementlist.append["--no-discordant"]
        if self.params["Bowtie_no_dovetail"] == "true":
            statementlist.append["--no-dovetail"]
        if self.params["Bowtie_no_contain"] == "true":
            statementlist.append["--no-contain"]
        if self.params["Bowtie_no_overlap"] == "true":
            statementlist.append["--no-overlap"]
        statementlist.append("-p {}".format(self.params["Bowtie_threads"]))
        
        return(" ".join(statementlist))


class FilterFromBam:

    def __init__(self,infile,outfile,seqdat,params):
        self.params = params
        self.seqdat = seqdat
        self.infile = infile
        self.outfile = outfile
        self.filtered = infile.strip(".mapped.bam")+".unmapped_reads.bam"
        self.pairsort = infile.strip(".mapped.bam")+".sorted_unmapped_reads.bam"
        self.outdir = os.path.dirname(outfile)
        self.statementlist = []
        self.unmapped()
        self.convertBam()

    def unmapped(self):
        if self.seqdat.paired == True:
            fflag,Fflag = self.params["Filtering_paired_pos"],self.params["Filtering_paired_neg"]
        else:
            fflag,Fflag = self.params["Filtering_un_pos"],self.params["Filtering_un_neg"]
        unmap = "samtools view -b -f {} -F {} {} >{}".format(fflag,Fflag,self.infile,self.filtered)
        self.statementlist.append(unmap)

    def convertBam(self):
        #if paired reads ensure they are sorted into pairs before conversion to FASTX
        if self.seqdat.paired == True:
            self.statementlist.append("samtools sort -n {} -o {}".format(
                self.filtered,self.pairsort))
            self.statementlist.append("rm {}".format(self.filtered))
            if self.seqdat.interleaved == False:
                self.statementlist.append("bedtools bamtofastq -i {} -fq {} -fq2 {}".format(
                self.pairsort, self.outfile.rstrip(".gz"), self.outdir+"/hostfiltered_"+self.seqdat.pairedname.rstrip(".gz")))
            else:
               self.statementlist.append("bedtools bamtofastq -i {} -fq {}".format(
                self.pairsort, self.outfile.rstrip(".gz"))) 
        else:
            self.statementlist.append("bedtools bamtofastq -i {} -fq {}".format(
                self.filtered, self.outfile.rstrip(".gz")))
        #compress
        self.statementlist.append("gzip {}".format(self.outfile.rstrip(".gz")))
        if self.seqdat.paired == True and self.seqdat.interleaved == False:
            self.statementlist.append("gzip {}".format(
                self.outdir+"/hostfiltered_"+self.seqdat.pairedname.rstrip(".gz")))
        #remove bam files
        self.statementlist.append("rm "+self.outdir+"/*.bam")

    def build(self):
        return(" && ".join(self.statementlist))
    

#function to clean up the directory structures at the end
def CleanUp(seqdat,outfile,params):
    statementlist = []
    rnadir = os.getcwd()+"/rrna_filter_out.dir/"
    gendir = os.getcwd()+"/genome_filter_out.dir/"
    outdir = os.path.dirname(outfile)
    #if only rRNA filtering symlink those files
    if params["General_rrna_filter"] == "true" and params["General_host_filter"] == "false":
        statementlist.append("ln -s {} {}".format(rnadir+seqdat.cleanname+"/other_"+seqdat.filename,outfile))
        if seqdat.paired == True and seqdat.interleaved == False:
            statementlist.append("ln -s {} {}".format(rnadir+seqdat.cleanname+"/other_"+seqdat.pairedname,outdir+"/filtered-"+seqdat.pairedname))
    #else symlink genome filtered files
    else:
        statementlist.append("ln -s {} {}".format(gendir+seqdat.cleanname+"/hostfiltered_"+seqdat.filename,outfile))
        if seqdat.paired == True and seqdat.interleaved == False:
            statementlist.append("ln -s {} {}".format(gendir+seqdat.cleanname+"/hostfiltered_"+seqdat.pairedname,outdir+"/filtered-"+seqdat.pairedname))
    return(" && ".join(statementlist))


#function to summarise reads filtered at each step
def CountReads(infile,outfile,params):
    def counter(seqfile,outfile):
        sdat=PipelineAssembly.SequencingData(seqfile)
        div=4
        if sdat.fileformat == "fasta":
            div=2
        return("zcat {} | wc -l | awk 'BEGIN {{ORS=\"\"}}; END {{x=$1/{}; print \"\\t\"x}}' >> {}".format(seqfile,div,outfile))
        
    original = PipelineAssembly.SequencingData(infile)
    call=['printf "File\\tInput\\tPost_rRNA_Filtering\\tPost_Genome_Filtering\\n{}" > {}'.format(original.cleanname,outfile)]
    ocount = counter(infile,outfile)
    rcount = 'printf "\\tNA" >> {}'.format(outfile)
    gcount = 'printf "\\tNA" >> {}'.format(outfile)
    rnadir = os.getcwd()+"/rrna_filter_out.dir/"
    gendir = os.getcwd()+"/genome_filter_out.dir/"
    if params["General_rrna_filter"] == "true": 
        rcount = counter(rnadir+original.cleanname+"/other_"+original.filename,outfile)
    if params["General_host_filter"] == "true":
        gcount = counter(gendir+original.cleanname+"/hostfiltered_"+original.filename,outfile)
    call.append(ocount)
    call.append(rcount)
    call.append(gcount)
    call.append('printf "\\n" >> {}'.format(outfile))
    return(" && ".join(call))
        
