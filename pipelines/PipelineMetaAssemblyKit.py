'''classes and utility functions for pipeline_MetaAssemblyKit.py

Different assembly tools will use different inputs. Some can take
fasta files whereas others will take fastq and in either case can be
paired-end (in the same or different files) or single end. Output formats
and locations will also differ slightly.

'''


import os
import re
import glob
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P

'''
object for sequencing data
    
On initialisation:
checks the format from the name
checks for compression from name
check for paired-endness

Additional commands:
counts the number of reads in the file
'''    
class SequencingData:
    def __init__(self,infile):
        self.fileformat =  None
        self.paired = False
        self.interleaved = False
        self.compressed = False
        self.pairedname = ""
        self.filepath = infile
        self.filename = os.path.basename(infile)
        self.openfile = None
        self.readcount = None
        self.cleanname = None
        self.head = []
        '''check file can be opened on init & capture header (first 5 lines used for interleave and format checks)'''
        try:
            self.openfile = IOTools.openFile(self.filepath)
        except FileNotFoundError as e:
            raise Exception("cannot open file {}".format(self.filepath)) from e
        self.head = [self.openfile.readline().rstrip("\n") for x in range(5)]
        '''autocheck file format and pairedness, read count must be specified seperately'''
        self.getFormat()
        self.isPaired()
        self.cleanName()
        
        
    '''check it is fasta or fastq and if compressed'''    
    def getFormat(self):
        extensions=("fna","fa","fasta","fastq")
        for i in extensions:
            assert not(self.filepath.endswith((".2",".2.gz"))), "Read 2 file provided ({}) please use read 1 file".format(self.filename)
            if self.filepath.endswith((i,i+".1.gz",i+".gz",i+".1")):
                if i == "fastq":
                    self.fileformat=i
                else:
                    self.fileformat="fasta"
        if self.filepath.find(".gz")!=-1:
            self.compressed = True
        assert self.fileformat,"file {} is not of the correct format (fasta or fastq).".format(self.filename)
        if self.fileformat == "fasta":
            assert self.head[0][0] == ">", "invalid header on first line for fasta format"
        else:
            assert self.head[0][0] == "@", "invalid header on first line for fastq format"

            
    '''check if paired and if containts interleaved pairs or matching files'''
    def isPaired(self):
        if self.filepath.endswith((".1",".1.gz")):
            paired_name = self.filepath.replace(".1",".2")
            assert len(glob.glob(paired_name)) > 0, "cannot find read 2 file at location {} associated with read 1 file {}".format(paired_name,self.filename)
            paired_name = os.path.basename(paired_name)
            self.paired = True
            self.pairedname = paired_name
        elif self.isInterleaved(self.filepath):
            self.paired = True
            self.interleaved = True

            
    '''abstract out the interleaving check for clarity'''
    def isInterleaved(self,fpath):
        pairindex = 2
        if self.fileformat == "fastq":
            pairindex = 4
        if self.head[0].endswith("/1"):
            assert self.head[0].strip("/1") == self.head[pairindex].strip("/2"), "first read in file {} is named in the interleaved format (/1) but does not have a matching read 2 as expected".format(self.filename)
            return True
        else:
            return False
    '''name with no extension to use later'''
    def cleanName(self):
        seqfile_regex =r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)"
        self.cleanname=re.search(seqfile_regex,self.filename).group(1)

    '''count number of reads in the file using IOTools
    this is not run on initialisation'''
    def readCount(self):
        count = IOTools.getNumLines(self.filepath, ignore_comments=False)
        divisor = 2
        if self.fileformat == "fastq":
            divisor = 4
        reads = count/divisor
        self.readcount=reads
        
'''
general class for assembly algorithms
'''
class Assembler():
    def __init__(self,seqdatobj,outdir,paramsobj):
        self.seqdat = seqdatobj
        self.outd = outdir
        self.params = paramsobj
        self.command = ""
        self.incall = ""
        self.outcall = ""
        self.methcall = ""
        self.dstr = os.getcwd()+"/"

    #these functions are overriden in child classes specific to each assembler    
    def inBuild(self):
        pass
    def outBuild(self):
        pass
    def methBuild(self):
        pass

    #basic command construction same for each - command, input files, method specific options, output file
    def construct(self,methodcommand):
        #build sections 
        self.inBuild()
        self.outBuild()
        self.methBuild()
        statement = " ".join([methodcommand,self.incall,self.methcall,self.outcall])
        return statement

'''
megahit assembler class
'''    
class Megahit(Assembler):
    #build the part of the cmd that describes the input seq data
    def inBuild(self):
        instring = ""
        if self.seqdat.paired == True:
            if self.seqdat.interleaved == True:
                instring += "--12 {}".format(self.dstr+self.seqdat.filename)
            else:
                instring += "-1 {} -2 {}".format(self.dstr+self.seqdat.filename,self.dstr+self.seqdat.pairedname)
        else:
            instring += "-r {}".format(self.dstr+self.seqdat.filename)
        self.incall = instring

    #build the command that sets output, megahit does its own suffix so provide inputname and dir
    def outBuild(self):
        outstring = "--out-dir {} --out-prefix {}".format(self.dstr+self.outd+"/"+self.seqdat.cleanname,self.seqdat.cleanname)
        self.outcall = outstring
        
    #build the metahit parameters from ini file
    def methBuild(self):
        m="Megahit_"
        methlist=[]
        methlist.append("--min-count {}".format(self.params[m+"min_count"]))
        if self.params[m+"use_k_list"] == "true":
            methlist.append("--k-list {}".format(self.params[m+"k_list"]))
        else:
            methlist.append("--k-min {} --k-max {} --k-step {}".format(
                self.params[m+"k_min"],self.params[m+"k_max"],self.params[m+"k_step"]))
        if self.params[m+"no_mercy"] == "true":
            methlist.append("--no-mercy")
        methlist.append("--bubble-level {}".format(self.params[m+"bubble_level"]))
        methlist.append("--merge-level {}".format(self.params[m+"merge_level"]))
        methlist.append("--prune-level {}".format(self.params[m+"prune_level"]))
        methlist.append("--prune-depth {}".format(self.params[m+"prune_depth"]))
        methlist.append("--low-local-ratio {}".format(self.params[m+"low_local_ratio"]))
        methlist.append("--max-tip-len {}".format(self.params[m+"max_tip_len"]))
        if self.params[m+"no_local"] == "true":
            methlist.append("--no-local")
        if self.params[m+"kmin_1pass"] == "true":
            methlist.append("--kmin-1pass")
        if self.params[m+"presets"] != "false":
            methlist.append("--presets {}".format(self.params[m+"presets"]))
        methlist.append("--memory {}".format(self.params[m+"memory"]))
        methlist.append("--mem-flag {}".format(self.params[m+"mem_flag"]))
        if self.params[m+"num_cpu_threads"] != "false":
            methlist.append("--num-cpu-threads {}".format(self.params[m+"num_cpu_threads"]))
        methlist.append("--min-contig-len {}".format(self.params[m+"min_contig_len"]))
        if self.params[m+"keep_tmp_files"] != "false":
            methlist.append("--keep-tmp-files")
        self.methcall = " ".join(methlist)
        
    def build(self):
        return self.construct("megahit")
                
'''
Metaspades assembler class, all follow basic format as megahit
'''
class Metaspades(Assembler):
    def inBuild(self):
        instring=""
        assert self.seqdat.paired == True, "metaspades currently only supports paired-end reads."
        if self.seqdat.interleaved == True:
            instring += "--12 {}".format(self.dstr+self.seqdat.filename)
        else:
            instring += "-1 {} -2 {}".format(self.dstr+self.seqdat.filename,self.dstr+self.seqdat.pairedname)
        self.incall = instring

        
    def outBuild(self):
        outstring = "-o {}".format(self.dstr+self.outd+"/"+self.seqdat.cleanname)
        self.outcall = outstring
    
    def methBuild(self):
        m="Metaspades_"
        methlist=[]
        if self.params[m+"only_error_correction"]=="true":
            methlist.append("--only-error-correction")
        if self.params[m+"only_assembler"]=="true":
            methlist.append("--only-assembler")
        if self.params[m+"continue"]=="true":
            methlist.append("--continue")
        if self.params[m+"restart_from"]!="false":
            methlist.append("--restart-from {}".format(self.params[m+"restart_from"]))
        if self.params[m+"disable_gzip_output"]=="true":
            methlist.append("--disable-gzip-output")
        if self.params[m+"disable_rr"]=="true":
            methlist.append("--disable-rr")
        if self.params[m+"dataset"]!="false":
            methlist.append("--dataset {}".format(self.params[m+"dataset"]))
        methlist.append("--threads {}".format(self.params[m+"threads"]))
        if int(self.params[m+"threads"])>1:
            methlist.append("--memory {}".format(int(self.params[m+"memory"])*(int(self.params[m+"threads"])-1)))
        else:
            methlist.append("--memory {}".format(self.params[m+"memory"]))
        methlist.append("-k {}".format(self.params[m+"k"]))
        methlist.append("--phred-offset {}".format(self.params[m+"phred_offset"]))
        self.methcall = " ".join(methlist)

    def build(self):
        return self.construct("metaspades.py")

'''
IDBA-UD assembler class following construction pattern as previously
'''
class Idbaud (Assembler):
    def inBuild(self):
        self.seqdat.cleanname=re.search("(\S+).(interleaved)",self.seqdat.cleanname).group(1)
        instring = "--read {}".format(self.dstr+self.outd+"/"+self.seqdat.cleanname+"/"+self.seqdat.filename)
        self.incall = instring
        
    def outBuild(self):
        outstring = "--out {}".format(self.dstr+self.outd+"/"+self.seqdat.cleanname+"/")
        self.outcall = outstring

    def methBuild(self):
        m="IDBAUD_"
        methlist=[]
        methlist.append("--mink {}".format(self.params[m+"mink"]))
        methlist.append("--maxk {}".format(self.params[m+"maxk"]))
        methlist.append("--step {}".format(self.params[m+"step"]))
        methlist.append("--inner_mink {}".format(self.params[m+"inner_mink"]))
        methlist.append("--inner_step {}".format(self.params[m+"inner_step"]))
        methlist.append("--prefix {}".format(self.params[m+"prefix"]))
        methlist.append("--min_count {}".format(self.params[m+"min_count"]))
        methlist.append("--min_support {}".format(self.params[m+"min_support"]))
        methlist.append("--num_threads {}".format(self.params[m+"clus_threads"]))
        methlist.append("--seed_kmer {}".format(self.params[m+"seed_kmer"]))
        methlist.append("--min_contig {}".format(self.params[m+"min_contig"]))
        methlist.append("--similar {}".format(self.params[m+"similar"]))
        methlist.append("--max_mismatch {}".format(self.params[m+"max_mismatch"]))
        methlist.append("--min_pairs {}".format(self.params[m+"min_pairs"]))
        if self.params[m+"no_bubble"] == "true":
            methlist.append("--no_bubble")
        if self.params[m+"no_local"] == "true":
            methlist.append("--no_local")
        if self.params[m+"no_coverage"] == "true":
            methlist.append("--no_coverage")
        if self.params[m+"no_correct"] == "true":
            methlist.append("--no_correct")
        if self.params[m+"pre_correction"] == "true":
            methlist.append("--pre_correction")
        self.methcall = " ".join(methlist)
            
    def build(self):
        return self.construct("idba_ud")
'''
Function to handle file conversion before IDBA-UD
'''
def IdbaudInterleave(seqdat,indir,outfile):
    infile=indir+"/"+seqdat.filename
    outfile=indir+"/"+outfile
    statementlist=[]
    #if not paired print message
    if seqdat.paired == False:
        print("IDBA-UD cannot assemble {} as it requires paired-end data".format(seqdat.filename))
        return None
    #if already interleaved ensure uncompressed fasta
    elif seqdat.interleaved == True:
        if seqdat.fileformat == "fasta":
            if seqdat.compressed == False:
                #file already correct format so just symlink to new filename
                statementlist.append("ln -s {} {}".format(infile,outfile))
            else:
                #just decompress to new location
                statementlist.append("zcat {} >{}".format(infile,outfile))
        else:
            if seqdat.compressed == False:
                #just convert to fasta using idbas function
                statementlist.append("fq2fa --paired --filter {} {}".format(infile,outfile))
            else:
                #decompress and convert to fasta
                statementlist.append("zcat {} >{} ".format(infile,outfile+".gztemp"))
                statementlist.append("fq2fa --paired --filter {} {}".format(outfile+".gztemp",outfile))
                statementlist.append("rm {}".format(outfile+".gztemp"))
    #if not interleaved do appropriate conversions and interleaving            
    else:
        pairedfile=indir+"/"+seqdat.pairedname
        file1=infile
        file2=pairedfile
        statementlist.append("seqtk mergepe {} {} >{} ".format(file1,file2,outfile+"_temp"))
        if seqdat.fileformat == "fastq":
            statementlist.append("fq2fa --paired --filter {} {}".format(outfile+"_temp",outfile))
            statementlist.append("rm {}".format(outfile+"_temp"))
        else:
            statementlist.append("mv {} {}".format(outfile+"_temp",outfile))
            
    statement = " && ".join(statementlist)
    return statement
            
            
"""
Function to generate summary statistics for the contigs in a fasta file
"""
def SummariseContigs(infile,outfile):
    dstr=os.getcwd()+"/"
    opath=dstr+outfile
    ipath=dstr+infile
    #call BBMap stats command, generates N50, counts, hist etc.
    command = "stats.sh in={} extended=t format=5 shist={} gchist={} overwrite=true >{}".format(
        ipath,opath+".shist",opath+".gchist",opath)
    return command

    
