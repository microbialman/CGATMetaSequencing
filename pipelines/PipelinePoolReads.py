'''
classes and utility functions for pipeline_poolreads.py

'''

import os
import PipelineMetaAssemblyKit
import gzip

'''
function to generate concatenation statement
'''

def poolReads(ftype,infiles,outname):
    statementlist = []
    if ftype.paired == True and ftype.interleaved == False:
        statementlist.append("touch {} && touch {}".format(outname+".1",outname+".2"))
    else:
        statementlist.append("touch {}".format(outname))
    #concatenate the reads as appropriate
    for i in infiles:
        curfile = PipelineMetaAssemblyKit.SequencingData(i)
        if ftype.paired == True and ftype.interleaved == False:
            statementlist.append("zcat -f {} >> {} && zcat -f {} >> {}".format(i,outname+".1",curfile.pairedname,outname+".2"))
        else:
            statementlist.append("zcat -f {} >> {}".format(i,outname))
    #if compressed in, compress the output
    if ftype.compressed == True:
        if ftype.paired == True and ftype.interleaved == False:
            statementlist.append("gzip {} && gzip {}".format(outname+".1",outname+".2"))
        else:
            statementlist.append("gzip {}".format(outname))
    return(" && ".join(statementlist))

'''
Function retunring the pooled file location from the input file list, required to get around ruffus merge limits
'''
def pooledName(infiles,PARAMS):
    ftype = PipelineMetaAssemblyKit.SequencingData(infiles[0])
    pooledname = "pooled.dir/"+PARAMS["General_output_prefix"]+"."+ftype.fileformat
    if ftype.paired == True and ftype.interleaved == False:
        pooledname += ".1"
    if ftype.compressed == True:
        pooledname += ".gz"
    return(PipelineMetaAssemblyKit.SequencingData(pooledname))


'''
function to rename the headers within a pooled reads file to ensure unqiue names
'''
def renameNums(pfile,path):
    #open files for editing
    if pfile.compressed == True:
        fileopen = gzip.open(path+pfile.filename,'r')
    else:
        fileopen = open(path+pfile.filename,'r')
    if pfile.paired == True and pfile.interleaved == False:
        if pfile.compressed == True:
            pairedopen = gzip.open(path+pfile.pairedname,'r')
        else:
            pairedopen = open(path+pfile.pairedname,'r')

    #open tempfiles to write to
    fileout = open(path+"tempout",'w')
    if pfile.paired == True and pfile.interleaved == False:
        pairedout = open(path+"tempout2",'w')

    #edit the files appropriately
    if pfile.fileformat == "fastq":
        ch = "@"
    else:
        ch = ">"
    if pfile.compressed == True:
        isbyte = True
    else:
        isbyte = False
        
    numFun(fileopen,fileout,ch,isbyte,pfile.interleaved)
    if pfile.paired == True and pfile.interleaved == False:
        numFun(pairedopen,pairedout,ch,isbyte,pfile.interleaved)

    #close the files
    fileopen.close()
    fileout.close()
    if pfile.paired == True and pfile.interleaved == False:
        pairedopen.close()
        pairedout.close()
    
    #generate commands to remove originals and rename temp files
    statementlist=[]
    statementlist.append("rm {} && mv {} {}".format(path+pfile.filename,path+"tempout",path+pfile.filename.rstrip(".gz")))
    if pfile.compressed == True:
        statementlist.append("gzip {}".format(path+pfile.filename.rstrip(".gz")))
    if pfile.paired == True and pfile.interleaved == False:
        statementlist.append("rm {} && mv {} {}".format(path+pfile.pairedname,path+"tempout2",path+pfile.pairedname.rstrip(".gz")))
        if pfile.compressed == True:
            statementlist.append("gzip {}".format(path+pfile.pairedname.rstrip(".gz")))
    return(" && ".join(statementlist))
    


#function to add numbers to file lines based on rules
def numFun(infile,outfile,char,isbyte,interleaved):
    counter = 1
    linecount = 0
    interpair = False

    def editline(line,counter,char):
        edited = line[1:]
        edited = char+str(counter)+edited
        return(edited)

    def interleaveCheck(interleaved,counter,interpair):
        if interleaved == False:
            counter+=1
        else:
            if interpair == True:
                counter+=1
                interpair=False
            else:
                interpair=True
        return(counter,interpair)
    
    for i in infile:
        if isbyte:
            i = i.decode()
        #different actions for fastq and fasta
        if char == "@":
            if linecount !=0 :
                outfile.write(i)
                linecount+=1
            if linecount == 0:
                outfile.write(editline(i,counter,char))
                linecount+=1
            if linecount==4:
                linecount=0
                counter, interpair = interleaveCheck(interleaved,counter,interpair)
                
        #for fasta
        else:
            if i[0] != ">":
                outfile.write(i)
            else:
                outfile.write(editline(i,counter,char))
                counter, interpair = interleaveCheck(interleaved,counter,interpair)
            
            

#function to derive the dereplication call to bbmap dedupe.sh
def derepReads(infile,path,params):
    statementlist = []
    dupelist = ["dedupe.sh"]
    if infile.paired == False:
        dupelist.append("in={}".format(path+infile.filename))
    elif infile.interleaved == True:
        dupelist.append("in={} interleaved=true".format(path+infile.filename))
    else:
        dupelist.append("in1={} in2={}".format(path+infile.filename,path+infile.pairedname))
    if infile.paired == True and infile.interleaved == False:
        dupelist.append("out={} ac=f".format(path+"tempdereplicated-"+infile.filename))
        statementlist.append("reformat.sh in={} out1={} out2={}".format(
            path+"tempdereplicated-"+infile.filename,
            path+"dereplicated-"+infile.filename,
            path+"dereplicated-"+infile.pairedname))
        statementlist.append("rm {}".format(path+"tempdereplicated-"+infile.filename))
    else:
        dupelist.append("out={} ac=f".format(path+"dereplicated-"+infile.filename))
    statementlist.insert(0," ".join(dupelist))
    return(" && ".join(statementlist))
