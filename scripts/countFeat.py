from argparse import ArgumentParser
import numpy as np
import os

def logCount(counter):
    c = counter+1
    if c%1000000 == 0:
        logfile.write(str(counter)+"\n")
        logfile.flush()
    return(c)

#get the input and output files
parser = ArgumentParser()
parser.add_argument("--orfinput", dest="infile", help="input file with orf counts")
parser.add_argument("--gtf", dest="gtf", help="gtf file with mapping from orf to feature")
parser.add_argument("--feature",dest="feat", help="feature to enumerate from the gtf and orf counts")
parser.add_argument("--output", dest="outfile", help="output file with one annotation per row, each annotation carries counts from original feature")
args = parser.parse_args()

#open the files
infile = open(args.infile,"rU")
gtffile = open(args.gtf,"rU")
logfile = open(os.path.dirname(args.outfile)+"/gtflog_"+os.path.basename(args.outfile),'w')
outfile = open(args.outfile,"w")
f = args.feat

orfs={}
feats={}

logfile.write("Parsing GTF\n\n")
logfile.flush()

counter=0

#parse the GTF to get the features
for i in gtffile:
    counter=logCount(counter)
    row = i.strip("\n").split("\t")
    start = int(row[3])
    end = int(row[4])
    length = end-start
    orfdat = row[-1].split(";")
    fval = ''
    for j in range(len(orfdat)-1):
        feature = orfdat[j].split('"')[0].replace(" ","")
        if feature == "gene_id":
            gid = orfdat[j].split('"')[1]
            gid = gid.replace('"','')
        if feature == f:
            fval = orfdat[j].split('"')[1]
            fval = fval.replace('"','')
    if fval != '':
        if gid not in orfs:
            orfs[gid]=[]
        for k in fval.split(","):
            orfs[gid].append(k)
            if k not in feats:
                feats[k]={"length":0}
            feats[k]["length"]+=length

logfile.write(str(counter))
logfile.write("\n\nParsing counts\n\n")
counter=0
            
#parse the counts file
header = False
for i in infile:
    counter=logCount(counter)
    if i[0] == '#':
        pass
    elif header == False:
        header = i.strip("\n").split("\t")
        samples = header[1:]
        #initiate an empty array of correct length in each feature
        for k in feats:
            feats[k]["counts"]=np.zeros(len(samples))
    #for each orf add counts to its matching features
    else:
        row = i.strip("\n").split("\t")
        orfid = row[0].replace('"','')
        if orfid in orfs:
            sampcounts = np.array(row[1:],dtype=float)
            for j in orfs[orfid]:
                feats[j]["counts"]+=sampcounts


#write data to output
outfile.write("{}\t{}\n".format(f,"\t".join(samples)))
for i in feats:
    strlist = np.char.mod('%f', feats[i]["counts"])
    outfile.write("{}\t{}\n".format(i,"\t".join(strlist)))
    outfile.flush()

logfile.write(str(counter))
    
outfile.close()
logfile.close()
                
