from argparse import ArgumentParser
import re

#get the input and output files
parser = ArgumentParser()
parser.add_argument("--orf_counts", dest="infile", help="input file with orf counts")
parser.add_argument("--gtf", dest="gtf", help="gtf file with mapping from orf to features")
parser.add_argument("--features",dest="feat", help="comma-seperated list of features to enumerate from the gtf and orf counts")
parser.add_argument("--outdir", dest="outdir", help="starting output directory, each feature will be given a child-directory with the sample-level counts")
parser.add_argument("--logfile", dest="logfile", help="destination for the log file")
args = parser.parse_args()

#open the files
infile = open(args.infile,"rU")
samplename = re.search("orf_counts.dir/(\S+).tsv",args.infile).group(1)
gtffile = open(args.gtf,"rU")
f = args.feat.split(",")

#parse the counts
orfcounts={}
header=0
for i in infile:
    if header<2:
       header+=1
    else:
        row=i.strip("\n").split()
        orfcounts[row[0]]=float(row[-1])

def parseAnnot(geneid,annot):
    if annot == "":
        return(None)
    else:
        annot=annot.split('"')
        feat=annot[0].strip(" ")
        val=annot[1].strip('"')
        if feat in f:
            if val == "":
                return(None)
            else:
                vals=val.split(',')
                for y in vals:
                    if y not in featurecounts[feat]:
                        featurecounts[feat][y]=0.0
                    featurecounts[feat][y]+=orfcounts[geneid]
                   
    
featurecounts={}
for i in f:
    featurecounts[i]={}
#parse the features
for i in gtffile:
    row=i.strip("\n").split("\t")
    annotations=row[-1].split(";")
    geneid=annotations[0].split(" ")[1].strip('"')
    for x in annotations[1:]:
        parseAnnot(geneid,x)

#write out the files
for i in featurecounts:
    outfile=open("{}{}/{}.tsv".format(args.outdir,i,samplename),"w")
    outfile.write("{}\tCount\n".format(i))
    for k in featurecounts[i]:
        outfile.write("{}\t{}\n".format(k,featurecounts[i][k]))
    outfile.close()

#write the log file
log=open(args.logfile,"w")
log.write("Counts aggregated for features {} for sample {}.".format(",".join(f),args.infile))
log.close()
