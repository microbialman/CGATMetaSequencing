from argparse import ArgumentParser
import re, gzip

#get the input and output files
parser = ArgumentParser()
parser.add_argument("--orf_counts", dest="infile", help="input file with orf counts")
parser.add_argument("--gtf", dest="gtf", help="gtf file with mapping from orf to features")
parser.add_argument("--feature_pairs",dest="feat", help="comma-seperated list of features to enumerate from the gtf and orf counts")
parser.add_argument("--multimethod",dest="multimethod",help="How to handle counts when a gene maps to multiple annotations of a feature? whole = add the gene's whole count to all of its annotations, split = divide the gene count evenly across its annotations")
parser.add_argument("--outdir", dest="outdir", help="starting output directory, each feature will be given a child-directory with the sample-level counts")
parser.add_argument("--logfile", dest="logfile", help="destination for the log file")
args = parser.parse_args()

#open the files
infile = gzip.open(args.infile,"rt")
samplename = re.search("orf_counts.dir/(\S+).tsv",args.infile).group(1)
gtffile = gzip.open(args.gtf,"rt")
fp = args.feat.split(",")
fp = [(x.split("_BY_")[0],x.split("_BY_")[1]) for x in fp]
multimethod = args.multimethod

#parse the counts
orfcounts={}
header=0
for i in infile:
    if header<2:
       header+=1
    else:
        row=i.strip("\n").split()
        orfcounts[row[0]]=float(row[-1])

def featCount(fname,j,c):
    if fname not in featurecounts[j]:
        featurecounts[j][fname]=0.0
    featurecounts[j][fname]+=c
    
        
        
#function to parse the feature names
def parsePairedAnnot(geneid,annot):
    if annot == "":
        return(None)
    else:
        subdic={}
        for i in annot:
            if i!= "":
                ann=i.split('"')
                feat=ann[0].strip(" ")
                val=ann[1].strip('"').split(",")
                subdic[feat]=val
        for j in fp:
            if j[0] in subdic and j[1] in subdic:
                obcount=orfcounts[geneid]
                if multimethod == "split":
                    obcount=obcount/(len(subdic[j[0]])*len(subdic[j[1]]))
                for k in subdic[j[0]]:
                    for l in subdic[j[1]]:
                        featname="{}_BY_{}".format(k,l)
                        featCount(featname,j,obcount)
            elif j[0] in subdic:
                obcount=orfcounts[geneid]
                if multimethod == "split":
                    obcount=obcount/len(subdic[j[0]])
                for k in subdic[j[0]]:
                    featname="{}_BY_UNANNOTATED".format(k)
                    featCount(featname,j,obcount)
            elif j[1] in subdic:
                obcount=orfcounts[geneid]
                if multimethod == "split":
                    obcount=obcount/len(subdic[j[1]])
                for l in subdic[j[1]]:
                    featname="UNANNOTATED_BY_{}".format(l)
                    featCount(featname,j,obcount)
            else:
                featname="UNANNOTATED_BY_UNANNOTATED"
                featCount(featname,j,orfcounts[geneid])
            
featurecounts={}
for i in fp:
    featurecounts[i]={}
    
#parse the features
for i in gtffile:
    row=i.strip("\n").split("\t")
    annotations=row[-1].split(";")
    geneid=annotations[0].split(" ")[1].strip('"')
    parsePairedAnnot(geneid,annotations)

   
#write out the files
for i in featurecounts:
    pairname="{}_BY_{}".format(i[0],i[1])
    outfile=open("{}{}/{}.tsv".format(args.outdir,pairname,samplename),"w")
    outfile.write("{}\tCount\n".format(pairname))
    for k in featurecounts[i]:
        outfile.write("{}\t{}\n".format(k,featurecounts[i][k]))
    outfile.close()

#write the log file
log=open(args.logfile,"w")
log.write("Counts aggregated for pairs {} for sample {}.".format(",".join(["_BY_".join(x) for x in fp]),args.infile))
log.close()
