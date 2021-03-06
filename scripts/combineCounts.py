from argparse import ArgumentParser
import re, glob

#get the input and output files
parser = ArgumentParser()
parser.add_argument("--feature", dest="feat", help="feature to combine counts for")
parser.add_argument("--countdir", dest="cdir", help="head directory contianing sub-directories with feature counts")
parser.add_argument("--outfile",dest="outfile", help="output destination for combined table")
parser.add_argument("--tpm",dest="tpm", help="true will  combine '.tpm.tsv' files, false '.tsv' files")
args = parser.parse_args()

#get all the files in the count directory
globs=glob.glob("{}/{}/*.tsv".format(args.cdir,args.feat))
#screen tpm as necessary
sampfiles=[]
for i in globs:
    if re.search("\.tpm\.tsv",i):
        if args.tpm == "true":
            sampfiles.append(i)
    else:
        if args.tpm != "true":
            sampfiles.append(i)

#set suffix based on tpm
suffix="\.tsv"
if args.tpm == "true":
    suffix="\.tpm\.tsv"

#dictionary to hold the counts
featcounts={}
samplenames=[]
featnames=[]

#parse the files
for i in range(len(sampfiles)):
    fileopen=open(sampfiles[i],'rU').readlines()
    sample=re.search("{}/(\S+){}".format(args.feat,suffix),sampfiles[i]).group(1)
    samplenames.append(sample)
    for j in fileopen[1:]:
        row=j.strip("\n").split("\t")
        if row[0] not in featcounts:
            featnames.append(row[0])
            featcounts[row[0]]=[0.0]*len(sampfiles)
        featcounts[row[0]][i]=float(row[1])

#write the output
outfile=open(args.outfile,"w")
outfile.write("{}\t{}\n".format(args.feat,"\t".join(samplenames)))
for i in featnames:
    outfile.write("{}\t{}\n".format(i,"\t".join([str(x) for x in featcounts[i]])))
outfile.close()
