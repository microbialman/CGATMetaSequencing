import argparse
import glob
import gzip

parser=argparse.ArgumentParser()
parser.add_argument("--indir",help="Input file made during the annotation.")
parser.add_argument("--outerind",help="Index for the higher level feature.")
parser.add_argument("--innerind",help="Index for the inner level feature.")
parser.add_argument("--minsetsize",help="Minimum set size to retain in the output.")
parser.add_argument("--outfile",help="Output file.")

args=parser.parse_args()
terms={}

files=glob.glob(args.indir+"/*.annotations.gz")
fcount=0
for i in files:
    fcount+=1
    print("File {}/{}".format(fcount,len(files)))
    with gzip.open(i,"rt") as f:
        for j in f:
            row=j.strip("\n").split("\t")
            inner=row[int(args.innerind)]
            outer=row[int(args.outerind)]
            if inner=="" or outer=="":
                pass
            else:
                inner=inner.split(",")
                outer=outer.split(",")
                for k in outer:
                    if k not in terms:
                        terms[k]=[]
                    terms[k]+=inner
                    
outfile=open(args.outfile,"w")
for i in terms:
    terms[i]=set(terms[i])
    if len(terms[i]) < int(args.minsetsize):
        pass
    else:
        outfile.write("{}\tNA\t{}\n".format(i,"\t".join(terms[i])))
outfile.close()



