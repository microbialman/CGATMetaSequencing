#takes a feature counts output table and generates a version with an additional TPM column
from argparse import ArgumentParser

#get the input and output files
parser = ArgumentParser()
parser.add_argument("--orf_counts", dest="infile", help="input file with orf counts")
parser.add_argument("--output", dest="outfile", help="output file that will have TPM counts")
args = parser.parse_args()

#open the files
infile=open(args.infile,"rU")
outfile=open(args.outfile,"w")

#parse the counts
orfrpk=[]
rows=[]
header=0
for i in infile:
    rows.append(i)
    if header<2:
        header+=1
    else:
        row=i.strip("\n").split()
        count=int(row[-1])
        length=int(row[-2])/1000
        orfrpk.append(count/length)

#convert to TPM and write out
scaling=sum(orfrpk)/1000000
header=0
rowc=0
for i in rows:
    if header<1:
        outfile.write(i)
        header+=1
    elif header<2:
        outfile.write(i.strip("\n")+"\ttpm_counts\n")
        header+=1
    else:
        tpmval=str(orfrpk[rowc]/scaling)
        outfile.write(i.strip("\n")+"\t{}\n".format(tpmval))
        rowc+=1
outfile.close()
        
        
        
        
