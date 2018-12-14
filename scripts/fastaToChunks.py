from argparse import ArgumentParser
import re

#get the orfs, function and taxonomy annotations from commadn line
parser = ArgumentParser()
parser.add_argument("--input", dest="inputfile", help="input fasta")
parser.add_argument("--output_prefix", dest="out", help="output prefix")
parser.add_argument("--chunk_size", dest="chunk", help="chunk size in no. of reads per chunk")
args = parser.parse_args()

#open the input file
fasta = open(args.inputfile,'rU')
outpre = args.out
chunk = int(args.chunk)

#initiate counts
current_chunk=1
read_count=0

#open first chunk file
chunkfile=open(outpre+".{}.chunk".format(current_chunk),"w")

#write chunks
for i in fasta:
    if i[0]==">":
        if read_count == chunk:
            chunkfile.close()
            current_chunk+=1
            chunkfile=open(outpre+".{}.chunk".format(current_chunk),"w")
            read_count = 0
        else:
            read_count+=1
    chunkfile.write(i)

chunkfile.close()
