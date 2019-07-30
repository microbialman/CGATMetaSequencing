from argparse import ArgumentParser
import re
import gzip

#get the orfs, function and taxonomy annotations from commadn line
parser = ArgumentParser()
parser.add_argument("--input", dest="inputfile", help="input fasta")
parser.add_argument("--output_prefix", dest="out", help="output prefix")
parser.add_argument("--chunk_size", dest="chunk", help="chunk size in no. of reads per chunk")
args = parser.parse_args()

#open the input file
fasta = gzip.open(args.inputfile)
outpre = args.out
chunk = int(args.chunk)

#initiate counts
current_chunk=1
read_count=0

#open first chunk file
chunkfile=open(outpre+".{}.chunk".format(current_chunk),"w")
chunklog=open(outpre+".{}.chunk.log".format(current_chunk),"w")
chunklog.write("Chunk no. {}.".format(current_chunk))

#write chunks
for i in fasta:
    i=i.decode()
    print(i)
    if i[0]==">":
        if read_count == chunk:
            chunkfile.close()
            chunklog.close()
            current_chunk+=1
            chunkfile=open(outpre+".{}.chunk".format(current_chunk),"w")
            chunklog.write("Chunk no. {}.".format(current_chunk))
            read_count = 0
        else:
            read_count+=1
    chunkfile.write(i)

chunkfile.close()
chunklog.close()
