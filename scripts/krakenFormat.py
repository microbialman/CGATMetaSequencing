from argparse import ArgumentParser

#get the orfs, function and taxonomy annotations from command line
parser = ArgumentParser()
parser.add_argument("--orfs", dest="orfs", help="fasta containing ORF sequences from prodigal")
parser.add_argument("--contig-taxonomy", dest="contigs", help="Contig taxonomic assignments from Kraken")
parser.add_argument("--orf-taxonomy-output", dest="output", help="Destination for output of converted taxonomic assignments")
args = parser.parse_args()

#load in the ORFs
orfdic={}
orffile=open(args.orfs,"rU")
for i in orffile:
    if i[0] == ">":
        orfname=i.split(" #")[0].strip(">")
        contigname="_".join(orfname.split("_")[:-1])
        orfdic[orfname]=contigname
        
#load in the contig taxonomies
contigdic={}
contigfile=open(args.contigs,"rU")
for i in contigfile:
    row=i.strip("\n").split()
    contigdic[row[0]]=row[1].split("|")


outfile=open(args.output,"w")
#write the output
for i in orfdic:
    row=[i," "]
    if orfdic[i] in contigdic:
        for j in contigdic[orfdic[i]]:
            row.append(j)
            row.append(" 100")
        outfile.write(";".join(row)+";\n")
        
    
        
    
