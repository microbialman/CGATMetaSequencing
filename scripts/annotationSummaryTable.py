from argparse import ArgumentParser
import re

#get the orfs, function and taxonomy annotations from commadn line
parser = ArgumentParser()
parser.add_argument("--gtf", dest="gtf", help="GTF containing the function and taxonomic data")
parser.add_argument("--output", dest="outfile", help="Output filename")
args = parser.parse_args()

#open the files
gtf = open(args.gtf,'rU')
outfile = open(args.outfile,'w')
outfile.write("ORF\tContig\tSize\tTaxa\tNo_GO_Terms\tNo_KEGG_KOs\tNo_OGs\tNo_COG_Cats\n")

#go through gtf files and tabulate ORFs
for i in gtf:
    row = i.split("\t")
    contig = row[0]
    size = str(int(row[4])-int(row[3]))
    vals = row[8].strip("\n").split(";")
    fundic = {}
    for j in vals:
        feat = j.split('=')
        fundic[feat[0]]=feat[1]
    orf = fundic["gene_id"]
    if "phylum" in fundic:
        taxa="1"
    else:
        taxa="NA"
    if "GO_terms" in fundic:
        GOcount = str(fundic["GO_terms"].count("GO"))
        KEGGcount = str(fundic["KEGG_KO"].count("K"))
        EGGcount = str(fundic["Matching_OGs"].count("@"))
        if fundic["COG_functional_categories"] == '""':
            COGcount=str(0)
        else:
            COGcount=str(1+fundic["COG_functional_categories"].count(","))
    else:
        GOcount = "NA"
        KEGGcount = "NA"
        EGGcount = "NA"
        COGcount = "NA"
    vals = [orf, contig, size, taxa, GOcount, KEGGcount, EGGcount,COGcount]
    outfile.write("\t".join(vals).replace('"','')+"\n")
outfile.close()
        
        
