#translates NCBI taxids in default kraken2 output to names in 'mpa' style using Taxonkit
from argparse import ArgumentParser
import subprocess

#get the default kraken2 output and the file to save translated names file to
parser = ArgumentParser()
parser.add_argument("--krakenout",dest="infile", help="Kraken output file to be translated")
parser.add_argument("--translatedout", dest="outfile", help="Output file name")
parser.add_argument("--taxdatadir", dest="tddir", help="Directory containing names.dmp and nodes.dmp from NCBI taoxnomy (for TaxonKit, uses taxonkit default if not set).", )
args = parser.parse_args()

#open the files
infile = open(args.infile, 'rU')
outfile = open(args.outfile, 'w')
        
#read name to taxid mapping
readnames=[]
taxids=[]

#unique taxids to find the taxonomic names for
uniqueids=[]

#go through each kraken result and get taxids
for i in infile:
    row = i.split("\t")
    readname = row[1]
    taxid = row[2]
    readnames.append(readname)
    taxids.append(taxid)
    if taxid not in uniqueids:
        uniqueids.append(taxid)

#dictionary to store found ids
taxdic={}
        
#get the full ineage names for unique taxids using Taxonkit
if args.tddir:
    taxonkit = subprocess.check_output("echo '{}' | taxonkit lineage --data-dir {} | taxonkit reformat --data-dir {}".format("\n".join(uniqueids)), shell=True)
else:
    taxonkit = subprocess.check_output("echo '{}' | taxonkit lineage | taxonkit reformat".format("\n".join(uniqueids)), shell=True)
taxonkit=taxonkit.decode().split("\n")

#taxonomic levels
levs=["d__","p__","c__","o__","f__","g__","s__"]

#function to generate the mpa name from a given taxid by calling Taxonkit
def formName(name):
    names = name.split(";")
    formnames = []
    for i in range(len(names)):
        if names[i] != "":
            formnames.append(levs[i]+names[i].replace(" ","_"))
    return("|".join(formnames))

#reformat to mpa style
for i in taxonkit:
    if i != "":
        row=i.split("\t")
        tid=row[0]
        mpaname=formName(row[2])
        if mpaname != "":
            taxdic[tid]=mpaname

#write the translated read to taxname file
for i in range(len(readnames)):
    t=taxids[i]
    if t in taxdic:
        outfile.write("{}\t{}\n".format(readnames[i],taxdic[t]))
