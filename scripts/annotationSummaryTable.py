#this scrpt generates a summary table of the annotations that can then be used in the report script
from argparse import ArgumentParser
import re, glob
import gzip

#get the orfs, function and taxonomy annotations from commadn line
parser = ArgumentParser()
parser.add_argument("--gtf-dir", dest="gtf", help="directory of the GTFs containing the combined functional and taxonomic annotations for each sample")
parser.add_argument("--annot-output", dest="outfile", help="Annotation summary output filename")
parser.add_argument("--orf-output", dest="orfout", help="ORF summary output filename")
args = parser.parse_args()

#get the list of gtfs (full not short) for each sample
gtfs = glob.glob(args.gtf+"/*.orf_annotations.gtf.gz")

featdic={}

taxlevs=["kingdom","phylum","class","order","family","genus","species"]
nonlistfun=["Predicted_protein_name"]
listfun=["Gene_Ontology_terms","EC_number","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","BRITE","BiGG_Reaction","eggNOG_OGs","COG_Functional_Category"]

#open the orf file
orffile=open(args.orfout,"w")
orffile.write("ORF\tContig\tSize\tSample\n")

#function to parse features, some are listed (i.e. one annotation conatins many feautes, eg. GO terms)
def featAdd(sample,feat,listed):
    def addDic(key,val):
        if key not in featdic:
            featdic[key]={"type":val,"samples":[0]*len(gtfs)}
        if featdic[key]["samples"][i]==0:
            featdic[key]["samples"][i]=1
    if listed == False:
        if feat[1] != "":
            addDic(feat[1],feat[0])
    else:
        for x in feat[1].split(","):
            if x != "":
                addDic(x,feat[0])

#go through gtf files and record occurences of each feature
for i in range(len(gtfs)):
    #open annotation
    gtf=gzip.open(gtfs[i],"rt")
    #parse
    for j in gtf:
        row=j.strip("\n").split("\t")
        #parse the orf
        contig=row[0]
        length=str(int(row[4])-int(row[3]))
        orf=row[-1].split(";")[0].split('"')[1]
        orffile.write("{}\t{}\t{}\t{}\n".format(orf,contig,length,gtfs[i]))
        #parse the annotations
        annotations=[x.split('"') for x in row[-1].split(";")]
        for k in annotations[1:]:
            k[0]=k[0].strip(" ")
            if k[0] in nonlistfun or k[0] in taxlevs:
                featAdd(i,k,False)
            elif k[0] in listfun:
                featAdd(i,k,True)
orffile.close()

#write the output file
outfile=open(args.outfile,"w")
outfile.write("feat_name\tfeat_type\t"+"\t".join(gtfs)+"\n")
for i in featdic:
    outfile.write("{}\t{}\t{}\n".format(i,featdic[i]["type"],"\t".join([str(x) for x in featdic[i]["samples"]])))
outfile.close()
