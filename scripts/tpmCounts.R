library(data.table)

#load the counts table
counts=fread(commandArgs(trailingOnly=T)[1],sep="\t",header=T)
counts=as.data.frame(counts)

#write the untransformed table out, removing excess cols
write.table(counts[,c(1,7:ncol(counts))],commandArgs(trailingOnly=T)[2],sep="\t",row.names = F)
#tpm transform and write out
counts[,7:ncol(counts)]=apply(counts[,7:ncol(counts)],2,function(x) x/counts[,6])
csums=colSums(counts[,7:ncol(counts)])
counts[,7:ncol(counts)]=sweep(counts[,7:ncol(counts)],2,csums,`/`)
counts[,7:ncol(counts)]=counts[,7:ncol(counts)]*1000000
write.table(counts[,c(1,7:ncol(counts))],commandArgs(trailingOnly=T)[3],sep="\t",row.names = F)
