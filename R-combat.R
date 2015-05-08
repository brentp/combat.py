#!/usr/bin/env R

source("http://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("bladderbatch")
library("sva")
options(stringsAsFactors=FALSE)

library(bladderbatch)
data(bladderdata)

pheno = pData(bladderEset)
# add fake age variable for numeric
pheno$age = c(1:7, rep(1:10, 5))
write.table(data.frame(cel=rownames(pheno), pheno), row.names=F, quote=F, sep="\t", file="bladder-pheno.txt")

edata = exprs(bladderEset)
write.table(edata, row.names=T, quote=F, sep="\t", file="bladder-expr.txt")
# use dataframe instead of matrix
mod = model.matrix(~as.factor(cancer) + age, data=pheno)
t = Sys.time()
cdata = ComBat(dat=edata, batch=as.factor(pheno$batch), mod=mod, numCov=match("age", colnames(mod)))
print(Sys.time() - t)
print(cdata[1:5, 1:5])

write.table(cdata, "r-batch.txt", sep="\t", quote=F)
