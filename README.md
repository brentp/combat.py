ComBat is an R package for removing batch effects from data.
This is a python version that matches the output from the ComBat function
in SVA (http://www.bioconductor.org/packages/release/bioc/html/sva.html).

To test, run this R code (requires sva and bladderbatch from bioconductor):

```R
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


```

Then, from the same directory, run

```Shell

python combat.py

```

and you should see that it runs much faster (> 10X) and outputs identical results
to 4 decimal places.


