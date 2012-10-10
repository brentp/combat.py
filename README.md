ComBat is an R package for removing batch effects from data.
This is a python version that matches the output from the ComBat function
in SVA (http://www.bioconductor.org/packages/release/bioc/html/sva.html).
This code is completely copied from the ComBat function in that package.

To test, run this R code (requires sva and bladderbatch from bioconductor):

```Shell

Rscript R-combat.R

```

Then, from the same directory, run

```Shell

python combat.py

```

and you should see that it runs much faster (> 10X on my machine) and outputs identical results
to 4 decimal places.

The python version is usable as a module, the function has the signature:

```Python

   combat(dat, batch, mod, numCovs=None)

```

which is the same as the R function except the non-parametric version is not supported.

 + dat is the expression/methylation data.
 + batch is a list containing the batch variable
 + mod is the model matrix (can use patsy for this from python)
 + numCovs is a list like ["age", "height"], that gives the column name or number
   of numeric variables in batch (otherwise they will be converted to factors).

