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

Performance
===========

On an identical dataset, of 30K rows * 190 samples, this python version finishes in 10.008s
as measured by unix `time`.
The R version takes 4m0.681s with output identical to 3 decimal places. This is a speed-up
of about *24x.*



Function
========

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

Read
====

    Johnson WE, Rabinovic A, Li C (2007). Adjusting batch effects in microarray
    expression data using Empirical Bayes methods. Biostatistics 8:118-127.  

    Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker, Andrew E. Jaffe
    and John D. Storey (). sva: Surrogate Variable Analysis. R package
    version 3.4.0.

