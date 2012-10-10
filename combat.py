import pandas as pa
import patsy
import sys
import numpy.linalg as la
import numpy as np



def adjust_nums(numCovs, drop_idxs):
    # if we dropped some values, have to adjust those with a larger index.
    if numCovs is None: return drop_idxs
    return [nc - sum(nc < di for di in drop_idxs) for nc in numCovs]

def design_mat(mod, numCovs):
    design = patsy.dmatrix("~ 0 + C(batch)", mod, return_type="dataframe")
    mod = mod.drop(["batch"], axis=1)
    print >>sys.stderr, "found %i batches" % design.shape[1]
    other_cols = [c for i, c in enumerate(mod.columns) if not i in numCovs]
    factor_matrix = mod[other_cols]
    design = pa.concat((design, factor_matrix), axis=1)
    if numCovs is not None:
        print >>sys.stderr, "found %i numerical covariates..." % len(numCovs)
        for i, nC in enumerate(numCovs):
            cname = mod.columns[nC]
            print >>sys.stderr, "\t", cname
            design[cname] = mod[mod.columns[nC]]
    print >>sys.stderr, "found %i categorical variables:" % len(other_cols)
    print >>sys.stderr, "\t" + ", ".join(other_cols)
    return design

def combat(dat, batch, mod, numCovs=None):
    if not isinstance(numCovs, (list, tuple)):
        numCovs = [numCovs]
    
    mod["batch"] = batch

    batch_info = [v for k, v in mod.groupby("batch").groups.items()]
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info])
    n_array = float(sum(n_batches))

    # drop intercept
    drop_cols = [cname for cname, inter in  ((mod == 1).all()).iterkv() if inter == True]
    drop_idxs = [list(mod.columns).index(cdrop) for cdrop in drop_cols]
    mod = mod[[c for c in mod.columns if not c in drop_cols]]
    numCovs = [list(mod.columns).index(c) if isinstance(c, basestring) else c for c in numCovs]

    design = design_mat(mod, numCovs)

    B_hat = np.dot(np.dot(la.inv(np.dot(design.T, design)), design.T), dat.T)
    grand_mean = np.dot((n_batches / n_array).T, B_hat[:n_batch,:])

    var_pooled = np.dot(((dat - np.dot(design, B_hat).T)**2), np.ones((n_array, 1)) / n_array)

    stand_mean = np.dot(grand_mean.T.reshape((len(grand_mean), 1)), np.ones((1, n_array)))
    tmp = np.array(design.copy())
    tmp[:,:n_batch] = 0
    stand_mean  += np.dot(tmp, B_hat).T

    s_data = ((dat - stand_mean) / np.dot(np.sqrt(var_pooled), np.ones((1, n_array))))

    batch_design = design[design.columns[:n_batch]]
    gamma_hat = np.dot(np.dot(la.inv(np.dot(batch_design.T, batch_design)), batch_design.T), s_data.T)

    delta_hat = []

    for i, batch_idxs in enumerate(batch_info):
        #batches = [list(mod.columns).index(b) for b in batches] 
        delta_hat.append(s_data[batch_idxs].var(axis=1))

    gamma_bar = gamma_hat.mean(axis=1) 
    t2 = gamma_hat.var(axis=1)
   

    a_prior = map(aprior, delta_hat)
    b_prior = map(bprior, delta_hat)

    print("finding parametric adjustments")
    gamma_star, delta_star = [], []
    for i, batch_idxs in enumerate(batch_info):
        #print '18 20 22 28 29 31 32 33 35 40 46'
        #print batch_info[batch_id]

        temp = it_sol(s_data[batch_idxs], gamma_hat[i],
                     delta_hat[i], gamma_bar[i], t2[i], a_prior[i], b_prior[i])
        # NOTE: OK to here!
        gamma_star.append(temp[0])
        delta_star.append(temp[1])

    print("adjusting data")
    bayesdata = s_data
    gamma_star = np.array(gamma_star)
    delta_star = np.array(delta_star)


    for j, batch_idxs in enumerate(batch_info):

        dsq = np.sqrt(delta_star[j,:])
        dsq = dsq.reshape((len(dsq), 1))
        #print dsq
	denom =  np.dot(dsq, np.ones((1, n_batches[j])))
        numer = np.array(bayesdata[batch_idxs] - np.dot(batch_design.ix[batch_idxs], gamma_star).T)

    	#print batch_design.ix[batch_idxs]
        bayesdata[batch_idxs] = numer / denom
   
    vpsq = np.sqrt(var_pooled).reshape((len(var_pooled), 1))
    bayesdata = bayesdata * np.dot(vpsq, np.ones((1, n_array))) + stand_mean
 
    return bayesdata

def it_sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=0.0001):
    n = (1 - np.isnan(sdat)).sum(axis=1)
    g_old = g_hat.copy()
    d_old = d_hat.copy()

    change = 1
    count = 0
    while change > conv:
        #print g_hat.shape, g_bar.shape, t2.shape
        g_new = postmean(g_hat, g_bar, n, d_old, t2)
        sum2 = ((sdat - np.dot(g_new.reshape((g_new.shape[0], 1)), np.ones((1, sdat.shape[1])))) ** 2).sum(axis=1)
        d_new = postvar(sum2, n, a, b)
       
        change = max((abs(g_new - g_old) / g_old).max(), (abs(d_new - d_old) / d_old).max())
        g_old = g_new
        d_old = d_new
        count = count + 1
    adjust = (g_new, d_new)
    return adjust 

    

def aprior(gamma_hat):
    m = gamma_hat.mean()
    s2 = gamma_hat.var()
    return (2 * s2 +m**2) / s2

def bprior(gamma_hat):
    m = gamma_hat.mean()
    s2 = gamma_hat.var()
    return (m*s2+m**3)/s2

def postmean(g_hat, g_bar, n, d_star, t2):
    return (t2*n*g_hat+d_star * g_bar) / (t2*n+d_star)

def postvar(sum2, n, a, b):
 
    return (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)


if __name__ == "__main__":
    # NOTE: run this first to get the bladder batch stuff written to files.
    """   
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
   """

    pheno = pa.read_table('bladder-pheno.txt', index_col=0)
    dat = pa.read_table('bladder-expr.txt', index_col=0)

    mod = patsy.dmatrix("~ age + cancer", pheno, return_type="dataframe")
    import time
    t = time.time()
    ebat = combat(dat, pheno.batch, mod, "age")
    print "%.2f seconds" % (time.time() - t)
 
    print ebat.ix[:5, :5]
