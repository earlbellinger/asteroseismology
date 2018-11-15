#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import emcee 
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn import tree
from sklearn.model_selection import cross_val_score
import re
import sys
from math import *

#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import corner

np.random.seed(seed=0)

print("Loading simulations")
DF = pd.read_table('simulations.dat', sep='\t')
DF['age'] = DF['age'] / 10**9

print("Loading KIC 7970740")
KIC = pd.read_table('../regression/perturb/bigG/7970740_perturb.dat')
KIC = KIC.drop(['Dnu0', 'dnu02', 'r02', 'r01', 'r10'], axis=1)
KIC = KIC.drop([i for i in KIC.columns 
    if re.search('nu_|r01_', i)], axis=1)
KIC_columns = KIC.columns.values
KIC_columns[1] = 'Fe_H'
KIC.columns = KIC_columns

y_star = KIC.iloc[0]

Sigma = np.cov(KIC.iloc[:,2:].T)
Sigma = np.vstack([np.zeros(Sigma.shape[0]), np.zeros(Sigma.shape[0]), Sigma])
Sigma = np.insert(Sigma, 0, np.zeros(Sigma.shape[0]), axis=1)
Sigma = np.insert(Sigma, 0, np.zeros(Sigma.shape[0]), axis=1)
Sigma[0,0] = 77**2
Sigma[1,1] = 0.1**2

Sigma_inv = np.linalg.inv(Sigma)

y_columns = KIC.columns.tolist() 
X_columns = ['age', 'M', 'Y', 'Z', 'alpha', 'beta']
X_labels = [r'$\tau/$Gyr', "$M/\mathrm{M}_\odot$", "$Y_0$", "$Z_0$", 
    r'$\alpha_{\mathrm{MLT}}$', r'$\beta$']
data_ = DF.loc[:,X_columns+y_columns].dropna()

X = data_.loc[:, X_columns]
ys = data_.loc[:, [i for i in y_columns if i not in X_columns]]

#ntrees = [1, 2, 4, 8, 16, 32, 64]
#scores = [np.median(cross_val_score(RandomForestRegressor(n_estimators=n_tree),
#              X, ys, cv=10)) 
#          for n_tree in ntrees]

forest = ExtraTreesRegressor(n_estimators=8)
print(cross_val_score(forest, X, ys, cv=3))
print("Training forest")
forest.fit(X, ys)

X_max = X.max()
X_min = X.min()

def lnprior(theta):
    if any(theta > X_max) or any(theta < X_min): #X.max()
        return -np.inf
    return 0.0

def lnlike(theta):
    y = forest.predict([theta])[0]
    resid = y-y_star
    chi2 = np.dot(np.dot(resid, Sigma_inv), resid)
    return -chi2/2
    return(-np.log(chi2)/2)

def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)

ndim, nwalkers, ntemps = 6, 100, 20
sampler = emcee.PTSampler(ntemps, nwalkers, ndim, lnlike, lnprior, threads=4)
p0 = np.array(
    [[[10., 0.77, 0.26, 0.005, 1.8, 0.] + 1e-3*np.random.randn(ndim) 
      for i in range(nwalkers)]
     for j in range(ntemps)])

def loadingBar(i, N, size):
    percent = float(i) / float(N)
    sys.stdout.write("\r"
                     + str(int(i)).rjust(3, '0')
                     + "/"
                     + str(int(N)).rjust(3, '0')
                     + ' ['
                     + '='*ceil(percent*size)
                     + ' '*floor((1-percent)*size)
                     + ']')

print("Burning in MCMC sampler")
iterations = 10000
i = 1
for p, lnp, lnl in sampler.sample(p0, iterations=iterations):
    loadingBar(i, iterations, size=50)
    i = i+1
    #pass

sampler.reset()

print()
print("Running MCMC")
f = open("chain.dat", "w")
f.close()
iterations = 10000
i = 1
for p, lnp, lnl in sampler.sample(p, 
        lnprob0=lnp, lnlike0=lnl, iterations=iterations, thin=10):
    loadingBar(i, iterations, size=50)
    i = i+1
    f = open("chain.dat", "a")
    for k in range(p.shape[0]):
        f.write("{0:4d} {1:s}\n".format(k, " ".join(p[k])))
    f.close()
    #pass

print()

"""
fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
samples = sampler.chain 
labels = X_labels
for i in range(ndim):
    ax = axes[i]
    arr = sampler.chain[..., i]
    arr = arr.reshape(-1, arr.shape[-1])
    ax.plot(arr, "k", alpha=0.3)
    ax.set_xlim(0, len(arr))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");
fig.savefig('chains.pdf')
#fig.show()
"""

arr = sampler.chain
arr = arr.reshape(-1, arr.shape[-1])

fig = corner.corner(arr, labels=X_labels,
    show_titles=True,
    title_fmt='.4f',
    truths=np.percentile(arr, 50, axis=0),
    title_kwargs={"fontsize": 12})
fig.savefig('corner.pdf')
#fig.show()

list(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
    zip(*np.percentile(arr, [16, 50, 84], axis=0))))

