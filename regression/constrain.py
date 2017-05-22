#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from time import time
from re import search
from multiprocessing import Pool

import numpy as np
from pandas import read_csv
from sklearn.ensemble import ExtraTreesRegressor#, RandomForestRegressor
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from scipy.optimize import minimize

import matplotlib as mpl 
mpl.use("agg")
from matplotlib import pyplot as plt 
mpl.rc('font', family='serif') 
mpl.rc('text', usetex='true') 
mpl.rc('text', dvipnghack='true') 
mpl.rcParams.update({'font.size': 16}) 
import pylab as P

## Load 16 Cyg A & B
exclude_pattern = "mass|nu_max|Dnu|radial_velocity|H|diffusion|overshoot"
cygA = read_csv('perturb/16CygA_perturb.dat', sep='\t')
cygA = cygA.drop([i for i in cygA.columns 
    if search(exclude_pattern, i)], axis=1)
cygB = read_csv('perturb/16CygB_perturb.dat', sep='\t')
cygB = cygB.drop([i for i in cygB.columns 
    if search(exclude_pattern, i)], axis=1)

## Load data
data = read_csv('../forward/simulations.dat', sep='\t')
data = data.drop([i for i in data.columns
    if search(exclude_pattern, i)], axis=1)
X = data.drop(["M", "Y", "Z", "alpha", "age"],axis=1)
#X = data.drop(["M", "Y", "Z", "alpha", "diffusion", "overshoot", "age"],axis=1)
ys = data.drop(X.columns, axis=1)
X = X.drop([name for name in X.columns if name not in cygA.columns], axis=1)
print(X.head())

## Scale variables
X_trfm = MinMaxScaler()#StandardScaler()
new_Xs = X_trfm.fit_transform(X)
y_trfm = StandardScaler()#MinMaxScaler()
new_ys = y_trfm.fit_transform(ys)

# Train regressor 
"""
rfr = Pipeline(steps=[
    ('scaler', MinMaxScaler()), 
    ('rf', ExtraTreesRegressor(n_estimators=100, n_jobs=-1, verbose=1,
            oob_score=1, bootstrap=1))])
"""
start = time()
rfr = ExtraTreesRegressor(n_estimators=100, n_jobs=1, verbose=1,
    oob_score=1, bootstrap=1)
rfr.fit(new_ys, new_Xs)
end = time()
print("%.5g seconds to train regressor" % (end-start))
print(rfr.oob_score_)

## Let x be a variable defined as
## x = [Y_0, Z_0, age, M_A, alpha_A, M_B, alpha_B]
x_latex = ("$Y_0$", "$Z_0$", "\mathrm{Age}", "$M_{\mathrm{A}}$", 
    "$\\alpha_{\mathrm{A}}$", "$M_{\mathrm{B}}$", "$\\alpha_{\mathrm{B}}$")

def make_x(x0A, x0B):
    return np.hstack((x0A[1:3], x0A[4], x0A[0], x0A[3], x0B[0], x0B[3]))

def transform(A=[1.07, 0.291, 0.0253, 1.84, 6.577657], 
              B=[1.02, 0.291, 0.0253, 1.84, 6.577657]):
    x0A = y_trfm.transform(A)
    x0B = y_trfm.transform(B)
    return make_x(x0A, x0B)

def inverse_transform(x):
    A = y_trfm.inverse_transform(np.hstack((x[3], x[:2], x[4], x[2])))
    B = y_trfm.inverse_transform(np.hstack((x[5], x[:2], x[6], x[2])))
    return make_x(A, B)

x0 = transform()
"""
x0 = np.array([0.2912798, 
               0.02531405, 
               6.577657, 
               1.07, 
               1.84, 
               1.02, 
               1.84])
x0_trfm = y_trfm.transform(x0)
"""

mins = new_ys.min(axis=0)
maxs = new_ys.max(axis=0)
M_bounds = (mins[0], maxs[0])#(0.8, 1.2)
Y_bounds = (mins[1], maxs[1])#(0.22, 0.34)
Z_bounds = (mins[2], maxs[2])#(0.0001, 0.04)
A_bounds = (mins[3], maxs[3])#(1.5, 2.3)
t_bounds = (mins[4], maxs[4])#(0, 13.9)
parm_bounds = [Y_bounds, Z_bounds, t_bounds, 
               M_bounds, A_bounds, M_bounds, A_bounds]
#parm_bounds = [(0,1)]*7

start = time()

eps = np.max([np.percentile(np.diff(np.unique(new_ys[:,i])), 10) 
    for i in range(new_ys.shape[1])])
print(eps)

def objective(x, cygA_i, cygB_i):
    return sum(
      (rfr.predict(np.hstack((x[3:5], x[:3]))) - cygA_i)[0]**2 +
      (rfr.predict(np.hstack((x[5:7], x[:3]))) - cygB_i)[0]**2)

def minimization(ii):
    result = minimize(objective, x0, bounds=parm_bounds, method='SLSQP',
        options={'eps': eps, 'maxiter': 1000}, 
        args=(X_trfm.transform(cygA.loc[ii].values), 
              X_trfm.transform(cygB.loc[ii].values)))
    return result#.x

monte_carlo = Pool(None).map(minimization, range(len(cygA)))
out = np.array([inverse_transform(a.x) for a in monte_carlo if a.success])
print(sum(1 for a in monte_carlo if a.success) / len(monte_carlo))
print(out.mean(axis=0))
print(np.median(out, axis=0))
print(np.std(out, axis=0))
weights = np.array([1/a.fun for a in monte_carlo if a.success])
print(np.average(out, axis=0, weights=weights))

#out = np.array(list(Pool(None).map(minimization, range(len(cygA)))))

"""
#pool = Pool()
out = np.array(list(Pool(None).map(lambda ii: minimize(lambda x: sum(
      (rfr.predict(np.hstack((x[3:5], x[:3]))) - cygA.loc[ii].values)[0]**2 +
      (rfr.predict(np.hstack((x[5:7], x[:3]))) - cygB.loc[ii].values)[0]**2),
    x0, method="Nelder-Mead").x,
  range(5))))
"""
"""
out = []
for ii in range(len(cygA)):
    objective = lambda x: sum(
        (rfr.predict(np.hstack((x[3:5], x[:3]))) - cygA.loc[ii].values)[0]**2 +
        (rfr.predict(np.hstack((x[5:7], x[:3]))) - cygB.loc[ii].values)[0]**2)
    res = minimize(objective, x0, method="Nelder-Mead").x
    #bounds=parm_bounds, method='L-BFGS-B').x
    if out == []:
        out = res
    else:
        out = np.vstack((out, res))
"""

end = time()
print("%.3g seconds to predict Cyg A and B" % (end-start))

## Plot it!
plt.figure(figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
plt.suptitle("16 Cyg")
outstr = ""
for (col_j, name) in enumerate(x_latex):
    if (col_j == 2):
        ax = plt.subplot2grid((4,4), (0,1), colspan=2)
    else:
        jj = col_j if col_j < 2 else col_j-1
        row = int(jj/2)+1
        col = (jj%2)*2
        ax = plt.subplot2grid((4,4), (row, col), colspan=2)
    
    (m, s) = (np.median(out[:,col_j]), np.std(out[:,col_j]))
    outstr += "%s\t%.3g +/- %.3g\n" % (name, m, s)
    
    n, bins, patches = P.hist(out[:,col_j], 50, normed=1, 
    histtype='stepfilled', color='white')
    y = P.normpdf(bins, m, s)
    P.plot(bins, y, 'b--', linewidth=1.5)
    
    ax.annotate(r"$\epsilon = %.3g\%%$" % (s/m*100),
    xy=(0.5, 0.125), xycoords='axes fraction',
    horizontalalignment='center', verticalalignment='center')
    
    P.xlabel(name)
    ax.set_yticklabels('',visible=False)
    
    P.locator_params(axis='x', nbins=3)
    xticks = [m-3*s, m, m+3*s]
    ax.set_xlim([m-4*s, m+4*s])
    ax.set_xticks(xticks)
    ax.set_xticklabels(['%.3g'%xtick for xtick in xticks])
    #ax.minorticks_on()
    #"""
    plt.tight_layout()

print(outstr)
plt.subplots_adjust(top=0.9)
plt.savefig(os.path.join('learn_plots', '16Cyg-c.png'))
plt.close()

