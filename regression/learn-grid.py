#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from time import time
import re
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA

import matplotlib as mpl 
mpl.use("agg")
from matplotlib import pyplot as plt 
mpl.rc('font', family='serif') 
mpl.rc('text', usetex='true') 
mpl.rc('text', dvipnghack='true') 
mpl.rcParams.update({'font.size': 16}) 
import pylab as P
from sys import argv
import corner

np.random.seed(seed=0) # for reproducibility

if len(argv) > 1:
    simulations_filename = argv[1]
else:
    simulations_filename = os.path.join('..', 'grid', 'simulations.dat')

data = pd.read_csv(simulations_filename, sep='\t')
exclude = "nu_max_classic|X_surf|Z"
data = data.drop([i for i in data.columns if re.search(exclude, i)], axis=1)

maxs = data.max()
mins = data.min()

test = pd.read_csv(os.path.join('data', 'yrec-hares', 'yrec_hares.dat'), 
    sep='\t')

test = test.drop([i for i in test.columns if i not in data.columns], axis=1)

X_cols = ['Teff', 'Fe/H', 'Dnu0', 'nu_max']

X_train = data[X_cols]
#data.drop([i for i in data.columns if i not in X_cols], axis=1)
y_train = data.drop([i for i in data.columns if i in X_cols], axis=1)
X_test = test[X_cols]
#test.drop([i for i in test.columns if i not in X_cols], axis=1)
y_test = test.drop([i for i in test.columns if i in X_cols], axis=1)
y_train = data[y_test.columns]
#y_train.drop([i for i in y_train.columns if i not in y_test], axis=1)


forest = Pipeline(steps=[
    ('forest', ExtraTreesRegressor(
        #RandomForestRegressor(
        n_estimators=1024, 
        n_jobs=16,
        oob_score=True, bootstrap=True))])
start = time()
forest.fit(X_train, y_train)
end = time()
print(forest.steps[0][1].oob_score_, end-start)

unc = [[0, 0, 0, 0],
       [50, 0.05, 0.1, 50],
       [100, 0.1, 1, 100],
       [150, 0.15, 2, 150]]

for ii in range(len(unc)):
    unci = unc[ii]
    print(unci)
    X_test2 = X_test.copy()
    for jj in range(len(unci)):
        if unci[jj] > 0:
            X_test2[X_cols[jj]] = np.random.normal(X_test[X_cols[jj]], unci[jj])
    y_predict = pd.DataFrame(forest.predict(X_test2), columns=y_test.columns)
    rel_diff = (y_test-y_predict)/y_test
    rel_diff.to_csv('yrec_hares/'+str(ii)+'.dat', sep='\t')



