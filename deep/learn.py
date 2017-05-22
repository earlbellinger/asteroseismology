#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from time import time
import re
import numpy as np
import pandas as pd
#from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sys import argv

np.random.seed(seed=0) # for reproducibility

if len(argv) > 1:
    simulations_filename = argv[1]
else:
    simulations_filename = os.path.join('training_mat_abs.dat')

### Load grid of models 
training = pd.read_csv(simulations_filename, sep='\t')

ys = training.filter(regex=('df'))
X = training.loc[:, [i for i in training.columns if i not in ys.columns]]

drop = X.index[X.apply(np.isnan).any(1).nonzero()]
X = X.drop(drop)
ys = ys.drop(drop)

drop = ys.index[ys.apply(np.isnan).any(1).nonzero()]
X = X.drop(drop)
ys = ys.drop(drop)    

forest = Pipeline(steps=[
    ('forest', ExtraTreesRegressor(
        n_estimators=256, 
        n_jobs=32,
        oob_score=True, bootstrap=True))])
start = time()
forest.fit(X, ys)#new_ys)
end = time()
print(forest.steps[0][1].oob_score_, end-start)
print()
print("%.5g seconds to train regressor" % (end-start))
print()

test = pd.read_csv('test_mat_r.dat', sep='\t')

forest.predict(test)
np.mean(forest.predict(test), 0)[0:100]
np.std(forest.predict(test), 0)[0:100]

