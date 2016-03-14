#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from itertools import combinations
from sys import argv

np.random.seed(seed=0) # for reproducibility

if len(argv) > 1:
    simulations_filename = argv[1]
else:
    simulations_filename = os.path.join('..', 'forward', 'simulations.dat')

bname = os.path.basename(simulations_filename).split('.')[0]
table_dir = 'learn_tables-'+bname

if not os.path.exists(table_dir):
    os.makedirs(table_dir)

### Load grid of models 
data = pd.read_csv(simulations_filename, sep='\t')
exclude = "nu_max|radial_velocity|Dnu_|slope"#|mass_cc"#|H|mass|X|surf"#|H|He"
data = data.drop([i for i in data.columns if re.search(exclude, i)], axis=1)

Xs = ["Teff", "Fe/H", "log_g",
      "Dnu0_median", "dnu02_median", "r_sep02_median", 
      "r_avg10_median", "r_avg01_median"]
y_init = ['M', 'Y', 'Z', 'alpha', 'overshoot', 'diffusion']
y_curr = ['age', 'X_c', 'log_g', 'L', 'radius', 'Y_surf', 'mass_cc']

forest = ExtraTreesRegressor(n_estimators=128, 
    n_jobs=62, oob_score=True, bootstrap=True)

def get_score(X_names, y_names, n_trees=128):
    X = data.loc[:, [i for i in X_names]]
    ys = data.loc[:, [i for i in y_names]]
    forest.fit(X, ys)
    return(forest.oob_score_)

### try with different combinations of number of tracks and models per track
# implement me!

### try with every possible combination of inputs 
for n_combinations in range(1, 8):
    sets = combinations(Xs, n_combinations)
    for var_set in sets:
        print(var_set, get_score(var_set, y_init+y_curr))


