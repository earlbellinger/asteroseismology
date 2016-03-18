#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.metrics import r2_score
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

def get_forest(X_names=Xs, y_names=y_init+y_curr, n_trees=128, data=data):
    X = data.loc[:, [i for i in X_names]]
    ys = data.loc[:, [i for i in y_names]]
    return(forest.fit(X, ys))

def get_score(X_names=Xs, y_names=y_init+y_curr, n_trees=128, data=data):
    return(get_forest(X_names, y_names, n_trees, data).oob_score_)

### try with different combinations of number of tracks and models per track
points_per_track = sum(data['M']==data.loc[0][0])
indices = np.arange(0, len(data), points_per_track)
print('n_tracks m_points score')
for n_tracks in [2**n for n in range(2, 
        int(np.log2(len(data)/points_per_track))+1)]:
    # pick out n+1 track ranges
    ranges = np.floor(np.linspace(0, len(indices), n_tracks+1))
    
    # pick n points from those ranges
    track_idxs = [np.random.randint(ranges[i], ranges[i+1]) 
                  for i in range(len(ranges)-1)]
    
    # build validation set without these tracks
    validation = data.copy()
    drop_idxs = []
    for i in track_idxs:
        drop_idxs += list(range(indices[i], indices[i]+points_per_track))
    validation = validation.drop(validation.index[drop_idxs])
    
    # grab m points from each track 
    for m_points in np.floor(np.linspace(5, points_per_track, 10)):
        tracks = pd.DataFrame()
        for idx in track_idxs:
            track = data[indices[idx]:(indices[idx]+points_per_track)]
            subset = track.iloc[np.floor(
                np.linspace(0, points_per_track-1, m_points))]
            tracks = tracks.append(subset)
        rfr = get_forest(data=tracks)
        X_test = validation.loc[:, [i for i in Xs]]
        y_true = validation.loc[:, [i for i in y_init+y_curr]]
        y_pred = rfr.predict(X_test)
        print(n_tracks, int(m_points), r2_score(y_true, y_pred, 
            multioutput='variance_weighted'))

print('\n\n')
print('(Variable set), score')
### try with every possible combination of inputs 
for n_combinations in range(1, 8):
    sets = combinations(Xs, n_combinations)
    for var_set in sets:
        print(var_set, get_score(var_set))


