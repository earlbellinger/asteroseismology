#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.metrics import r2_score, explained_variance_score
from itertools import combinations
from sys import argv
from time import time

np.random.seed(seed=0) # for reproducibility

if len(argv) > 1:
    simulations_filename = argv[1]
else:
    simulations_filename = os.path.join('..', 'forward', 'simulations.dat')

### Load grid of models 
raw = pd.read_csv(simulations_filename, sep='\t')
exclude = "nu_max|radial_velocity"#|Dnu|dnu"#|Dnu_|dnu|slope"
#|mass_cc"#|H|mass|X|surf"#|H|He"
data = raw.drop([i for i in raw.columns if re.search(exclude, i)], axis=1)

Xs = ["Teff", "Fe/H", "log_g",
      "Dnu0_median", "dnu02_median", 
      "r_sep02_median", "r_avg10_median", "r_avg01_median"]
ys = ['M', 'Y', 'Z', 'alpha', 'overshoot', 'diffusion',
      'age', 'X_c', 'mass_cc', 'Y_surf', 'L', 'radius']

num_trials = 20
points_per_track = sum(data['M']==data.loc[0][0])
indices = np.arange(0, len(data), points_per_track)

#forest = ExtraTreesRegressor(#RandomForestRegressor(#
#    n_estimators=128, n_jobs=62, oob_score=True, bootstrap=True)

def get_forest(X_names=Xs, y_names=ys, num_trees=128, data=data):
    forest = ExtraTreesRegressor(#RandomForestRegressor(#
        n_estimators=num_trees, n_jobs=62, bootstrap=True)
    X = data.loc[:, [i for i in X_names]]
    y = data.loc[:, [i for i in y_names]]
    start = time()
    rfr = forest.fit(X, y)#np.ravel(y)))
    end = time()
    return(rfr, end-start)

def train_test_set(n_tracks=8192, m_points=32):
    # pick out n+1 track ranges
    ranges = np.floor(np.linspace(0, len(indices), n_tracks+1))
    
    # pick n points from those ranges
    track_idxs = [np.random.randint(ranges[i], ranges[i+1]) 
                  for i in range(len(ranges)-1)]
    
    # get training set 
    tracks = pd.DataFrame()
    for idx in track_idxs:
        track = data[indices[idx]:(indices[idx]+points_per_track)]
        subset = track.iloc[np.floor(
            np.linspace(0, points_per_track-1, m_points))]
        tracks = tracks.append(subset)
    
    # build validation set without these tracks
    validation = data.copy()
    drop_idxs = []
    for i in track_idxs:
        drop_idxs += list(range(indices[i], 
                                indices[i]+points_per_track))
    validation = validation.drop(validation.index[drop_idxs])
    
    return (tracks, validation)

def make_tests(var_name, rfr, validation):
    Xs_test = validation.loc[:, [x for x in Xs]]
    ys_true = validation.loc[:, [y for y in ys]]
    ys_pred = rfr.predict(Xs_test)
    for y_i in range(len(ys)):
        y_pred = ys_pred[:,y_i]
        y_true = ys_true[ys[y_i]]
        print(var_name, ys[y_i], 
            r2_score(y_true, y_pred),
            explained_variance_score(y_true, y_pred))

### try with different amounts of tracks 
print('n_tracks variable cv ev')
for n_tracks in [2**n for n in range(1, 
        int(np.log2(len(data)/points_per_track))+1)]:
    for trial_i in range(num_trials):
        (tracks, validation) = train_test_set(n_tracks=n_tracks)
        (rfr, train_time) = get_forest(data=tracks)
        make_tests(n_tracks, rfr, validation)

### try with different models per track
print('m_points variable cv ev')
for m_points in [2**n for n in range(1, 1+int(np.log2(points_per_track)))]:
    for trial_i in range(num_trials):
        (tracks, validation) = train_test_set(m_points=m_points)
        (rfr, train_time) = get_forest(data=tracks)
        make_tests(m_points, rfr, validation)

### try with different amounts of trees 
print('num_trees variable cv ev')
for num_trees in [2**n for n in range(0, 8)]:
    for trial_i in range(num_trials):
        (tracks, validation) = train_test_set()
        (rfr, train_time) = get_forest(data=tracks, num_trees=num_trees)
        make_tests(num_trees, rfr, validation)

### try with different amounts of trees for training time
print('num_trees train_time')
for num_trees in [2**n for n in range(0, 8)]:
    for trial_i in range(num_trials):
        (tracks, validation) = train_test_set()
        (rfr, train_time) = get_forest(data=tracks, num_trees=num_trees)
        print(num_trees, train_time)
