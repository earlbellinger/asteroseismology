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

#import sys
#sys.stdout = open('subsets/outfile.dat', 'w')

np.random.seed(seed=0) # for reproducibility

if len(argv) > 1:
    simulations_filename = argv[1]
else:
    simulations_filename = os.path.join('..', 'forward', 'simulations.dat')

out_dir = os.path.join('subsetter')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

### Load grid of models 
raw = pd.read_csv(simulations_filename, sep='\t')
exclude = "nu_max|radial_velocity"#|Dnu|dnu"#|Dnu_|dnu|slope"
#|mass_cc"#|H|mass|X|surf"#|H|He"
data = raw.drop([i for i in raw.columns if re.search(exclude, i)], axis=1)
#data = data[data.M<1.2] ## only low mass

Xs = ["Teff", "Fe/H", "log_g", "Dnu0", "dnu02", "r02", "r10", "r01"]
ys = ['M', 'Y', 'Z', 'alpha', 'overshoot', 'diffusion',
      'age', 'X_c', 'mass_cc', 'Y_surf', 'L', 'radius']

num_trials = 25
points_per_track = sum(data['M']==data.loc[0][0])
indices = np.arange(0, len(data), points_per_track)
num_tracks = 2**(int(np.log2(len(data)/points_per_track)))

#forest = ExtraTreesRegressor(#RandomForestRegressor(#
#    n_estimators=128, n_jobs=62, oob_score=True, bootstrap=True)

def get_forest(X_names=Xs, y_names=ys, num_trees=128, data=data):
    forest = ExtraTreesRegressor(#RandomForestRegressor(#
        n_estimators=num_trees, n_jobs=32, bootstrap=True)
    X = data.loc[:, [i for i in X_names]]
    y = data.loc[:, [i for i in y_names]]
    start = time()
    rfr = forest.fit(X, y)#np.ravel(y)))
    end = time()
    return(rfr, end-start)

def train_test_set(n_tracks=num_tracks, m_points=points_per_track):
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

def percent_difference_score(y_true, y_pred):
    #print(abs(y_true - y_pred) / y_true)
    keep = np.logical_and(y_true != 0, y_pred != 0)
    if sum(keep) == 0: return 0
    y_true = y_true[keep]
    y_pred = y_pred[np.where(keep)]
    return np.median(abs(y_true - y_pred) / y_pred)

def accuracy_per_precision_score(y_true, y_pred, y_stds):
    if np.any(y_stds == 0):
        for y_std in np.where(y_stds == 0)[0]:
            if y_true.iloc[y_std] == y_pred[y_std]:
                y_stds[y_std] = 1
    return np.median(abs(y_true - y_pred) / y_stds)

def absolute_difference_score(y_true, y_pred):
    return np.median(abs(y_true - y_pred))

def make_tests(var_name, rfr, validation):
    Xs_test = validation.loc[:, [x for x in Xs]]
    ys_true = validation.loc[:, [y for y in ys]]
    ys_pred = rfr.predict(Xs_test)
    
    estimates = [estimator.predict(Xs_test) for estimator in rfr.estimators_]
    #ys_pred = np.mean(estimates, 0)
    ys_stds = np.std(estimates, 0)
    
    result = ""
    for y_i in range(len(ys)):
        y_pred = ys_pred[:,y_i]
        y_stds = ys_stds[:,y_i]
        y_true = ys_true[ys[y_i]]
        
        result += str(var_name) + " " + ys[y_i] + " " + \
            str(r2_score(y_true, y_pred)) + " " + \
            str(explained_variance_score(y_true, y_pred)) + " " + \
            str(accuracy_per_precision_score(y_true, y_pred, y_stds)) + " " + \
            str(absolute_difference_score(y_true, y_pred)) + "\n"
        
    return(result)

def write(result, f):
    print(result)
    f.write(result)
    f.flush()
    os.fsync(f.fileno())

### try with different amounts of tracks
fname = os.path.join(out_dir, "num_tracks.dat")
f = open(fname, 'w') 
header = 'num_tracks variable r2 ev sigma diff'
print(header)
f.write(header + "\n")
for n_tracks in [2**n for n in range(1, 
        int(np.log2(len(data)/points_per_track))+1)]:
    for trial_i in range(num_trials):
        (tracks, validation) = train_test_set(n_tracks=n_tracks)
        (rfr, train_time) = get_forest(data=tracks)
        result = make_tests(n_tracks, rfr, validation)
        write(result, f)
f.close()

### try with different models per track
fname = os.path.join(out_dir, "num_points.dat")
f = open(fname, 'w') 
header = 'num_points variable r2 ev sigma diff'
print(header)
f.write(header + "\n")
for m_points in [2**n for n in range(2, 1+int(np.log2(points_per_track)))]:
    for trial_i in range(num_trials):
        (tracks, validation) = train_test_set(m_points=m_points)
        (rfr, train_time) = get_forest(data=tracks)
        result = make_tests(m_points, rfr, validation)
        write(result, f)
f.close()

### try with different amounts of trees 
fname = os.path.join(out_dir, "num_trees.dat")
f = open(fname, 'w') 
header = 'num_trees variable r2 ev sigma diff'
print(header)
f.write(header + "\n")
for num_trees in [2**n for n in range(0, 8)]:
    for trial_i in range(num_trials):
        (tracks, validation) = train_test_set()
        (rfr, train_time) = get_forest(data=tracks, num_trees=num_trees)
        result = make_tests(num_trees, rfr, validation)
        write(result, f)
f.close()

### try with different amounts of trees for training time
#print('num_trees train_time')
#for num_trees in [2**n for n in range(0, 8)]:
#    for trial_i in range(num_trials):
#        (tracks, validation) = train_test_set()
#        (rfr, train_time) = get_forest(data=tracks, num_trees=num_trees)
#        print(num_trees, train_time)

