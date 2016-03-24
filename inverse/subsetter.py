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
raw = pd.read_csv(simulations_filename, sep='\t')
exclude = "nu_max|radial_velocity"#|Dnu|dnu"#|Dnu_|dnu|slope"
#|mass_cc"#|H|mass|X|surf"#|H|He"
data = raw.drop([i for i in raw.columns if re.search(exclude, i)], axis=1)

Xs = ["Teff", "Fe/H", "log_g",
      "Dnu0_median", "dnu02_median", 
      "r_sep02_median", "r_avg10_median", "r_avg01_median"]
y_init = ['M', 'Y', 'Z', 'alpha', 'overshoot', 'diffusion']
y_curr = ['age', 'X_c', 'mass_cc', 'Y_surf', 'L', 'radius']

num_trials = 20
points_per_track=sum(data['M']==data.loc[0][0])

forest = ExtraTreesRegressor(#RandomForestRegressor(#
    n_estimators=128, n_jobs=62, oob_score=True, bootstrap=True)

def get_forest(X_names=Xs, y_names=y_init+y_curr, num_trees=128, data=data):
    if num_trees != 128:
        forest = ExtraTreesRegressor(#RandomForestRegressor(#
            n_estimators=num_trees, n_jobs=62, oob_score=True, bootstrap=True)
    X = data.loc[:, [i for i in X_names]]
    ys = data.loc[:, [i for i in y_names]]
    return(forest.fit(X, np.ravel(ys)))

def train_test_set(n_tracks=4096, m_points=32,
        indices=np.arange(0, len(data), points_per_track)):
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


### try with different amounts of tracks 
points_per_track = sum(data['M']==data.loc[0][0])
indices = np.arange(0, len(data), points_per_track)
print('n_tracks m_points variable cv oob')
for n_tracks in [2**n for n in range(2, 
        int(np.log2(len(data)/points_per_track))+1)]:
        
        # grab 32 points from each track 
        m_points = 32
        for trial_i in range(num_trials):
            
            # pick out n+1 track ranges
            ranges = np.floor(np.linspace(0, len(indices), n_tracks+1))
            
            # pick n points from those ranges
            track_idxs = [np.random.randint(ranges[i], ranges[i+1]) 
                          for i in range(len(ranges)-1)]
            
            # build validation set without these tracks
            validation = data.copy()
            drop_idxs = []
            for i in track_idxs:
                drop_idxs += list(range(indices[i], 
                                        indices[i]+points_per_track))
            validation = validation.drop(validation.index[drop_idxs])
            
            tracks = pd.DataFrame()
            for idx in track_idxs:
                track = data[indices[idx]:(indices[idx]+points_per_track)]
                subset = track.iloc[np.floor(
                    np.linspace(0, points_per_track-1, m_points))]
                tracks = tracks.append(subset)
            for yy in y_init+y_curr:
                rfr = get_forest(data=tracks, y_names=[yy])
                X_test = validation.loc[:, [i for i in Xs]]
                y_true = validation.loc[:, yy] #[i for i in y_init+y_curr]]
                y_pred = rfr.predict(X_test)
                print(n_tracks, m_points, yy, 
                    r2_score(y_true, y_pred), 
                    rfr.oob_score_)


### try with different models per track
points_per_track = sum(data['M']==data.loc[0][0])
indices = np.arange(0, len(data), points_per_track)
print('n_tracks m_points variable cv oob')
#[2048]:#128, 256, 512, 1024, 2048]:#
n_tracks = 4096
for m_points in range(2, points_per_track+1):
    for trial_i in range(num_trials):
        # pick out n+1 track ranges
        ranges = np.floor(np.linspace(0, len(indices), n_tracks+1))
        
        # pick n points from those ranges
        track_idxs = [np.random.randint(ranges[i], ranges[i+1]) 
                      for i in range(len(ranges)-1)]
        
        # build validation set without these tracks
        validation = data.copy()
        drop_idxs = []
        for i in track_idxs:
            drop_idxs += list(range(indices[i], 
                                    indices[i]+points_per_track))
        validation = validation.drop(validation.index[drop_idxs])
        
        tracks = pd.DataFrame()
        for idx in track_idxs:
            track = data[indices[idx]:(indices[idx]+points_per_track)]
            subset = track.iloc[np.floor(
                np.linspace(0, points_per_track-1, m_points))]
            tracks = tracks.append(subset)
        for yy in y_init+y_curr:
            rfr = get_forest(data=tracks, y_names=[yy])
            X_test = validation.loc[:, [i for i in Xs]]
            y_true = validation.loc[:, yy] #[i for i in y_init+y_curr]]
            y_pred = rfr.predict(X_test)
            print(n_tracks, m_points, yy, 
                r2_score(y_true, y_pred), 
                rfr.oob_score_)


### try with different amounts of trees 
points_per_track = sum(data['M']==data.loc[0][0])
indices = np.arange(0, len(data), points_per_track)
print('n_tracks m_points num_trees variable cv oob')
n_tracks = 4096
m_points = 32
for num_trees in [2**n for n in range(0, 12)]:
    for trial_i in range(num_trials):
        
        # train and test 
        for yy in y_init+y_curr:
            rfr = get_forest(data=tracks, y_names=[yy], num_trees=num_trees)
            X_test = validation.loc[:, [i for i in Xs]]
            y_true = validation.loc[:, yy] #[i for i in y_init+y_curr]]
            y_pred = rfr.predict(X_test)
            print(n_tracks, m_points, num_trees yy, 
                r2_score(y_true, y_pred), 
                rfr.oob_score_)

# ### try with different combinations of number of tracks and models per track
# points_per_track = sum(data['M']==data.loc[0][0])
# indices = np.arange(0, len(data), points_per_track)
# print('n_tracks m_points variable cv oob')
# #[2048]:#128, 256, 512, 1024, 2048]:#
# for n_tracks in [2**n for n in range(8, 
        # int(np.log2(len(data)/points_per_track))+1)]:
        
        # # grab m points from each track 
        # for m_points in [2**n for n in range(2, 
                    # int(np.log2(points_per_track)+1))]:
            # #np.floor(np.linspace(5, points_per_track, 10)):
            # #cv_scores = {}
            # #oob_scores = {}
            # for trial_i in range(num_trials):
                
                # # pick out n+1 track ranges
                # ranges = np.floor(np.linspace(0, len(indices), n_tracks+1))
                
                # # pick n points from those ranges
                # track_idxs = [np.random.randint(ranges[i], ranges[i+1]) 
                              # for i in range(len(ranges)-1)]
                
                # # build validation set without these tracks
                # validation = data.copy()
                # drop_idxs = []
                # for i in track_idxs:
                    # drop_idxs += list(range(indices[i], 
                                            # indices[i]+points_per_track))
                # validation = validation.drop(validation.index[drop_idxs])
                
                # tracks = pd.DataFrame()
                # for idx in track_idxs:
                    # track = data[indices[idx]:(indices[idx]+points_per_track)]
                    # subset = track.iloc[np.floor(
                        # np.linspace(0, points_per_track-1, m_points))]
                    # tracks = tracks.append(subset)
                # for yy in y_init+y_curr:
                    # rfr = get_forest(data=tracks, y_names=[yy])
                    # X_test = validation.loc[:, [i for i in Xs]]
                    # y_true = validation.loc[:, yy] #[i for i in y_init+y_curr]]
                    # y_pred = rfr.predict(X_test)
                    # print(n_tracks, m_points, yy, 
                        # r2_score(y_true, y_pred), 
                        # rfr.oob_score_)
                    # #cv_scores[yy] += [r2_score(y_true, y_pred)]
                        # #multioutput='variance_weighted')]
                    # #oob_scores += [rfr.oob_score_]
            
            # #print(n_tracks, m_points, 
            # #      np.mean(cv_scores), np.std(cv_scores),
            # #      np.mean(oob_scores), np.std(oob_scores))


#print('\n\n')
#print('(Variable set), score')
#### try with every possible combination of inputs 
#for n_combinations in range(1, 8):
#    sets = combinations(Xs, n_combinations)
#    for var_set in sets:
#        print(var_set, get_forest(var_set).oob_score_)
