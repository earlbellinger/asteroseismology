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
import corner

plot_dir = 'learn_plots'
perturb_dir = 'perturb'
perturb_pattern = '.+_perturb.dat'

if not os.path.exists(perturb_dir):
    os.makedirs(perturb_dir)

### Load grid of models 
data = pd.read_csv('../forward/simulations.dat', sep='\t')
exclude = "nu_max|radial_velocity"#|H|He"
data = data.drop([i for i in data.columns if re.search(exclude, i)], axis=1)
#data = data.loc[data['M'] >= 0.8]

################################################################################
### Helpers ####################################################################
################################################################################
y_latex = {
    "M": "Mass M/M$_\\odot$", 
    "Y": "Initial helium Y$_0$", 
    "Z": "Initial metallicity Z$_0$", 
    "alpha": "Mixing length $\\alpha_{\mathrm{MLT}}$", 
    "age": "Age $\\tau$/Gyr", 
    "radius": "Radius R/R$_\\odot$", 
    "H": "Hydrogen fraction X", 
    "He": "Helium fraction Y",
    "Hc": "Core-hydrogen fraction X$_c$",
    "log_g": "Surface gravity log g/dex", 
    "L": "Luminosity L/L$_\\odot$",
    "mass_cc": "Convective core mass fraction M$_{cc}$"
}

y_latex2 = {
    "M": "M/M$_\\odot$", 
    "Y": "Y$_0$", 
    "Z": "Z$_0$", 
    "alpha": "$\\alpha_{\mathrm{MLT}}$", 
    "age": "$\\tau$/Gyr", 
    "radius": "R/R$_\\odot$", 
    "H": "X", 
    "He": "Y",
    "Hc": "X$_c$",
    "log_g": "log g/dex", 
    "L": "L/L$_\\odot$",
    "mass_cc": "M$_{cc}$"
}

y_show = ['M', 'Y', 'Z', 'age', 'radius']

def train_regressor(data, X_columns, y_exclude="median|slope|intercept"):
    X = data.loc[:,X_columns]#data.drop(list(ys.columns) + drop_cols, axis=1)
    
    ys = data.drop(X_columns, axis=1)
    drop_cols = [i for i in ys.columns 
            if re.search(y_exclude, i) if i in ys.columns]
    ys = ys.drop(drop_cols, axis=1)
    y_trfm = Pipeline(steps=[('scaler', MinMaxScaler())])#, ('pca', PCA())])
    new_ys = y_trfm.fit_transform(ys)
    
    for n_trees in [n for n in range(2,65)]:
        #forest = RandomForestRegressor(n_estimators=n_trees, n_jobs=-1,#,
        #    oob_score=True, bootstrap=True)
        forest = Pipeline(steps=[
            ('scaler', StandardScaler()), 
            #('pca', PCA()),
            ('forest', RandomForestRegressor(n_estimators=n_trees, n_jobs=-1,
                oob_score=True, bootstrap=True))])
        start = time()
        forest.fit(X, new_ys)
        end = time()
        print(n_trees, forest.steps[1][1].oob_score_, end-start)
    
    #forest = ExtraTreesRegressor(n_estimators=1000, n_jobs=-1, verbose=1,
    #    oob_score=1, bootstrap=1)
    start = time()
    forest.fit(X, new_ys)
    end = time()
    print()
    print("%.5g seconds to train regressor" % (end-start))
    #print("out-of-bag score: %.5g" % forest.oob_score_)
    print('\t\t'.join(ys.columns))
    #return [forest, ys.columns, X.columns]#, y_trfm]
    return [forest, ys.columns, X.columns, y_trfm]

def get_rc(num_ys):
    rows = (num_ys+(1 if num_ys%2 else 0)) // 2
    cols = 4 if num_ys%2 else 2
    if num_ys%3==0:
        rows = num_ys//3
        cols = num_ys//rows
    sqrt = np.sqrt(num_ys)
    if int(sqrt) == sqrt:
        rows, cols = int(sqrt), int(sqrt)
    return rows, cols, sqrt

def plot_star(star, predict, y_names, out_dir=plot_dir):
    ## Corner plot
    figure = corner.corner(
        predict[:,[i for i,y in enumerate(y_names) if y in y_show]], 
        labels=[y_latex2[y_name] for y_name in y_names
                if y_name in y_show],
        show_titles=True, title_args={"fontsize": 16})
    figure.gca().annotate(star, 
        xy=(0.5, 1.0), xycoords="figure fraction",
        xytext=(0, -5), textcoords="offset points",
        ha="center", va="top")
    plt.savefig(os.path.join(out_dir, star + '-corner.png'), dpi=400)
    plt.close()

    
    ## Regular histograms
    middles = np.median(predict, 0)
    stds = np.std(predict, 0)
    
    outstr = star
    num_ys = predict.shape[1]
    
    rows, cols, sqrt = get_rc(num_ys)
    
    plt.figure(figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.suptitle(star)
    for (pred_j, name) in enumerate(y_names[0:num_ys]):
        (m, s) = (middles[pred_j], stds[pred_j])
        outstr += "\t%.3g +/- %.3g" % (m, s)
        
        if num_ys%2 == 0 or num_ys%3==0 or int(sqrt)==sqrt:
            ax = plt.subplot(rows, cols, pred_j+1)
        elif pred_j%2 and pred_j == num_ys:
            ax = plt.subplot2grid((rows, cols), (rows, 1), colspan=2)
        else:
            ax = plt.subplot2grid((rows, cols), (pred_j//2, (pred_j%2)*2),
                colspan=2)
        
        n, bins, patches = P.hist(predict[:,pred_j], 50, normed=1, 
            histtype='stepfilled', color='white')
        y = P.normpdf(bins, m, s)
        P.plot(bins, y, 'b--', linewidth=1.5)
        
        ax.annotate(r"$\epsilon = %.3g\%%$" % (s/m*100),
            xy=(0.5, 0.125), xycoords='axes fraction', #fontsize=16,
            horizontalalignment='center', verticalalignment='center')
        
        P.xlabel(y_latex[y_names[pred_j]])
        P.locator_params(axis='x', nbins=3)
        xticks = [m-3*s, m, m+3*s]
        ax.set_xticks(xticks)
        ax.set_xlim([m-4*s, m+4*s])
        ax.set_xticklabels(['%.3g'%xtick for xtick in xticks])
        ax.set_yticklabels('',visible=False)
        ax.minorticks_on()
        plt.tight_layout()
    
    print(outstr)
    plt.subplots_adjust(top=0.9)
    plt.savefig(os.path.join(out_dir, star + '.png'), dpi=400)
    plt.close()

def process_dir(directory=perturb_dir, perturb_pattern=perturb_pattern):
    stars = [os.path.join(perturb_dir, f) 
             for f in os.listdir(perturb_dir) if re.match(perturb_pattern, f)]
    if directory != perturb_dir:
        stars = [os.path.join(perturb_dir, directory, f) 
                 for f in os.listdir(os.path.join(perturb_dir, directory)) 
                 if re.match(perturb_pattern, f)] + stars
    
    out_dir = os.path.join(plot_dir, directory)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    forest = None
    for star_fname in stars:
        star = os.path.split(star_fname)[-1].split("_")[0]
        #star_fname = os.path.join(perturb_dir, directory, star)
        star_data = pd.read_csv(star_fname, sep='\t')
        star_data = star_data.drop([i for i in star_data.columns 
            if re.search(exclude, i)], axis=1).dropna()
        
        if forest is None: # train the regressor
            out = train_regressor(data, star_data.columns)
            forest, y_names, X_names, y_trfm = out
            #forest, y_names, X_names = out#, y_trfm = out
        
        star_data = star_data.drop([i for i in star_data.columns 
            if i not in X_names], axis=1)
        #predict = forest.predict(star_data)
        predict = y_trfm.inverse_transform(forest.predict(star_data))
        plot_star(star, predict, y_names, out_dir)

################################################################################
### Start ######################################################################
################################################################################
process_dir()
for directory in [f for f in os.listdir(perturb_dir) 
                  if not re.match(perturb_pattern, f)]:
    process_dir(directory)

