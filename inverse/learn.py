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
cov_dir = 'learn_covs'
perturb_dir = 'perturb'
perturb_pattern = '.+_perturb.dat'

if not os.path.exists(perturb_dir):
    os.makedirs(perturb_dir)

if not os.path.exists(cov_dir):
    os.makedirs(cov_dir)

### Load grid of models 
data = pd.read_csv('../forward/simulations.dat', sep='\t')
exclude = "nu_max|radial_velocity|Dnu"#|H|mass|X|surf"#|H|He"
data = data.drop([i for i in data.columns if re.search(exclude, i)], axis=1)
#data = data.loc[data['M'] >= 0.8]

################################################################################
### Helpers ####################################################################
################################################################################
y_latex = {
    "M": r"Mass M$/$M$_\odot$", 
    "Y": r"Initial helium Y$_0$", 
    "Z": r"Initial metallicity Z$_0$", 
    "alpha": r"Mixing length $\alpha_{\mathrm{MLT}}$", 
    "diffusion": r"Diffusion coefficient D",
    "overshoot": r"Overshoot f",
    "age": r"Age $\tau/$Gyr", 
    "radius": r"Radius R$/$R$_\odot$", 
#    "H": r"Hydrogen fraction X", 
#    "He": r"Helium fraction Y",
#    "Hc": r"Core-hydrogen fraction X$_c$",
    "mass_X": r"Hydrogen mass-fraction X", 
    "mass_Y": r"Helium mass-fraction Y",
    "X_surf": r"Surface hydrogen X$_{\mathrm{surf}}$", 
    "Y_surf": r"Surface helium Y$_{\mathrm{surf}}$",
    "X_c": r"Core-hydrogen fraction X$_c$",
    "log_g": r"Surface gravity log g (cgs)", 
    "L": r"Luminosity L$/$L$_\odot$",
    "mass_cc": r"Convective core mass-fraction M$_{\mathrm{cc}}$"
}

y_latex2 = {
    "M": r"M$/$M$_\odot$", 
    "Y": r"Y$_0$", 
    "Z": r"Z$_0$", 
    "alpha": r"$\alpha_{\mathrm{MLT}}$", 
    "diffusion": r"D",
    "overshoot": r"f",
    "age": r"$\tau/$Gyr", 
    "radius": r"R$/$R$_\odot$", 
#    "H": r"X", 
#    "He": r"Y",
#    "Hc": r"X$_\mathrm{c}$",
    "mass_X": r"X", 
    "mass_Y": r"Y",
    "X_surf": r"X", 
    "Y_surf": r"Y",
    "X_c": r"X$_c$",
    "log_g": r"log g (cgs)", 
    "L": r"L$/$L$_\odot$",
    "mass_cc": r"M$_{\mathrm{cc}}$"
}

#y_show = ['M', 'Y', 'Z', 'age', 'radius']
#y_show = ['M', 'Y', 'Z', 'alpha', 'diffusion', 'overshoot', 'age']
#all_ys = ['M', 'Y', 'Z', 'alpha', 'diffusion', 'overshoot', 'age', 
#    'log_g', 'L', 'radius', 'Y_surf', 'X_c', 'mass_cc']
y_init = ['M', 'Y', 'Z', 'alpha', 'diffusion', 'overshoot']
y_curr = ['age', 'log_g', 'L', 'radius', 'Y_surf', 'X_c', 'mass_cc']

def train_regressor(data, X_columns, y_show=y_init+y_curr):
    X = data.loc[:,X_columns]
    ys = data.loc[:, [i for i in y_show if i not in X_columns]]
    #ys = data.drop(X_columns, axis=1).loc[:, y_show]
    
    #ys = ys.drop([i for i in y_show if re.search(ys.columns, i)], axis=1)
    #data.drop(X_columns, axis=1)
    #drop_cols = [i for i in ys.columns 
    #        if re.search(y_exclude, i) if i in ys.columns]
    #ys = ys.drop(drop_cols, axis=1)
    #y_trfm = Pipeline(steps=[('scaler', MinMaxScaler())])#, ('pca', PCA())])
    #new_ys = y_trfm.fit_transform(ys)
    
    print()
    for n_trees in [1024]:#[n for n in range(1, 2048)]:
        forest = Pipeline(steps=[
            ('forest', ExtraTreesRegressor(n_estimators=n_trees, 
                n_jobs=min(n_trees, 62),
                oob_score=True, bootstrap=True))])
        start = time()
        forest.fit(X, ys)#new_ys)
        end = time()
        print(n_trees, forest.steps[0][1].oob_score_, end-start)
    
    print()
    print("%.5g seconds to train regressor" % (end-start))
    print()
    
    y_names = ys.columns
    X_names = X.columns
    return [forest, y_names, X_names]#, y_trfm]

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

def print_star(star, predict, y_names):
    middles = np.mean(predict, 0)
    stds = np.std(predict, 0)
    outstr = star
    num_ys = predict.shape[1]
    for (pred_j, name) in enumerate(y_names[0:num_ys]):
        (m, s) = (middles[pred_j], stds[pred_j])
        outstr += r" & %.3g $\pm$ %.2g" % (m, s)
    print(outstr + r' \\')

def plot_star(star, predict, y_names, out_dir=plot_dir, y_show=y_init):
    ## Corner plot
    #truths = [np.NaN for i in range(len(y_names))]
    #if star == 'Sun' or star == 'Tagesstern':
    #    truths[0] = 1
    #    truths[6] = 4.57
    figure = corner.corner(
        predict[:,[i for i,y in enumerate(y_names) if y in y_show]], 
        labels=[y_latex2[y_name] for y_name in y_names
                if y_name in y_show],
        title_fmt='.2g',
        #truths=truths,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True, title_args={"fontsize": 16})
    plt.savefig(os.path.join(out_dir, star + '-corner.pdf'))
    plt.close()
    
    ## Regular histograms
    middles = np.mean(predict, 0)
    stds = np.std(predict, 0)
    
    #outstr = star
    num_ys = predict.shape[1]
    
    rows, cols, sqrt = get_rc(num_ys)
    
    plt.figure(figsize=(6.97522*2, 4.17309*2), dpi=400, 
        facecolor='w', edgecolor='k')
    #plt.suptitle(star)
    for (pred_j, name) in enumerate(y_names[0:num_ys]):
        (m, s) = (middles[pred_j], stds[pred_j])
        #outstr += "\t%.3g +/- %.3g" % (m, s)
        #outstr += r" & %.3g $\pm$ %.2g" % (m, s)
        
        if num_ys%2==0 or num_ys%3==0 or int(sqrt)==sqrt:
            ax = plt.subplot(rows, cols, pred_j+1)
        elif pred_j%2==0 and pred_j == num_ys-1:
            ax = plt.subplot2grid((rows, cols), (pred_j//2, 1), colspan=2)
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
        xticks = [max(0, m-3*s), m, m+3*s]
        ax.set_xticks(xticks)
        ax.set_xlim([max(0, m-4*s), m+4*s])
        ax.set_xticklabels(['%.3g'%xtick for xtick in xticks])
        ax.set_yticklabels('',visible=False)
        ax.minorticks_on()
        plt.tight_layout()
    
    #print(outstr + r' \\')
    #plt.subplots_adjust(top=0.9)
    plt.savefig(os.path.join(out_dir, star + '.png'), dpi=400)
    plt.close()

def process_dir(directory=perturb_dir, perturb_pattern=perturb_pattern):
    stars = [os.path.join(perturb_dir, f) 
             for f in os.listdir(perturb_dir) if re.match(perturb_pattern, f)]
    if directory != perturb_dir:
        others = [os.path.join(perturb_dir, directory, f) 
                  for f in os.listdir(os.path.join(perturb_dir, directory)) 
                  if re.match(perturb_pattern, f)]
        names = [os.path.split(a)[-1].split("_")[0] for a in others]
        if all([x.isdigit() for x in names]):
            names = [int(x) for x in names]
        sort = [b for (a,b) in sorted(zip(names, others))]
        stars = sort + stars
    
    out_dir = os.path.join(plot_dir, directory)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    cov_subdir = os.path.join(cov_dir, directory)
    if not os.path.exists(cov_subdir):
        os.makedirs(cov_subdir)
    
    didnt_work = []
    forest = None
    for star_fname in stars:
        star = os.path.split(star_fname)[-1].split("_")[0]
        #star_fname = os.path.join(perturb_dir, directory, star)
        star_data = pd.read_csv(star_fname, sep='\t')
        star_data = star_data.drop([i for i in star_data.columns 
            if re.search(exclude, i)], axis=1).dropna()
        
        if forest is None: # train the regressor
            X_columns = star_data.columns
            out = train_regressor(data, X_columns)
            forest, y_names, X_names = out#, y_trfm = out
            
            ## Plot importances
            est = forest.steps[0][1]
            importances = est.feature_importances_
            indices = np.argsort(importances)[::-1]
            import_dist = np.array([tree.feature_importances_ 
                for tree in est.estimators_])
            
            np.savetxt(os.path.join(cov_dir, 'feature-importance-'+star+'.dat'),
                import_dist, header=" ".join(X_names[indices]), comments='')
            
            import_dist = import_dist.T[indices][::-1].T
            
            print("Feature ranking:")
            for f in range(len(X_names)):
                print("%d. %s (%.3g)" % (f+1, X_names[indices[f]], 
                    importances[indices[f]]))
            
            mpl.rc('text', usetex='false') 
            plt.figure(figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
            plt.boxplot(import_dist, vert=0)
            plt.yticks(range(1,1+len(X_names)),
                       np.array(X_names)[indices][::-1])
            plt.xlabel("Feature importance")
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, 'feature_importance' + \
                star + '.pdf'))
            plt.close()
            
            print()
            print(r"\colhead{Name} & "+\
                  ' & '.join([r"\colhead{" + y_latex2[yy] + r"}"
                              for yy in y_names]))
            #forest, y_names, X_names = out#, y_trfm = out
        
        #print()
        
        #star_data = star_data.drop([i for i in star_data.columns 
        #    if i not in X_names], axis=1)
        if not set(X_names).issubset(set(star_data.columns)):
            didnt_work += [star]
            continue
        star_X = star_data.loc[:,X_names]
        predict = forest.predict(star_X)
        #predict = y_trfm.inverse_transform(forest.predict(star_X))
        np.savetxt(os.path.join(cov_subdir, star+'.dat'), predict,
            header=" ".join(y_names), comments='')
        print_star(star, predict, y_names)
        plot_star(star+"_init", predict, y_names, out_dir)
        plot_star(star+"_curr", predict, y_names, out_dir, y_show=y_curr)
    print("\ncouldn't process", didnt_work)

################################################################################
### Start ######################################################################
################################################################################
process_dir()
for directory in [f for f in os.listdir(perturb_dir) 
                  if not re.match(perturb_pattern, f)]:
    process_dir(directory)

