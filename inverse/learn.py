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
#import corner

np.random.seed(seed=0) # for reproducibility

if len(argv) > 1:
    simulations_filename = argv[1]
else:
    simulations_filename = os.path.join('..', 'forward', 'simulations.dat')

bname = os.path.basename(simulations_filename).split('.')[0]
plot_dir = 'learn_plots-'+bname
cov_dir = 'learn_covs-'+bname
table_dir = 'learn_tables-'+bname
perturb_dir = 'perturb'
perturb_pattern = '.+_perturb.dat'

if not os.path.exists(perturb_dir):
    os.makedirs(perturb_dir)

if not os.path.exists(cov_dir):
    os.makedirs(cov_dir)

if not os.path.exists(table_dir):
    os.makedirs(table_dir)

### Load grid of models 
data = pd.read_csv(simulations_filename, sep='\t')
exclude = "nu_max|radial_velocity|Dnu_|slope"#|mass_cc"#|H|mass|X|surf"#|H|He"
data = data.drop([i for i in data.columns if re.search(exclude, i)], axis=1)
#data = data.loc[data['M'] >= 0.8]

maxs = data.max()
mins = data.min()

################################################################################
### Helpers ####################################################################
################################################################################
y_latex = {
    "M": r"Mass M$/$M$_\odot$", 
    "Y": r"Initial helium Y$_0$", 
    "Z": r"Initial metallicity Z$_0$", 
    "alpha": r"Mixing length $\alpha_{\mathrm{MLT}}$", 
    "diffusion": r"Diffusion factor D",
    "overshoot": r"Overshoot $\alpha_{\mathrm{ov}}$", 
    "age": r"Age $\tau/$Gyr", 
    "radius": r"Radius R$/$R$_\odot$", 
    "mass_X": r"Hydrogen mass-fraction X", 
    "mass_Y": r"Helium mass-fraction Y",
    "X_surf": r"Surface hydrogen X$_{\mathrm{surf}}$", 
    "Y_surf": r"Surface helium Y$_{\mathrm{surf}}$",
    "X_c": r"Core-hydrogen fraction X$_c$",
    "log_g": r"Surface gravity log g (cgs)", 
    "L": r"Luminosity L$/$L$_\odot$",
    "mass_cc": r"Convective core mass-fraction M$_{\mathrm{cc}}$"
}

y_latex_short = {
    "M": r"M$/$M$_\odot$", 
    "Y": r"Y$_0$", 
    "Z": r"Z$_0$", 
    "alpha": r"$\alpha_{\mathrm{MLT}}$", 
    "diffusion": r"D",
    "overshoot": r"$\alpha_{\mathrm{ov}}$",
    "age": r"$\tau/$Gyr", 
    "radius": r"R$/$R$_\odot$", 
    "mass_X": r"X", 
    "mass_Y": r"Y",
    "X_surf": r"X$_{\mathrm{surf}}$", 
    "Y_surf": r"Y$_{\mathrm{surf}}$",
    "X_c": r"X$_{\mathrm{c}}$/M$_*$",
    "log_g": r"log g", 
    "L": r"L$/$L$_\odot$",
    "mass_cc": r"M$_{\mathrm{cc}}$"
}

#y_show = ['M', 'Y', 'Z', 'age', 'radius']
#y_show = ['M', 'Y', 'Z', 'alpha', 'diffusion', 'overshoot', 'age']
#all_ys = ['M', 'Y', 'Z', 'alpha', 'diffusion', 'overshoot', 'age', 
#    'log_g', 'L', 'radius', 'Y_surf', 'X_c', 'mass_cc']
y_init = ['M', 'Y', 'Z', 'alpha', 'overshoot', 'diffusion']
y_curr = ['age', 'X_c', 'log_g', 'L', 'radius', 'Y_surf', 'mass_cc']

def train_regressor(data, X_columns, y_show=y_init+y_curr):
    X = data.loc[:,X_columns]
    ys = data.loc[:, [i for i in y_show if i not in X_columns]]
    
    print()
    for n_trees in [1024]:
    #list(range(4, 16)) + [18,20] + [2**n for n in range(4, 12)]:
    #[n for n in range(4, 64)]:#[2**n for n in range(1, 12)]:
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
    return [forest, y_names, X_names]

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

def print_star(star, predict, y_names, table_curr, table_init):
    middles = np.mean(predict, 0)
    stds = np.std(predict, 0)
    outstr = star
    initstr = "\n" + star
    currstr = "\n" + star
    for (pred_j, name) in enumerate(y_names):
        (m, s) = (middles[pred_j], stds[pred_j])
        outstr += "\t%.3g±%.2g" % (m, s)
        if name in y_init:
            initstr += r" & %.3g $\pm$ %.2g" % (m, s)
        if name in y_curr:
            currstr += r" & %.3g $\pm$ %.2g" % (m, s)
    print(outstr)
    table_curr.write(currstr + r' \\')
    table_init.write(initstr + r' \\')

def plot_star(star, predict, y_names, out_dir=plot_dir, y_show=y_init):
    ## Corner plot
    predict_cols = [i for i,y in enumerate(y_names) if y in y_show]
    #truths = [np.NaN for i in range(len(predict_cols))]
    #if star[:3] == 'Sun' or star[:10] == 'Tagesstern':
    #    if y_show==y_init:
    #        truths[0] = 1
    #    else:
    #        truths[0] = 4.57
    figure = corner.corner(
        predict[:,predict_cols], 
        labels=[y_latex_short[y_name] for y_name in y_names
                if y_name in y_show],
        title_fmt='.2g',
        #truths=truths,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True, 
        fill_contours=True, ret=True, 
        bins=50, smooth=1.0,
        title_args={"fontsize": 12})
    figure.set_size_inches(2*6.97522, 4*4.17309)
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
    for (pred_j, name) in enumerate(y_names[0:num_ys]):
        (m, s) = (middles[pred_j], stds[pred_j])
        
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
            xy=(0.5, 0.125), xycoords='axes fraction',
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
    
    basename = os.path.basename(directory)
    table_curr_fname = os.path.join(table_dir, basename+"_curr.dat")
    table_init_fname = os.path.join(table_dir, basename+"_init.dat")
    table_curr = open(table_curr_fname, 'w')
    table_init = open(table_init_fname, 'w')
    
    wrong_cols = []
    outside = []
    run_times = []
    forest = None
    for star_fname in stars:
        star = os.path.split(star_fname)[-1].split("_")[0]
        star_data = pd.read_csv(star_fname, sep='\t')
        star_data = star_data.drop([i for i in star_data.columns 
            if re.search(exclude, i)], axis=1).dropna()
        
        if forest is None: # train the regressor
            X_columns = star_data.columns
            out = train_regressor(data, X_columns)
            forest, y_names, X_names = out
            
            ## Plot importances
            est = forest.steps[0][1]
            importances = est.feature_importances_
            indices = np.argsort(importances)[::-1]
            import_dist = np.array([tree.feature_importances_ 
                for tree in est.estimators_])
            
            np.savetxt(os.path.join(cov_dir, 'feature-importance-' + \
                    basename + '.dat'),
                import_dist, header=" ".join(X_names), comments='')
            
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
            plt.savefig(os.path.join(plot_dir, 'feature_importance-' + \
                basename + '.pdf'))
            plt.close()
            
            print()
            print("Name\t"+"\t".join(y_names))
            table_curr.write(
                r"\tablehead{\colhead{Name} & "+\
                ' & '.join([r"\colhead{" + y_latex_short[yy] + r"}"
                            for yy in y_names if yy in y_curr]) +\
                r"}\startdata" )
            table_init.write(
                r"\tablehead{\colhead{Name} & "+\
                ' & '.join([r"\colhead{" + y_latex_short[yy] + r"}"
                            for yy in y_names if yy in y_init]) +\
                r"}\startdata" )
        
        ## check that it has all of the right columns
        if not set(X_names).issubset(set(star_data.columns)):
            wrong_cols += [star]
            continue
        
        ## check if all the params are within the grid
        out_of_range = False
        for X_name in set(X_names):
            upper = maxs[X_name]
            lower = mins[X_name]
            X_vals = star_data.loc[:, X_name]
            if np.any(X_vals > upper) or np.any(X_vals < lower):
                print(X_name, "is out of range")
                print(star, X_vals.min(), lower, X_vals.max(), upper)
                out_of_range = True
                break
        if out_of_range:
            outside += [(star, X_name)]
            continue
        
        star_X = star_data.loc[:,X_names]
        start = time()
        predict = forest.predict(star_X)
        end = time()
        run_times += [end-start]
        np.savetxt(os.path.join(cov_subdir, star+'.dat'), predict,
            header=" ".join(y_names), comments='')
        print_star(star, predict, y_names, table_curr, table_init)
        #plot_star(star+"_init", predict, y_names, out_dir)
        #plot_star(star+"_curr", predict, y_names, out_dir, y_show=y_curr)
    
    table_curr.close()
    table_init.close()
    print("\ntotal prediction time:", sum(run_times))
    print("time per star:", sum(run_times)/len(run_times))
    print("time per perturbation:", sum(run_times)/len(run_times)/10000)
    print("\ncouldn't process", wrong_cols)
    print("out of bounds", outside)

################################################################################
### Start ######################################################################
################################################################################
process_dir()
for directory in [f for f in os.listdir(perturb_dir) 
                  if not re.match(perturb_pattern, f)]:
    process_dir(directory)

