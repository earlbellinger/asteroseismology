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
#import pickle 
from sklearn.externals import joblib
import hashlib
import gc

from sklearn.cluster import MeanShift
from scipy.stats import gaussian_kde

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

drop_oob = False #True # False
enforce_bounds = True #True #False# 
always_retrain = True# False#
use_pickle = False #True #False#
make_plots = False#True# 
write_covs = True 
find_mode = False#True#
#fill = False#True

print(argv)

if len(argv) > 1:
    simulations_filename = argv[1]
else:
    simulations_filename = os.path.join('..', 'grid', 'SG_US_step.dat')
    #simulations_filename = os.path.join('..', 'forward', 'simulations.dat')
    #simulations_filename = os.path.join('..', 'grid', 'simulations.dat')
    #simulations_filename = os.path.join('classical', 'models.dat')

if len(argv) > 2:
    allname = argv[2]
else:
    allname = 'feh' #'inversion' #'kages' #

bname = os.path.basename(simulations_filename).split('.')[0]
plot_dir = os.path.join('learn-'+allname, 'plots-'+bname)
cov_dir = os.path.join('learn-'+allname, 'covs-'+bname)
table_dir = os.path.join('learn-'+allname, 'tables-'+bname)
pickle_dir = 'pickle'#os.path.join('learn-'+allname, 'pickles-'+bname)
perturb_dir = 'perturb'
perturb_pattern = '.+_perturb.dat'

if not os.path.exists('learn'):
    os.makedirs('learn')

if not os.path.exists(perturb_dir):
    os.makedirs(perturb_dir)

if not os.path.exists(cov_dir):
    os.makedirs(cov_dir)

if not os.path.exists(table_dir):
    os.makedirs(table_dir)

if use_pickle and not os.path.exists(pickle_dir):
    os.makedirs(pickle_dir)

### Load grid of models 
data = pd.read_csv(simulations_filename, sep='\t')
data['density'] = pd.Series(data['M'] / data['radius']**3, index=data.index)
#data = data.loc[data['ev_stage']<=3]
#exclude = "radial_velocity|I_mean_mag|V_mean_mag"#"nu_max|radial_velocity|mass_cc|delta_nu_asym|delta_Pg_asym"#|Dnu|dnu"
#include = r"Teff|^Fe|^nu_max$|^Dnu0$|^M$|^radius$|^age$|log_g"
#|slope|mass_cc"#|Dnu"#|dnu"#|mass_cc"

#exclude_models = '^' + '$|^'.join(['Dnu0', 'dnu02', 'r02', 'r01', 'r13', 'r10', 
#    'epsilon_p', 'undershoot', 'under_exp', 'over_exp', 'nu_max']) + '$'

#exclude_models = '^' + '$|^'.join(['Dnu0', 'dnu02', 'dnu13', 
#    'undershoot', 'under_exp', 'over_exp', 'nu_max']) + '$'

exclude_models = '|'.join(['r02_', 'r01_', 'r13_', 'r10_', 
    #'Dnu0', 'dnu02', 'dnu13', 
    #'Dnu0', 
    'epsilon_p', 
    #'undershoot', 
    'under_exp', 'over_exp', 'nu_max'])

#exclude_models = '^Dnu0$|^dnu02$|^r02$|^r01$|^dnu13$|^r13$|^r10$|epsilon_p|^over_exp$|' # "$^"
data = data.drop([i for i in data.columns 
    if re.search(exclude_models, i)], axis=1)
#data = data.drop([i for i in data.columns if not re.match(include, i)], axis=1)
data = data.replace([np.inf, -np.inf], np.nan)#.dropna()
#exclude = "radius|log_g"

#exclude_data = "^radius$|^L$|radial_velocity|I_mean_mag|V_mean_mag|" + exclude_models
exclude_data = "radial_velocity|I_mean_mag|V_mean_mag|" + exclude_models

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
    "undershoot": r"Undershoot $\alpha_{\mathrm{us}}$",
    "age": r"Age $\tau/$Gyr", 
    "radius": r"Radius R$/$R$_\odot$", 
    "mass_X": r"Hydrogen mass X", 
    "mass_Y": r"Helium mass Y",
    "X_surf": r"Surface hydrogen X$_{\mathrm{surf}}$", 
    "Y_surf": r"Surface helium Y$_{\mathrm{surf}}$",
    "X_c": r"Core-hydrogen X$_{\mathrm{c}}$",
    "log_g": r"Surface gravity log g (cgs)", 
    "L": r"Luminosity L$/$L$_\odot$",
    "Teff": r"Effective Temperature T$_{\mathrm{eff}}$/K",
    "mass_cc": r"Convective-core mass M$_{\mathrm{cc}}$",
    "log_R": r"Radius $\log($R$/$R$_\odot)$",
    "logR": r"Radius $\log($R$/$R$_\odot)$",
    "logL": r"Luminosity $\log($L$/$L$_\odot)$",
    "I_mean_mag": r"I Mag",
    "V_mean_mag": r"V Mag",
    "V_M0": r"V",
    "I_M0": r"I",
    "W": r"W",
    "density": r"Mean density $\rho$"
}

y_latex_short = {
    "M": r"M$/$M$_\odot$", 
    "Y": r"Y$_0$", 
    "Z": r"Z$_0$", 
    "alpha": r"$\alpha_{\mathrm{MLT}}$", 
    "diffusion": r"D",
    "overshoot": r"$\alpha_{\mathrm{ov}}$",
    "undershoot": r"$\alpha_{\mathrm{us}}$",
    "age": r"$\tau/$Gyr", 
    "radius": r"R$/$R$_\odot$", 
    "mass_X": r"X", 
    "mass_Y": r"Y",
    "X_surf": r"X$_{\mathrm{surf}}$", 
    "Y_surf": r"Y$_{\mathrm{surf}}$",
    "X_c": r"X$_{\mathrm{c}}$",
    "log_g": r"log g", 
    "L": r"L$/$L$_\odot$",
    "Teff": r"T$_{\mathrm{eff}}$/K",
    "mass_cc": r"M$_{\mathrm{cc}}$",
    "log_R": r"$\log($R$/$R$_\odot)$",
    "logR": r"$\log($R$/$R$_\odot)$",
    "logL": r"$\log($L$/$L$_\odot)$",
    "I_mean_mag": r"I",
    "V_mean_mag": r"V",
    "V_M0": r"V",
    "I_M0": r"I",
    "W": r"W",
    "density": r"$\rho$"
}

y_init = ['M', 'Y', 'Z', 'alpha', 'overshoot', 'diffusion']
y_curr = ['age', 'X_c', 'log_g', 'L', 'radius', 'Y_surf', 'Teff']
#y_sg = ['M', 'age', 'log_g', 'radius']
#y_classical = ['M', 'Y', 'Z', 'logR', 'logL', 'Teff', 'I_mean_mag', 'V_mean_mag']
#y_cep = ['M', 'logL', 'Teff', 'logR', 'V_M0', 'I_M0', 'W']
#y_cep = ['M', 'logL', 'Teff', 'logR', 'V_M0', 'I_M0', 'W']
#y_rrab = ['M', 'Y', 'Z', 'logL', 'logR', 'Teff', 'V_M0', 'I_M0', 'W']

y_init = ['M', 'Y', 'Z', 'alpha', 'overshoot', 'undershoot', 'diffusion']
y_curr = ['age', 'X_c', 'log_g', 'density', 'L', 'radius', 'Y_surf', 'Teff']
y_show = y_init+y_curr #y_classical # 
#y_show = y_cep

#y_show = ['M', 'Y', 'Z', 'alpha', 'age', 'radius']#, 'L']
#y_curr = y_show
#y_init = []

def train_regressor(data, X_columns, y_show=y_show, n_trees=256):
    data_ = data.loc[:,list(set(list(X_columns)+y_show))].dropna()
    X = data_.loc[:, X_columns]
    ys = data_.loc[:, [i for i in y_show if i not in X_columns]]
    
    y_names = ys.columns
    X_names = X.columns
    
    #print(X.max())
    #print(X.min())
    #print(ys.max())
    #print(ys.min())
    
    print()
    
    start = time()
    has_pickle = False 
    forest = None 
    gc.collect() 
    if use_pickle:
        pkl_filename = os.path.basename(simulations_filename) + '-' + \
            str(n_trees) + '-' + \
            '_'.join(sorted(X.columns)) + '-' + \
            '_'.join(sorted(ys.columns)) 
        hashed = hashlib.md5(pkl_filename.encode()).hexdigest()
        pkl_path = os.path.join(pickle_dir, hashed)
        print(pkl_filename)
        print(pkl_path)
        if os.path.exists(pkl_path):
            print("Loading pickle file")
            forest = joblib.load(pkl_path)
            has_pickle = True 
    if not has_pickle: 
        np.random.seed(seed=0) # for reproducibility
        gc.collect()
        forest = Pipeline(steps=[
            ('forest', ExtraTreesRegressor(
                n_estimators=n_trees, 
                n_jobs=min(n_trees, int(os.environ["OMP_NUM_THREADS"])),
                oob_score=True, bootstrap=True))])
        forest.fit(X, ys)#new_ys)
        if use_pickle:
            print("Saving pickle file")
            joblib.dump(forest, pkl_path)
    end = time()
    
    print(n_trees, forest.steps[0][1].oob_score_, end-start)
    
    print()
    print("%.5g seconds to train/load regressor" % (end-start))
    print()
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

def weighted_avg_and_std(values, weights):
    average = np.average(values, axis=0, weights=weights)
    variance = np.average((values-average)**2, axis=0, weights=weights) 
    return (average, np.sqrt(variance))

def gumr(xn, xu):
    z2 = np.trunc(np.log10(xu))+1
    z1 = np.around(xu/(10**z2), 3)
    y1 = np.around(xn*10**(-z2), 3)
    value = y1*10**z2
    uncert = z1*10**z2
    return('%g'%value, '%g'%uncert)

def print_star(star, predict, y_names, table_curr, table_init, table_ascii,
        table_modes):
    predict = predict[~np.isnan(predict).any(axis=1)]
    middles = np.mean(predict, 0)
    stds = np.std(predict, 0)
    #middles, stds = weighted_avg_and_std(predict, 1/stds)
    outstr = star
    initstr = "\n" + star
    currstr = "\n" + star
    asciistr = "\n" + star
    for (pred_j, name) in enumerate(y_names):
        (m, s) = (middles[pred_j], stds[pred_j])
        m, s = gumr(m, s)
        outstr += "\t%s±%s" % (m, s)
        asciistr += "\t%s\t%s" % (m, s)
        #strend = "\t" if pred_j < len(y_names)-1 else "\\"
        #table_ascii.write(m + "\t" + s + strend)
        if name in y_init:
            initstr += r" & %s $\pm$ %s" % (m, s)
        if name in y_curr:
            currstr += r" & %s $\pm$ %s" % (m, s)
    print(outstr)
    table_curr.write(currstr + r' \\')
    table_init.write(initstr + r' \\')
    table_ascii.write(asciistr)# + "\n")
    
    # find mode of the distribution using the mean shift algorithm 
    if find_mode: 
        mean_shift = MeanShift()
        mean_shift.fit(predict)
        modes = mean_shift.cluster_centers_
        kde = gaussian_kde(predict.T)
        mode = modes[np.argmax([kde.pdf(mode) for mode in modes])]
        table_modes.write("\n" + star + '\t' + '\t'.join(map(str, mode)))

def write_table_headers(y_names, table_curr, table_init, table_ascii, 
        table_modes):
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
    table_ascii.write('Name\t' +\
        '\t'.join([yy + '\t' 'e_'+yy for yy in y_names]))
    if find_mode:
        table_modes.write('Name\t' + '\t'.join(y_names))

def plot_line(value, n, bins, plt, style):
    cands = bins > value
    if any(cands):
        height = n[np.where(bins > value)[0][0]-1]
    else:
        height = n[0]
    plt.plot((value, value), (0, height), style)

def plot_star(star, predict, y_names, out_dir=plot_dir, nbins=10):
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
        
        n, bins, patches = ax.hist(predict[:,pred_j], nbins, normed=1, 
            histtype='stepfilled', color='white')
        
        if star == 'Sun' or star == 'Tagesstern' or star == '5774694':
            if name is 'age':
                plot_line(4.572, n, bins, plt, 'r--')
            if name is 'M':
                plot_line(1, n, bins, plt, 'r--')
        
        #if y_central is not None:
        #    mean = np.mean(y_central[:,pred_j])
        #    std = np.std(y_central[:,pred_j])
        #    
        #    plot_line(mean, n, bins, plt, 'r--')
        #    plot_line(mean+std, n, bins, plt, 'r-.')
        #    plot_line(mean-std, n, bins, plt, 'r-.')
        
        q_16, q_50, q_84 = corner.quantile(predict[:,pred_j], [0.16, 0.5, 0.84])
        q_m, q_p = q_50-q_16, q_84-q_50
        plot_line(q_50, n, bins, plt, 'k--')
        plot_line(q_16, n, bins, plt, 'k-.')
        plot_line(q_84, n, bins, plt, 'k-.')
        
        # Format the quantile display.
        fmt = "{{0:{0}}}".format(".3g").format
        title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
        title = title.format(fmt(q_50), fmt(q_m), fmt(q_p))
        
        ax.annotate(r"$\epsilon = %.3g\%%$" % (s/m*100),
            xy=(0.99, 0.12), xycoords='axes fraction',
            horizontalalignment='right', verticalalignment='right')
        
        P.xlabel(y_latex[y_names[pred_j]] + " = " + title)
        P.locator_params(axis='x', nbins=3)
        
        xs = [m-4*s, m+4*s]#[max(0, m-4*s), m+4*s]
        
        xticks = [m-3*s, m, m+3*s]#[max(0, m-3*s), m, m+3*s]
        ax.set_xticks(xticks)
        ax.set_xlim(xs)
        ax.set_xticklabels(['%.3g'%xtick for xtick in xticks])
        ax.set_yticklabels('',visible=False)
        
        ax.set_frame_on(False)
        ax.get_xaxis().tick_bottom()
        ax.axes.get_yaxis().set_visible(False)
        
        #xmin, xmax = ax1.get_xaxis().get_view_interval()
        ymin, ymax = ax.get_yaxis().get_view_interval()
        ax.add_artist(mpl.lines.Line2D(xs, (ymin, ymin), 
            color='black', linewidth=2))
        
        #ax.minorticks_on()
        plt.tight_layout()
    
    plt.savefig(os.path.join(out_dir, star + '.pdf'), dpi=400)
    plt.close()

def plot_importances(forest, cov_dir, basename, X_names):
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
    print()
    
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

def process_dir(directory=perturb_dir, perturb_pattern=perturb_pattern,
        test_sun=True):
    stars = [os.path.join(perturb_dir, f) 
             for f in os.listdir(perturb_dir) if re.match(perturb_pattern, f)]
    if not test_sun:
        stars = []
    if directory != perturb_dir:
        others = [os.path.join(perturb_dir, directory, f) 
                  for f in os.listdir(os.path.join(perturb_dir, directory)) 
                  if re.match(perturb_pattern, f)]
        names = [os.path.split(a)[-1].split("_")[0] for a in others]
        if all([x.isdigit() for x in names]):
            names = [int(x) for x in names]
        sort = [b for (a,b) in sorted(zip(names, others))]
        stars = sort + stars
    
    col_map = {}
    for star_fname in stars:
        star = os.path.split(star_fname)[-1].split("_")[0]
        star_data = pd.read_csv(star_fname, sep='\t')
        star_data.rename(columns={'Fe/H':'Fe_H'}, inplace=True)
        star_data = star_data.drop([i for i in star_data.columns 
            if re.search(exclude_data, i)], axis=1).dropna()
        star_data = star_data.drop([i for i in star_data.columns 
            if i not in data.columns], axis=1)
        if star_data.isnull().values.any() or star_data.shape[0] <= 0:
            continue
        
        if drop_oob:
            for X_name in set(star_data.columns):
                upper = maxs[X_name]
                lower = mins[X_name]
                X_vals = star_data.loc[:, X_name]
                if np.any(X_vals > upper) or np.any(X_vals < lower):
                    print("dropping", X_name)
                    star_data = star_data.drop(X_name, axis=1)
        
        X_columns = star_data.columns
        col_key = '_'.join(sorted(X_columns))
        if col_key not in col_map:
            col_map[col_key] = [star_fname]
        else:
            col_map[col_key] += [star_fname]
    
    stars = []
    sorted_col_map = sorted(col_map.keys(), 
        key=lambda k: len(col_map[k]), reverse=True)
    for key in sorted_col_map:
        stars += col_map[key]
    
    #stars = stars[:30]
    
    out_dir = os.path.join(plot_dir, directory)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    cov_subdir = os.path.join(cov_dir, directory)
    if not os.path.exists(cov_subdir):
        os.makedirs(cov_subdir)
    
    basename = os.path.basename(directory)
    table_curr_fname = os.path.join(table_dir, basename+"_curr-latex.dat")
    table_init_fname = os.path.join(table_dir, basename+"_init-latex.dat")
    table_ascii_fname = os.path.join(table_dir, basename+'.dat')
    table_modes_fname = os.path.join(table_dir, basename+'_modes.dat')
    table_curr = open(table_curr_fname, 'w')
    table_init = open(table_init_fname, 'w')
    table_ascii = open(table_ascii_fname, 'w')
    
    table_modes = None 
    if find_mode:
        table_modes = open(table_modes_fname, 'w')
    
    wrong_cols = []
    outside = []
    run_times = []
    forest = None
    X_columns_old = None
    y_names_old = None
    for star_fname in stars:
        star = os.path.split(star_fname)[-1].split("_")[0]
        star_data = pd.read_csv(star_fname, sep='\t')
        #if 'Fe/H' in star_data.columns:
        #    star_data['Fe/H']
        star_data.rename(columns={'Fe/H':'Fe_H'}, inplace=True)
        star_data = star_data.drop([i for i in star_data.columns 
            if re.search(exclude_data, i)], axis=1).dropna()
        star_data = star_data.drop([i for i in star_data.columns 
            if i not in data.columns], axis=1)
        
        #star_data = star_data.drop([i for i in data.columns 
        #                            if re.search('Teff', i)], axis=1)
        
        # check if there are any NaNs in the data 
        if star_data.isnull().values.any() or star_data.shape[0] <= 0:
            outside += [star]
            continue
        
        if drop_oob:
            for X_name in set(star_data.columns):
                upper = maxs[X_name]
                lower = mins[X_name]
                X_vals = star_data.loc[:, X_name]
                if np.any(X_vals > upper) or np.any(X_vals < lower):
                    print("dropping", X_name)
                    star_data = star_data.drop(X_name, axis=1)
        
        X_columns = star_data.columns
        different = X_columns_old is None or list(X_columns) != X_columns_old
        X_columns_old = list(X_columns)
        
        if always_retrain and different or forest is None: # train the regressor
            np.random.seed(seed=0)
            
            forest, y_names, X_names = train_regressor(data, X_columns)
            #forest, y_names, X_names = out
            
            if make_plots: 
                plot_importances(forest, cov_dir, basename, X_names)
            
            if y_names_old is None: #or list(y_names) != list(y_names_old):
                write_table_headers(y_names, table_curr, table_init, 
                    table_ascii, table_modes)
                y_names_old = y_names 
        
        ## check that it has all of the right columns
        if not set(X_names).issubset(set(star_data.columns)):
            wrong_cols += [star]
            continue
        
        ## predict values
        star_X = star_data.loc[:,X_names]
        start = time()
        predict = forest.predict(star_X)
        #estimates = [estimator.predict(star_X) 
        #             for estimator in forest.steps[0][1].estimators_]
        #stds = np.std(estimates, 0)
        end = time()
        run_times += [end-start]
        
        ## check if all the params are within the grid
        if enforce_bounds:
            #out_of_range = False
            for X_name in set(X_names):
                upper = maxs[X_name]
                lower = mins[X_name]
                X_vals = star_data.loc[:, X_name]
                if np.any(X_vals > upper) or np.any(X_vals < lower):
                    #print(X_name, "is out of range", 
                    #    sum(X_vals>upper), sum(X_vals<lower))
                    #star_data[:, X_name]
                    predict[np.where(X_vals > upper)] = np.nan
                    predict[np.where(X_vals < lower)] = np.nan
                    star_data[X_vals > upper] = np.nan
                    star_data[X_vals < lower] = np.nan
                    #print(star, X_vals.min(), lower, X_vals.max(), upper)
                    #out_of_range = True
                    #break
            #if out_of_range:
            #    outside += [(star, X_name)]
            #    continue
        
        if write_covs:
            np.savetxt(os.path.join(cov_subdir, star+'.dat'), predict,
                header=" ".join(y_names), comments='')
        
        print_star(star, predict, y_names, table_curr, table_init, table_ascii,
            table_modes)
        
        if make_plots:
            plot_star(star, predict, y_names, out_dir)#, y_central)
        
    
    table_curr.close()
    table_init.close()
    table_ascii.close()
    if find_mode:
        table_modes.close()
    print("\ntotal prediction time:", sum(run_times))
    print("time per star:", np.mean(run_times), "+/-", np.std(run_times))
    #sum(run_times)/len(run_times))
    print("time per perturbation:", np.mean(np.array(run_times) / 10000), "+/-",
         np.std(np.array(run_times) / 10000))
    #sum(run_times)/len(run_times)/10000)
    print("\ncouldn't process", wrong_cols)
    print("out of bounds", outside)

################################################################################
### Start ######################################################################
################################################################################
process_dir(allname, test_sun=False)
#process_dir()
#process_dir('feh', test_sun=True)
#process_dir('inversions', test_sun=False)
#process_dir('LMC_CEP', test_sun=False)
#process_dir('sun')
#process_dir('inversions')
#process_dir()
#process_dir('classical')
#process_dir('sg-basu')
#process_dir('inversions')
#process_dir('benard')
#process_dir('legacyRox2')
#for directory in [f for f in os.listdir(perturb_dir) 
#                  if not re.match(perturb_pattern, f)]:
#    process_dir(directory)
#process_dir('Dnu')
#process_dir('legacy')
#process_dir('gangelou')
#process_dir('kages')
#process_dir('procyon')

