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
import corner

np.random.seed(seed=0) # for reproducibility

if len(argv) > 1:
    simulations_filename = argv[1]
else:
#    simulations_filename = os.path.join('..', 'giants', 'simulations.dat')
    simulations_filename = os.path.join('..', 'forward', 'simulations.dat')

perturb_dir = 'perturb'
perturb_pattern = '.+_perturb.dat'

bname = os.path.basename(simulations_filename).split('.')[0]
dir_name = 'learn-rox'#'learn-giants'
plot_dir = os.path.join(dir_name, 'plots-'+bname)
cov_dir = os.path.join(dir_name, 'covs-'+bname)
table_dir = os.path.join(dir_name, 'tables-'+bname)

for dir_name in ['learn-giants', plot_dir, perturb_dir, cov_dir, table_dir]:
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

### Load grid of models 
data = pd.read_csv(simulations_filename, sep='\t')
exclude = "radial_velocity|mass_cc"#|r02|r10|r01"#|Dnu|dnu"
#|slope|mass_cc"#|Dnu"#|dnu"#|mass_cc"
data = data.drop([i for i in data.columns if re.search(exclude, i)], axis=1)
#data = data.drop([i for i in data.columns 
#    if re.search('nu_cross', i) and 
#   not re.search('nu_cross.1.(1|2|3)$', i)], axis=1)

#nu_max = data[['M']].values * 3090 / (data[['radius']].values**2 * np.sqrt(data[['Teff']]/5777))
#data[['nu_cross.1.1']] = data[['nu_cross.1.1']].values/data[['Dnu0']].values
#data[['nu_cross.1.2']] = data[['nu_cross.1.2']].values/data[['Dnu0']].values
#data[['nu_cross.1.3']] = data[['nu_cross.1.3']].values/data[['Dnu0']].values
#data[['nu_cross.1.1']] = data[['nu_cross.1.1']].values/nu_max
#data[['nu_cross.1.2']] = data[['nu_cross.1.2']].values/nu_max
#data[['nu_cross.1.3']] = data[['nu_cross.1.3']].values/nu_max

maxs = data.max()
mins = data.min()

forests = {}

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
    "mass_X": r"Hydrogen mass X", 
    "mass_Y": r"Helium mass Y",
    "X_surf": r"Surface hydrogen X$_{\mathrm{surf}}$", 
    "Y_surf": r"Surface helium Y$_{\mathrm{surf}}$",
    "X_c": r"Core-hydrogen X$_{\mathrm{c}}$",
    "log_g": r"Surface gravity log g (cgs)", 
    "L": r"Luminosity L$/$L$_\odot$",
    "Teff": r"Effective Temperature T$_{\mathrm{eff}}$/K",
    "mass_cc": r"Convective-core mass M$_{\mathrm{cc}}$",
    "tau_MS": r"Main-sequence age $\tau_{\mathrm{MS}}$"
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
    "X_c": r"X$_{\mathrm{c}}$",
    "log_g": r"log g", 
    "L": r"L$/$L$_\odot$",
    "Teff": r"T$_{\mathrm{eff}}$/K",
    "mass_cc": r"M$_{\mathrm{cc}}$",
    "tau_MS": r"$\tau_{\mathrm{MS}}$"
}

y_init = ['M', 'Y', 'Z', 'alpha', 'overshoot', 'diffusion']
y_curr = ['age', 'X_c', 'log_g', 'L', 'radius', 'Y_surf', 'tau_MS']

def train_regressor(data, X_columns, y_show=y_init+y_curr):
    X = data.loc[:,[i for i in X_columns if i in data.columns]]
    ys = data.loc[:, [i for i in y_show if i not in X_columns]]
    
    X = X.replace([np.inf, -np.inf], np.nan)
    ys = ys.replace([np.inf, -np.inf], np.nan)
    
    drop = X.index[X.apply(np.isnan).any(1).nonzero()]
    X = X.drop(drop)
    ys = ys.drop(drop)
    
    print('Training with', X.shape[0], 'models')
    if X.shape[0] < 10:
        print("Too few models to train forest (", X.shape[0], ")")
        return [None, None, None]
    
    #print(X)
    #print()
    #print(ys)
    #print()
    for n_trees in [256]:
    #list(range(4, 16)) + [18,20] + [2**n for n in range(4, 12)]:
    #[n for n in range(4, 64)]:#[2**n for n in range(1, 12)]:
        forest = Pipeline(steps=[
            ('forest', ExtraTreesRegressor(
                #RandomForestRegressor(
                n_estimators=n_trees, 
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

def print_star(star, predict, y_names, table_curr, table_init):
    middles = np.mean(predict, 0)
    stds = np.std(predict, 0)
    #middles, stds = weighted_avg_and_std(predict, 1/stds)
    outstr = star
    initstr = "\n" + star
    currstr = "\n" + star
    for (pred_j, name) in enumerate(y_names):
        (m, s) = (middles[pred_j], stds[pred_j])
        m, s = gumr(m, s)
        outstr += "\t%s±%s" % (m, s)
        if name in y_init:
            initstr += r" & %s $\pm$ %s" % (m, s)
        if name in y_curr:
            currstr += r" & %s $\pm$ %s" % (m, s)
    print(outstr)
    table_curr.write(currstr + r' \\')
    table_init.write(initstr + r' \\')

def write_table_headers(y_names, table_curr, table_init):
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

def plot_line(value, n, bins, plt, style):
    cands = bins > value
    if any(cands):
        height = n[np.where(bins > value)[0][0]-1]
    else:
        height = n[0]
    plt.plot((value, value), (0, height), style)

def plot_star(star, predict, y_names, out_dir=plot_dir):
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
        
        n, bins, patches = ax.hist(predict[:,pred_j], 50, normed=1, 
            histtype='stepfilled', color='white')
        
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
        
        xs = [max(0, m-4*s), m+4*s]
        
        xticks = [max(0, m-3*s), m, m+3*s]
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

def process_dir(directory=perturb_dir, perturb_pattern=perturb_pattern):
    print("Processing", directory)
    
    """
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
    """
    stars = [os.path.join(perturb_dir, directory, f) 
             for f in os.listdir(os.path.join(perturb_dir, directory)) 
             if re.match(perturb_pattern, f)]
    
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
    write_headers = True
    
    wrong_cols = []
    outside = []
    run_times = []
    #forest = None
    for star_fname in stars:
        print(star_fname)
        star = os.path.split(star_fname)[-1].split("_")[0]
        star_data = pd.read_csv(star_fname, sep='\t')
        star_data = star_data.drop([i for i in star_data.columns 
            if re.search(exclude, i)], axis=1).dropna()
        
        #star_data = star_data.drop([i for i in data.columns 
        #                            if re.search('Teff', i)], axis=1)
        
        #Dnu = star_data[['Dnu0']].values
        #nu_max = star_data[['nu_max']].values #Dnu**(1/0.77)
        #if 'nu_cross.1.1' in star_data: 
        #    star_data[['nu_cross.1.1']] = star_data[['nu_cross.1.1']].values/nu_max
        #if 'nu_cross.1.2' in star_data:
        #    star_data[['nu_cross.1.2']] = star_data[['nu_cross.1.2']].values/nu_max
        #if 'nu_cross.1.3' in star_data:
        #    star_data[['nu_cross.1.3']] = star_data[['nu_cross.1.3']].values/nu_max
        
        X_columns = star_data.columns
        forest_name = tuple(sorted([i for i in X_columns if i in data.columns]))
        if forest_name in forests.keys():
            forest = forests[forest_name]
            if forest is None: continue
        else: # train the regressor
            print("Training a forest for", forest_name)
            out = train_regressor(data, sorted(X_columns))
            
            forest, y_names, X_names = out
            forests[forest_name] = forest
            
            if forest is None: continue
            
            plot_importances(forest, cov_dir, basename, X_names)
            if write_headers:
                write_table_headers(y_names, table_curr, table_init)
                write_headers = False
        
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
        
        ## predict values
        star_X = star_data.loc[:,X_names]
        start = time()
        predict = forest.predict(star_X)
        #estimates = [estimator.predict(star_X) 
        #             for estimator in forest.steps[0][1].estimators_]
        #stds = np.std(estimates, 0)
        end = time()
        run_times += [end-start]
        np.savetxt(os.path.join(cov_subdir, star+'.dat'), predict,
            header=" ".join(y_names), comments='')
        print_star(star, predict, y_names, table_curr, table_init)
        
        #plot_star(star, predict, y_names, out_dir)#, y_central)
    
    table_curr.close()
    table_init.close()
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
#process_dir()

process_dir('legacyRox')

#for mixed_test in ['x', '1', '(1|2)', '(1|2|3)']:
for mixed_test in [0,1]:
    continue
    print("Starting mixed test", mixed_test)
    bname = os.path.basename(simulations_filename).split('.')[0]
    cov_dir = os.path.join('learn-giants', 'covs-'+bname+str(mixed_test))
    
    for dir_name in ['learn-giants', plot_dir, perturb_dir, cov_dir, table_dir]:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
    
    ### Load grid of models 
    data = pd.read_csv(simulations_filename, sep='\t')
    exclude = "radial_velocity|mass_cc"
    data = data.drop([i for i in data.columns if re.search(exclude, i)], axis=1)
    data = data.drop([i for i in data.columns 
        if re.search('.3$', i) or re.search('.2$', i)], axis=1)
    if mixed_test:
        data = data.drop([i for i in data.columns if re.search('mm', i)], axis=1)
    #data = data.drop([i for i in data.columns 
    #    if re.search('nu_cross', i) and 
    #   not re.search('nu_cross.1.'+mixed_test+'$', i)], axis=1)
    
    #Dnu = data[['Dnu0']].values
    #nu_max = data[['M']].values * 3090 / (data[['radius']].values**2 * np.sqrt(data[['Teff']]/5777))
    #if 'nu_cross.1.1' in data: 
    #    data[['nu_cross.1.1']] = data[['nu_cross.1.1']].values/nu_max
    #if 'nu_cross.1.2' in data:
    #    data[['nu_cross.1.2']] = data[['nu_cross.1.2']].values/nu_max
    #if 'nu_cross.1.3' in data:
    #    data[['nu_cross.1.3']] = data[['nu_cross.1.3']].values/nu_max
    
    maxs = data.max()
    mins = data.min()
    
    forests = {}
    
    process_dir('sg-hares-mesa')



#process_dir('procyon')
#for directory in [f for f in os.listdir(perturb_dir) 
#                  if not re.match(perturb_pattern, f)]:
#    process_dir(directory)
#process_dir('Dnu')
#process_dir('legacy')
#process_dir('procyon')
