import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
#import pandas as pd 

n, l, nu, dnu = np.loadtxt('../regression/data/16CygB-freqs.dat', skiprows=1).T

def normalize(x):
    return (x-np.min(x))/(np.max(x)-np.min(x))

plt.figure()
xs = normalize(nu%234)
ys = normalize(nu)
plt.plot(xs, ys, 'ro')
plt.savefig('echelle.pdf')


ys = normalize(nu)
def ce(Dnu, xbins=10, ybins=10): 
    xs = normalize(nu%Dnu)
    bins, *_ = np.histogram2d(xs, ys, [xbins, ybins], [[0,1], [0,1]])
    size = len(xs)
    return np.sum((lambda p: p * np.log(np.sum(bins[i,:]) / size / p) \
                                 if p > 0 else 0)(bins[i][j] / size)
                  for i in np.arange(0, xbins)
                  for j in np.arange(0, ybins)) if size > 0 else np.PINF


#ces = [ce(Dnu) for Dnu in Dnus]

import multiprocessing

pool = multiprocessing.Pool(processes=16)
Dnus = np.arange(0.001, 300, 0.001)
ces = pool.map(ce, Dnus)

np.savetxt(np.array([Dnus, ces]), 'ces.dat')

plt.figure()
plt.plot(Dnus, ces, 'r-')
plt.savefig('ces.pdf')





from scipy.signal import lombscargle
time, _, _, flux, *_ = np.loadtxt('data/16CygB/kplr100002742-2011303113607_slc_wg1.dat').T
time = time[np.isfinite(flux)]
flux = flux[np.isfinite(flux)]
scaled_flux = (flux - flux.mean()) / flux.std()
freqs = np.arange(1, 5000, 0.1)
lsp = lombscargle(time*(60*60*24)/10**6, scaled_flux, freqs/(2*np.pi))

plt.figure()
plt.plot(freqs, lsp)
plt.savefig('lsp.pdf')




def CE(period, data, xbins=10, ybins=5):
    if period <= 0:
        return np.PINF

    r = rephase(data, period)
    bins, *_ = np.histogram2d(r[:,0], r[:,1], [xbins, ybins], [[0,1], [0,1]])
    size = r.shape[0]

# The following code was once more readable, but much slower.
# Here is what it used to be:
# -----------------------------------------------------------------------
#    return np.sum((lambda p: p * np.log(np.sum(bins[i,:]) / size / p) \
#                             if p > 0 else 0)(bins[i][j] / size)
#                  for i in np.arange(0, xbins)
#                  for j in np.arange(0, ybins)) if size > 0 else np.PINF
# -----------------------------------------------------------------------
# TODO: replace this comment with something that's not old code
    if size > 0:
        # bins[i,j] / size
        divided_bins = bins / size
        # indices where that is positive
        # to avoid division by zero
        arg_positive = divided_bins > 0

        # array containing the sums of each column in the bins array
        column_sums = np.sum(divided_bins, axis=1) #changed 0 by 1
        # array is repeated row-wise, so that it can be sliced by arg_positive
        column_sums = np.repeat(np.reshape(column_sums, (xbins,1)), ybins, axis=1)
        #column_sums = np.repeat(np.reshape(column_sums, (1,-1)), xbins, axis=0)


        # select only the elements in both arrays which correspond to a
        # positive bin
        select_divided_bins = divided_bins[arg_positive]
        select_column_sums  = column_sums[arg_positive]

        # initialize the result array
        A = np.empty((xbins, ybins), dtype=float)
        # store at every index [i,j] in A which corresponds to a positive bin:
        # bins[i,j]/size * log(bins[i,:] / size / (bins[i,j]/size))
        A[ arg_positive] = select_divided_bins \
                         * np.log(select_column_sums / select_divided_bins)
        # store 0 at every index in A which corresponds to a non-positive bin
        A[~arg_positive] = 0

        # return the summation
        return np.sum(A)
    else:
        return np.PINF
