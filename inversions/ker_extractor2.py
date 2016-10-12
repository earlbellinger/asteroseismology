import sys
import numpy as np
import os
import matplotlib as mpl 
mpl.use('Agg')
from matplotlib import pyplot as plt 
from kertools2 import load_fgong, kernel


def main(arguments):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                        help='directory containing eigenfunctions',
                        type=str)
    parser.add_argument('-f', '--frequencies', required=True,
                        help='file containing oscillation frequencies', 
                        type=str)
    parser.add_argument('-g', '--fgong', required=True,
                        help='FGONG file', 
                        type=str)
    parser.add_argument('-o', '--output', default='',
                        help='output directory',
                        type=str)
    parser.add_argument('-n', '--normalize', 
                        help="divide radius by max(radius)",
                        default=False, action='store_true')
    parser.add_argument('-d', '--suppress-dat', 
                        help="don't save kernel files", 
                        default=False, action='store_true')
    parser.add_argument('-p', '--suppress-plot', 
                        help="don't save kernel plots",
                        default=False, action='store_true')
    args = parser.parse_args(arguments)
    parse_dir(directory=args.input,
              freqs_fname=args.frequencies,
              fgong_fname=args.fgong,
              output_dir=args.output,
              normalize=args.normalize,
              save_dat=(not args.suppress_dat),
              save_plot=(not args.suppress_plot))


def parse_dir(directory, freqs_fname, fgong_fname, 
              output_dir='', normalize=False, 
              save_dat=True, save_plot=True, full=False):
    if (save_dat or save_plot) and output_dir != '' and \
        not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    if save_plot:
        mpl.rc('font', family='serif') 
        mpl.rc('text', usetex='true') 
        mpl.rc('text', dvipnghack='true') 
        mpl.rcParams.update({'font.size': 18}) 
    
    freqs = np.loadtxt(freqs_fname, usecols=range(3))
    fgong = load_fgong(fgong_fname, N=16)
    R = fgong['glob'][1]
    print('R =', R)
    
    r = fgong['var'][::-1,0]
    P = fgong['var'][::-1,3]                 # pressure
    rho = fgong['var'][::-1,4]               # density
    Gamma1 = fgong['var'][::-1,9]            # first adiabatic index
    cs2 = Gamma1*P/rho                       # square of the sound speed
    #np.savetxt('c2_rho.dat', np.vstack((r/R, cs2, rho)).T)
    
    for filename in list(os.walk(directory))[0][2]:
        if filename[-4:] == '.dat':
            # filename format: eig_l=ELL_n=N.dat
            parts = filename.split('_')
            ell = int(parts[1].split('=')[1])
            n = int(parts[2].split('=')[1].split('.')[0])
            idx = np.logical_and(freqs[:,0]==ell, freqs[:,1]==n)
            
            if not(any(idx)): 
                print("Could not find mode l =", ell, "n =", n)
                continue
            nu = freqs[idx][0][2] * 1e-6 # convert to Hz
            print("l =", ell, "; n =", n, '; nu =', nu)
            
            eig = np.loadtxt(os.path.join(directory, filename))
            x = eig[:,0]
            if normalize: x = x/max(x)
            
            kernels = kernel(ell, nu, eig, fgong)
            for var1, var2 in kernels:
                K1, K2 = kernels[var1, var2]
                save(x, K1, K2, R, var1, var2, ell, n,
                    normalize, output_dir, save_dat, save_plot)


def save(x, K1, K2, R, var1, var2, l, n, normalize, 
        output_dir, save_dat, save_plot):
    out_fname = '%s-%s_l=%d_n=%d'%(var1,var2,l,n)
    if save_dat:
        np.savetxt(os.path.join(output_dir, out_fname+'.dat'), 
            np.vstack((x, K1, K2)).T)
    if save_plot:
        varnames = {'c': r'c', 'c2': r'c^2', 'Gamma1': r'\Gamma_1', 'u': r'u',
                    'rho': r'\rho', 'Y': r'Y', 'psi': r'\psi'}
        plt.axhspan(0, 0, ls='dashed')
        latex = (varnames[var1], varnames[var2])
        plt.plot(x, R*K1, 'r-', label="$RK_{%s,%s}$"%latex)
        plt.plot(x, R*K2, 'b-', label="$RK_{%s,%s}$"%latex[::-1])
        plt.xlabel(r"r/R")
        plt.ylabel(r"$RK^{(n,\ell)}$")
        plt.title(r'Kernel for the $\ell=%d,\;n=%d$ mode'%(l,n))
        plt.legend(loc='upper left')
        #plt.ylim([min(0, min(R*K1), min(R*K2))-0.1,
        #          max(1, max(R*K1), max(R*K2))+0.1])
        plt.ylim([-5, 5])
        if normalize:
            plt.xlim([0, 1.01])
        else:
            plt.xlim([0, max(x)])
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, out_fname+'.png'))
        plt.close()


if __name__ == '__main__':
    exit(main(sys.argv[1:]))

