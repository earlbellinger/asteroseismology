import sys
import struct
import numpy as np
import os
import matplotlib as mpl 
from re import sub
mpl.use('Agg')
from matplotlib import pyplot as plt 

def main(arguments):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', 
                        default='basu/ker.u-Y.1',
                        help='binary file (.1)',
                        type=str)
    parser.add_argument('-o', '--output', default=None,
                        help='output directory',
                        type=str)
    parser.add_argument('-n', '--normalize', 
                        help="divide radius by max(radius)",
                        default=False, action='store_true')
    parser.add_argument('-d', '--suppress-dat', 
                        help="don't save files",
                        default=False, action='store_true')
    parser.add_argument('-p', '--suppress-plot', 
                        help="don't save plots",
                        default=False, action='store_true')
    parser.add_argument('-e', '--endianness', 
                        help="> or <", default=">")
    args = parser.parse_args(arguments)
    if args.output is None:
        args.output = sub('[^a-zA-Z/]', '', args.input)
    parse_eig(filename=args.input, 
              output_dir=args.output,
              normalize=args.normalize,
              save_dat=(not args.suppress_dat),
              save_plot=(not args.suppress_plot),
              e=args.endianness)

def parse_eig(filename, output_dir='', normalize=False, 
              save_dat=True, save_plot=True, e=">"):
    if (save_dat or save_plot) and output_dir != '' and \
        not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    if save_plot:
        mpl.rc('font', family='serif') 
        mpl.rc('text', usetex='true') 
        mpl.rc('text', dvipnghack='true') 
        mpl.rcParams.update({'font.size': 18}) 
    
    with open(filename, 'rb') as f:
        bin_file = f.read()
    
    # read(21)np,(rd(ji),ji=np,1,-1)
    # Followed by 1 line per mode
    # read(21,end=999)lll,nnn, ooo,(ak(ii),ii=np,1,-1)
    
    # number of mesh points and the mesh
    fmt = e+'2i'
    size = struct.calcsize(fmt)
    tmp, N = struct.unpack(fmt, bin_file[:size])
    bin_file = bin_file[size:]
    
    # read in the x-axis
    fmt = e+N*'f'
    size = struct.calcsize(fmt)
    x = np.array(struct.unpack(fmt, bin_file[:size]))
    if normalize:
        x = x/max(x)
    bin_file = bin_file[size:]
    
    while len(bin_file)>4:
        # fortran header
        fmt = e+'2i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt, bin_file[:size])
        bin_file = bin_file[size:]
        
        # obtain mode summary
        fmt = e+2*'i'
        size = struct.calcsize(fmt)
        l, n = struct.unpack(fmt, bin_file[:size])
        bin_file = bin_file[size:]
        
        fmt = e+'f'
        size = struct.calcsize(fmt)
        nu = struct.unpack(fmt, bin_file[:size])[0]
        bin_file = bin_file[size:]
        
        print("Extracting kernels for mode l=%d, n=%d, nu=%f"%(l,n,nu))
        
        # kernel
        fmt = e+N*'f'
        size = struct.calcsize(fmt)
        z = np.array(struct.unpack(fmt, bin_file[:size]))
        bin_file = bin_file[size:]
        
        # save
        save(l, n, np.vstack((x, z)).T, 
            normalize, output_dir, save_dat, save_plot)
    
    if len(bin_file) != 4:
        print("Error: failed to parse file")

def save(l, n, kernel, normalize, output_dir, save_dat, save_plot):
    out_fname = 'l=%d_n=%d'%(l,n)
    if save_dat:
        np.savetxt(os.path.join(output_dir, out_fname+'.dat'), kernel)
    if save_plot:
        plt.axhspan(0, 0, ls='dashed')
        plt.plot(kernel[:,0], kernel[:,1], 'r-', label='kernel')
        plt.xlabel("r/R")
        plt.ylabel("K")
        plt.legend(loc='upper center')
        plt.title('Kernels for the $\ell=%d,\;n=%d$ mode'%(l,n))
        plt.ylim([-5, 5])
        #plt.ylim([min(0, min(kernel[:,-2]), min(kernel[:,-1]))-0.1,
        #          max(1, max(kernel[:,-2]), max(kernel[:,-1]))+0.1])
        if normalize:
            plt.xlim([0, 1.01])
        else:
            plt.xlim([0, max(kernel[:,0])])
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, out_fname+'.png'))
        plt.close()

if __name__ == '__main__':
    exit(main(sys.argv[1:]))

