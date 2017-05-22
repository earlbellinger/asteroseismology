import sys
import struct
import numpy as np
import os
import matplotlib as mpl 
from matplotlib import pyplot as plt 

def main(arguments):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                        help='eigenfunction file (.amde)',
                        type=str)
    parser.add_argument('-o', '--output', default='',
                        help='output directory',
                        type=str)
    parser.add_argument('-n', '--normalize', 
                        help="divide radius by max(radius)",
                        default=False, action='store_true')
    parser.add_argument('-d', '--suppress-dat', 
                        help="don't save eigenfunction files",
                        default=False, action='store_true')
    parser.add_argument('-p', '--suppress-plot', 
                        help="don't save eigenfunction plots",
                        default=False, action='store_true')
    parser.add_argument('-f', '--full',
                        help="full set of eigenfunctions (nfmod0 = 1) or not" \
                             "(nfmod0 = 2 or 3)",
                        default=False, action='store_true')
    args = parser.parse_args(arguments)
    parse_eig(filename=args.input, 
              output_dir=args.output,
              normalize=args.normalize,
              save_dat=(not args.suppress_dat),
              save_plot=(not args.suppress_plot),
              full=args.full)

def parse_eig(filename, output_dir='', normalize=False, 
              save_dat=True, save_plot=True, full=False):
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
    
    if full:
        nfmod1(bin_file, normalize, output_dir, save_dat, save_plot)
    else:
        nfmod2or3(bin_file, normalize, output_dir, save_dat, save_plot)

def nfmod1(bin_file, normalize, output_dir, save_dat, save_plot):
    while len(bin_file)>0:
        # fortran header
        fmt = '<i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt, bin_file[:size])
        bin_file = bin_file[size:]
        
        # obtain mode summary
        fmt = '<'+50*'d'
        size = struct.calcsize(fmt)
        cs = struct.unpack(fmt, bin_file[:size])
        l = int(cs[17])
        n = int(cs[18])
        print("Extracting eigenfunction for mode l=%d, n=%d"%(l,n))
        bin_file = bin_file[size:]
        
        # obtain mesh
        fmt = '<i'
        size = struct.calcsize(fmt)
        N2, = struct.unpack(fmt, bin_file[:size])
        bin_file = bin_file[size:]
        
        # eigenfunctions
        N = 7
        fmt = '<'+N*N2*'d'
        size = struct.calcsize(fmt)
        z = np.array(struct.unpack(fmt, bin_file[:size]))
        x = z[0::N]
        y1 = z[1::N]
        y2 = z[2::N]
        y3 = z[3::N]
        y4 = z[4::N]
        y5 = z[5::N]
        y6 = z[6::N]
        bin_file = bin_file[size:]
        
        # save
        save(l, n, #np.vstack((x, y1, y2, y3, y4, y5, y6)).T,
            z.reshape(N2, N), 
            normalize, output_dir, save_dat, save_plot)
        
        # fortran footer
        fmt = '<i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt, bin_file[:size])
        bin_file = bin_file[size:]

def nfmod2or3(bin_file, normalize, output_dir, save_dat, save_plot):
    # number of mesh points and the mesh
    fmt = '<2i'
    size = struct.calcsize(fmt)
    N1, N2 = struct.unpack(fmt, bin_file[:size])
    bin_file = bin_file[size:]
    
    # read in the x-axis
    fmt = '<'+N2*'d'
    size = struct.calcsize(fmt)
    x = np.array(struct.unpack(fmt, bin_file[:size]))
    if normalize:
        x = x/max(x)
    bin_file = bin_file[size:]
    
    while len(bin_file)>4:
        # fortran header
        fmt = '<2i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt, bin_file[:size])
        bin_file = bin_file[size:]
        
        # obtain mode summary
        fmt = '<'+50*'d'
        size = struct.calcsize(fmt)
        cs = struct.unpack(fmt, bin_file[:size])
        l = int(cs[17])
        n = int(cs[18])
        print("Extracting eigenfunction for mode l=%d, n=%d"%(l,n))
        bin_file = bin_file[size:]
        
        # eigenfunctions
        N = 2
        fmt = '<'+N*N2*'d'
        size = struct.calcsize(fmt)
        z = np.array(struct.unpack(fmt, bin_file[:size]))
        y1 = z[0::N]
        y2 = z[1::N]
        bin_file = bin_file[size:]
        
        # save
        save(l, n, np.vstack((x, y1, y2)).T, 
            normalize, output_dir, save_dat, save_plot)
    
    if len(bin_file) != 4:
        print("Error: failed to parse eigenfunction")

def save(l, n, eigenf, normalize, output_dir, save_dat, save_plot):
    out_fname = 'eig_l=%d_n=%d'%(l,n)
    if save_dat:
        np.savetxt(os.path.join(output_dir, out_fname+'.dat'), eigenf)
    if save_plot:
        plt.axhspan(0, 0, ls='dashed')
        plt.plot(eigenf[:,0], eigenf[:,-2], 'r-', label='radial ($y_1$)')
        plt.plot(eigenf[:,0], eigenf[:,-1], 'b-', label='horizontal ($y_2$)')
        plt.xlabel("r/R")
        plt.ylabel("Displacement")
        plt.legend(loc='upper center')
        plt.title('Eigenfunctions for the $\ell=%d,\;n=%d$ mode'%(l,n))
        plt.ylim([min(0, min(eigenf[:,-2]), min(eigenf[:,-1]))-0.1,
                  max(1, max(eigenf[:,-2]), max(eigenf[:,-1]))+0.1])
        if normalize:
            plt.xlim([0, 1.01])
        else:
            plt.xlim([0, max(eigenf[:,0])])
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, out_fname+'.png'))
        plt.close()

if __name__ == '__main__':
    exit(main(sys.argv[1:]))

