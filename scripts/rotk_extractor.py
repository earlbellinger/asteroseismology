import sys
import struct
import numpy as np
import os

def main(arguments):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                        help='kernel file (.k)',
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
    parse_k(filename=args.input, 
              output_dir=args.output,
              normalize=args.normalize,
              save_dat=(not args.suppress_dat),
              save_plot=(not args.suppress_plot))

def parse_k(filename, output_dir='', normalize=False, 
              save_dat=True, save_plot=True):
    if (save_dat or save_plot) and output_dir != '' and \
        not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    if save_plot:
        import matplotlib as mpl 
        mpl.use('Agg')
        from matplotlib import pyplot as plt 
        mpl.rc('font', family='serif') 
        mpl.rc('text', usetex='true') 
        mpl.rc('text', dvipnghack='true') 
        mpl.rcParams.update({'font.size': 18}) 
    
    with open(filename, 'rb') as f:
        bin_file = f.read()
    
    x=0
    while len(bin_file)>0:
        # not sure what this is
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
        print("Extracting kernel for mode l=%d, n=%d"%(l,n))
        bin_file = bin_file[size:]
        
        # get number of grid points
        fmt = '<i'
        size = struct.calcsize(fmt)
        N2, = struct.unpack(fmt, bin_file[:size])
        bin_file = bin_file[size:]
        
        # kernel 
        N = 2
        fmt = '<'+N*N2*'d'
        size = struct.calcsize(fmt)
        z = np.array(struct.unpack(fmt, bin_file[:size]))
        x = z[0::N]
        if normalize:
            x = x/max(x)
        Knl = z[1::N]
        kernel = np.vstack((x, Knl)).T
        bin_file = bin_file[size:]
        
        # save 
        out_fname = 'l.%d_n.%d'%(l,n)
        if save_dat:
            np.savetxt(os.path.join(output_dir, out_fname), kernel[:,1],
                header=out_fname, comments='')
        if save_plot:
            plt.axhspan(0, 0, ls='dashed')
            plt.plot(x, Knl, 'r-')
            plt.xlabel("x (r/R)")
            plt.ylabel("$K_{nl}(x^n)$")
            plt.title('Kernel for the $\ell=%d,\;n=%d$ mode'%(l,n))
            plt.ylim([min(0, min(Knl))-0.1,
                      max(1, max(Knl))+0.1])
            if normalize:
                plt.xlim([0, 1.01])
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, out_fname+'.png'))
            plt.close()
        
        # not sure what this is either
        fmt = '<i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt, bin_file[:size])
        bin_file = bin_file[size:]
    
    np.savetxt(os.path.join(output_dir, 'x'), x, header='x', comments='')
    
    if len(bin_file) != 0:
        print("Error: failed to parse k file")

if __name__ == '__main__':
    exit(main(sys.argv[1:]))

