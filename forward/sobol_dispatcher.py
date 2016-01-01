#!usr/bin/env python

from sobol_lib import i4_sobol
import numpy as np
import subprocess
import argparse
from time import sleep
from math import log10

def main(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument('-M', default=[0.7, 1.3], nargs=2, type=float,
                        help='range of masses')
    parser.add_argument('-Y', default=[0.22, 0.34], nargs=2, type=float, 
                        help='range of helium values')
    parser.add_argument('-Z', default=[0.0001, 0.04], nargs=2, type=float,
                        help='range of metallicity values')
    parser.add_argument('-a', '--alpha', default=[1.5, 2.5], nargs=2,type=float,
                        help='range of mixing length parameter values')
    parser.add_argument('-N', default=1000, help='number of models to generate',
                        type=int)
    parser.add_argument('-s', '--skip', default=0, type=int,
                        help='offset for sobol numbers')
    parser.add_argument('-d', '--directory', default="deleter", type=str,
                        help='offset for sobol numbers')
    args = parser.parse_args(arguments)
    print(args)
    dispatch(ranges=np.vstack((args.M,args.Y,10**np.array(args.Z),args.alpha)),
             N=args.N, 
             logs=[0, 0, 1, 0], 
             directory=args.directory, 
             skip=args.skip)

def dispatch(ranges, N, logs, directory, skip=0):
    shift = ranges[:,0]
    scale = np.array([(b-a) for a,b in ranges])
    for i in range(skip, N+skip):
        vals = shift+np.array(i4_sobol(len(ranges), i)[0])*scale
        for j, val in enumerate(vals):
            if (logs[j]):
                vals[j] = log10(val)
        bash_cmd = "maybe_sub.sh -p dispatch.sh -d %s "\
            "-M %.6f -Y %.6f -Z %.6f -a %.6f"%\
            tuple([directory] + [val for val in vals])
        subprocess.Popen(bash_cmd.split(), shell=False)
        sleep(0.1)

if __name__ == '__main__':
    import sys
    exit(main(sys.argv[1:]))
