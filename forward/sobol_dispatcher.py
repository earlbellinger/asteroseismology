#!usr/bin/env python
#### Call dispatch.sh with quasi-random inputs 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

import numpy as np
import subprocess
import argparse
from time import sleep

from sys import path
path.append('../scripts')
from sobol_lib import i4_sobol

def main(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument('-M', default=[0.7, 1.6], nargs=2, type=float,
                        help='range of masses')
    parser.add_argument('-Y', default=[0.22, 0.34], nargs=2, type=float, 
                        help='range of helium values')
    parser.add_argument('-Z', default=[0.0004, 0.04], nargs=2, type=float,
                        help='range of metallicity values')
    parser.add_argument('-a', '--alpha', default=[1.5, 2.5], nargs=2,type=float,
                        help='range of mixing length parameter values')
    parser.add_argument('-D', '--diffusion', default=[0, 2], nargs=2,
                        type=float, 
                        help='range of diffusion coefficient values')
    parser.add_argument('-f', '--overshoot', default=[0, 0.5], nargs=2,
                        type=float, 
                        help='range of overshoot values')
    parser.add_argument('-N', default=1000, help='number of tracks to generate',
                        type=int)
    parser.add_argument('-s', '--skip', default=20000, type=int,
                        help='offset for sobol numbers')
    parser.add_argument('-d', '--directory', default="simulations", type=str,
                        help='offset for sobol numbers')
    parser.add_argument('-p', '--parallel', default=1, 
                        type=int, help='number of CPUs to use')
    parser.add_argument('-l', '--logs', default=[0, 0, 1, 0, 0, 0], 
                        type=list,
                        help='booleans of whether to log M, Y, Z, alpha, D, f')
    args = parser.parse_args(arguments)
    print(args)
    ranges = np.vstack((args.M, args.Y, args.Z, args.alpha,
        args.diffusion, args.overshoot))
    for i, val in enumerate(ranges):
        if (args.logs[i]):
            ranges[i] = np.log10(ranges[i])
    print(ranges)
    dispatch(ranges=ranges, N=args.N, logs=args.logs, 
        directory=args.directory, skip=args.skip, parallel=args.parallel)

def dispatch(ranges, N, logs, directory, skip=0, parallel=r"$OMP_NUM_THREADS"):
    shift = ranges[:,0]
    scale = np.array([(b-a) for a,b in ranges])
    for i in range(skip, N+skip):
        vals = shift+np.array(i4_sobol(len(ranges), i)[0])*scale
        for j, val in enumerate(vals):
            if (logs[j]):
                vals[j] = 10**val
        bash_cmd = "maybe_sub.sh -n -p %d dispatch.sh -d %s "\
            "-M %.6f -Y %.6f -Z %.6f -a %.6f -D %.6f -f %.6f -r 1"%\
            tuple([parallel, directory] + [val for val in vals])
        print(bash_cmd)
        #exit()
        subprocess.Popen(bash_cmd.split(), shell=False)
        sleep(0.05)

if __name__ == '__main__':
    import sys
    exit(main(sys.argv[1:]))

