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
    init = parser.add_argument_group('initial conditions')
    init.add_argument('-M', default=[0.70, 0.84], nargs=2, type=float,
                      help='range of masses (solar = 1)')
    init.add_argument('-Y', default=[0.22, 0.30], nargs=2, type=float, 
                      help='range of helium values (solar = 0.273)')
    init.add_argument('-Z', default=[0.001, 0.02], nargs=2, type=float,
                      help='range of metallicity values (solar = 0.0185)')
    init.add_argument('-a', '--alpha', default=[0.5, 2.5], nargs=2, type=float,
                      help='range of mixing length parameter values (solar = 1.83)')
    init.add_argument('-t', '--age', default=[5, 13.799], nargs=2, type=float,
                      help='range of ages')
    init.add_argument('-b', '--beta', default=[-0.1, 0.1], nargs=2, type=float,
                      help='range of gravitational power law parameters')
    init.add_argument('-l', '--logs', default=[0, 0, 1, 0, 0, 0], nargs=8,
                      type=list, 
                      help='booleans of whether to log M, Y, Z, alpha, age, beta')
    
    job = parser.add_argument_group('job')
    job.add_argument('-N', '--tracks', default=65536, type=int, 
                     help='number of tracks to generate')
    job.add_argument('-s', '--skip', default=20000, type=int,
                     help='offset for sobol numbers')
    job.add_argument('-d', '--directory', default="simulations", type=str,
                     help='offset for sobol numbers')
    job.add_argument('-p', '--parallel', default=1, 
                     type=int, help='number of CPUs to use')
    job.add_argument('-m', '--image', default=1709000, type=int,
                     help='set default image size for job')
    job.add_argument('-r', '--remove', default=False, action='store_true',
                     help='delete models upon completion')
    job.add_argument('-n', '--nice', default=False, action='store_true',
                     help='run as nice job')
    job.add_argument('-x', '--exclude', default=False, action='store_true',
                     help='exclude compute cores for extra niceness')
    
    args = parser.parse_args(arguments)
    print(args)
    
    ranges = np.vstack((args.M, args.Y, args.Z, args.alpha, 
        args.age, args.beta))
    
    for i, val in enumerate(ranges):
        if args.logs[i]:
            ranges[i] = np.log10(ranges[i])
    print(ranges)
    
    dispatch(ranges=ranges, tracks=args.tracks, 
        logs=args.logs, directory=args.directory, remove=args.remove,
        skip=args.skip, parallel=args.parallel, nice=args.nice, 
        image=args.image)

def dispatch(ranges, tracks, logs, directory, remove=0, skip=0, 
             parallel=r"$OMP_NUM_THREADS", nice=0, image=0):
    shift = ranges[:,0]
    scale = np.array([(b-a) for a,b in ranges])
    init_conds = []
    for i in range(skip, tracks+skip):
        vals = shift+np.array(i4_sobol(len(ranges), i)[0])*scale
        for j, val in enumerate(vals):
            if logs[j]:
                vals[j] = 10**val if not np.isnan(val) else 0
        init_conds += [[tmp for tmp in vals]]
        for j, val in enumerate(vals):
            if np.isnan(vals[j]):
                vals[j] = 0
        
        bash_cmd = "maybe_sub.sh -e %s%s-p %d ./dispatch.sh -d %s "\
            "-n %d "\
            "-M %.6f -Y %.6f -Z %.6f -a %.6f -t %.6fd9 -b %.6f "\
            "%s"%\
            tuple(["-n " if nice else ""] + 
                  ["-i %d "%image if image>0 else "" ] +
                  [parallel, directory] + 
                  [i] + 
                  [val for val in vals] + 
                  ["-r " if remove else ""])
        print(bash_cmd)
        #exit()
        process = subprocess.Popen(bash_cmd.split(), shell=False)
        process.wait()
    #np.savetxt('initial_conditions.dat', np.array(init_conds))

if __name__ == '__main__':
    import sys
    exit(main(sys.argv[1:]))

