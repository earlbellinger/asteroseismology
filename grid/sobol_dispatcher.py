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
    parser.add_argument('-M', default=[0.7, 1.8], nargs=2, type=float,
                        help='range of masses')
    parser.add_argument('-Z', default=[10**-4, 0.06], nargs=2, type=float,
                        help='range of metallicity values')
    parser.add_argument('-Y', default=1.405366, nargs=2, type=float,
                        help='delta Y / delta Z law')
    parser.add_argument('-N', default=1000, type=int, 
                        help='number of tracks to generate')
    parser.add_argument('-s', '--skip', default=20000, type=int,
                        help='offset for sobol numbers')
    parser.add_argument('-d', '--directory', default="simulations", type=str,
                        help='offset for sobol numbers')
    parser.add_argument('-p', '--parallel', default=1, 
                        type=int, help='number of CPUs to use')
    parser.add_argument('-L', '--light', default=False, action='store_true',
                        help='only calculate light element diffusion')
    parser.add_argument('-r', '--remove', default=False, action='store_true',
                        help='delete models upon completion')
    parser.add_argument('-m', '--memory', default=1709000, type=int,
                        help='set maximum memory consumption for job')
    parser.add_argument('-l', '--logs', default=[0, 1], 
                        type=list,
                        help='booleans of whether to log M, Z')
    parser.add_argument('-n', '--nice', default=False, action='store_true',
                        help='run as nice job')
    parser.add_argument('-c', '--cpu', default=0, type=str,
                        help='machine to run the job on')
    args = parser.parse_args(arguments)
    print(args)
    ranges = np.vstack((args.M, args.Z))
    for i, val in enumerate(ranges):
        if args.logs[i]:
            ranges[i] = np.log10(ranges[i])
    print(ranges)
    dispatch(ranges=ranges, YoverZ=args.Y,
        N=args.N, logs=args.logs,
        directory=args.directory, light=args.light, remove=args.remove,
        skip=args.skip, parallel=args.parallel, nice=args.nice, 
        memory=args.memory, cpu=args.cpu)

def dispatch(ranges, YoverZ, N, logs, 
             directory, light=0, remove=0, skip=0, 
             parallel=r"$OMP_NUM_THREADS", nice=0, memory=0, cpu=0):
    shift = ranges[:,0]
    scale = np.array([(b-a) for a,b in ranges])
    init_conds = []
    for i in range(skip, N+skip):
        vals = shift+np.array(i4_sobol(len(ranges), i)[0])*scale
        for j, val in enumerate(vals):
            if logs[j]:
                vals[j] = 10**val
        #init_conds += [[tmp for tmp in vals]]
        #if vals[0] >= diffusion_cutoff: vals[-1] = 0 # if M>1.5 then D=0
        M = vals[0]
        Z = vals[1]
        Y = 0.2463 + YoverZ*vals[1]
        vals = [M, Y, Z]
        print(np.log10( Z / (1-Y-Z) / 0.02293 ))
        bash_cmd = "maybe_sub.sh %s%s%s-p %d ./dispatch.sh -d %s "\
            "-M %.6f -Y %.6f -Z %.6f %s%s"%\
            tuple(["-n " if nice else ""] + 
                  ["-m %d "%memory if memory>0 else ""] +
                  ["-c %s "%cpu if cpu!=0 else ""] +
                  [parallel, directory] + 
                  [val for val in vals] + 
                  ["-L " if light else ""] +
                  ["-r " if remove else ""])
        print(bash_cmd)
        #exit()
        process = subprocess.Popen(bash_cmd.split(), shell=False)
        process.wait()
        #sleep(0.2)
    #np.savetxt('initial_conditions.dat', np.array(init_conds))

if __name__ == '__main__':
    import sys
    exit(main(sys.argv[1:]))

