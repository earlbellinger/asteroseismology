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
    init.add_argument('-M', default=[0.7, 3.0], nargs=2, type=float,
                      help='range of masses (solar = 1)')
    init.add_argument('-Y', default=[0.22, 0.34], nargs=2, type=float, 
                      help='range of helium values (solar = 0.273)')
    init.add_argument('-Z', default=[10**-4, 0.04], nargs=2, type=float,
                      help='range of metallicity values (solar = 0.0185)')
    init.add_argument('-a', '--alpha', default=[1, 3], nargs=2,type=float,
                      help='range of mixing length parameter values (solar = 1.83)')
    init.add_argument('-o', '--overshoot', default=[10**-4, 1], nargs=2,
                      type=float, help='range of step overshoot')
    init.add_argument('-oe', '--over_exp', default=[10**-4, 1], nargs=2,
                      type=float, help='range of exponential overshoot')
    init.add_argument('-u', '--undershoot', default=[10**-4, 1], nargs=2,
                      type=float, help='range of step undershoot')
    init.add_argument('-ue', '--under_exp', default=[10**-4, 3], nargs=2,
                      type=float, help='range of exponential undershoot')
    init.add_argument('-D', '--diffusion', default=[10**-4, 3], nargs=2,
                      type=float, help='range of diffusion factors')
    init.add_argument('-g', '--grav_sett', default=[10**-4, 3], nargs=2,
                      type=float, help='range of gravitational settling')
    init.add_argument('-e', '--eta', default=[10**-4, 2], nargs=2,
                      type=float, help="range of Reimer's mass loss parameter")
    init.add_argument('-l', '--logs', default=[0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1], nargs=8,
                      type=list, 
                      help='booleans of whether to log M, Y, Z, alpha, o, oe, u, ue, D, g, eta')
    init.add_argument('-T', '--threshold', 
                      type=list, help='consider as 0 if <= this value', nargs=8,
                      default=[0, 0, 0, 0, 10**-3, 10**-3, 10**-3, 10**-3, 10**-3, 10**-3, 10**-3])
    
    job = parser.add_argument_group('job')
    job.add_argument('-N', '--tracks', default=65536, type=int, 
                     help='number of tracks to generate')
    job.add_argument('-i', '--points', default=128, type=int, 
                     help='number of models per phase of track')
    job.add_argument('-s', '--skip', default=20000, type=int,
                     help='offset for sobol numbers')
    job.add_argument('-d', '--directory', default="simulations", type=str,
                     help='offset for sobol numbers')
    job.add_argument('-p', '--parallel', default=1, 
                     type=int, help='number of CPUs to use')
    job.add_argument('-m', '--memory', default=1709000, type=int,
                     help='set maximum memory consumption for job')
    job.add_argument('-r', '--remove', default=False, action='store_true',
                     help='delete models upon completion')
    job.add_argument('-n', '--nice', default=False, action='store_true',
                     help='run as nice job')
    job.add_argument('-rotk', '--rotk', default=False, action='store_true',
                     help='calculate rotation kernels')
    
    physics = parser.add_argument_group('physics')
    physics.add_argument('-L', '--light', default=False, action='store_true',
                         help='only calculate light element diffusion')
    physics.add_argument('-MS', '--mainseq', default=False, action='store_true',
                         help='stop after main sequence')
    physics.add_argument('-S', '--subgiant', default=False, action='store_true',
                         help='stop after subgiant phase')
    physics.add_argument('-t', '--taper', default=False, action='store_true',
                         help='turn down diffusion with increasing mass')
    physics.add_argument('-c', '--chem_ev', default=False, action='store_true',
        help='set Y = c*Z + 0.2463, where the -Y flag becomes the c range (solar c = 1.4276221707417)')
    physics.add_argument('-C', '--couple', default=False, action='store_true',
        help='couple diffusion to gravitational settling')
    
    args = parser.parse_args(arguments)
    print(args)
    
    ranges = np.vstack((args.M, args.Y, args.Z, args.alpha,
        args.overshoot, args.over_exp, args.undershoot, args.under_exp,
        args.diffusion, args.grav_sett, args.eta))
    
    for i, val in enumerate(ranges):
        if args.logs[i]:
            ranges[i] = np.log10(ranges[i])
    print(ranges)
    
    dispatch(ranges=ranges, tracks=args.tracks, points=args.points, 
        logs=args.logs, threshold=args.threshold,
        directory=args.directory, light=args.light, remove=args.remove,
        skip=args.skip, parallel=args.parallel, nice=args.nice, 
        memory=args.memory, mainseq=args.mainseq, subgiant=args.subgiant, 
        taper=args.taper, chem_ev=args.chem_ev, rotk=args.rotk, 
        couple=args.couple)

def dispatch(ranges, tracks, points, logs, threshold, directory, 
             light=0, remove=0, skip=0, 
             parallel=r"$OMP_NUM_THREADS", nice=0, memory=0, 
             mainseq=0, subgiant=0, taper=0, chem_ev=0, rotk=0, couple=0):
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
            if vals[j] <= threshold[j] or np.isnan(vals[j]):
                vals[j] = 0
        if chem_ev:
            vals[1] = vals[1] * vals[2] + 0.2463
        if couple:
            vals[9] = vals[8]
        
        bash_cmd = "maybe_sub.sh -e %s%s-p %d ./dispatch.sh -d %s "\
            "-n %d -N %d "\
            "-M %.6f -Y %.6f -Z %.6f -a %.6f "\
            "-o %.6f -oe %.6f -u %.6f -ue %.6f "\
            "-D %.6f -g %.6f -e %.6f "\
            "%s%s%s%s%s%s"%\
            tuple(["-n " if nice else ""] + 
                  ["-m %d "%memory if memory>0 else "" ] +
                  [parallel, directory] + 
                  [i] + [points] +
                  [val for val in vals] + 
                  ["-L " if light else ""] +
                  ["-r " if remove else ""] +
                  ["-MS " if mainseq else ""] +
                  ["-S " if subgiant else ""] +
                  ["-t " if taper else ""] +
                  ["-rotk " if rotk else ""])
        print(bash_cmd)
        #exit()
        process = subprocess.Popen(bash_cmd.split(), shell=False)
        process.wait()
    #np.savetxt('initial_conditions.dat', np.array(init_conds))

if __name__ == '__main__':
    import sys
    exit(main(sys.argv[1:]))

