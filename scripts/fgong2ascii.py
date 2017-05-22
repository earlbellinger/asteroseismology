#### Convert FGONG file to a sensible format 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

import numpy as np
import sys
import argparse

add_vars = 5
var_names = ('x', 'm', 'c2', 'u', 'Y', 
             'r', 'q', 'T', 'P', 'rho', 'X', 'L_r', 'kappa', 
             'eps', 'Gamma1', 'nabla_ad', 'delta', 'c_p', 'mu_e_inv', 
             'conv_stab', 'r_X', 'Z', 'R_min_r', 'eps_g', 'L_g', 
             'X_He3', 'X_C12', 'X_C13', 'X_N14', 'X_O16', 
             'Gamma1_rho', 'Gamma1_P', 'Gamma1_Y', 
             'X_H2', 'X_He4', 'X_Li7', 'X_Be7', 
             'X_N15', 'X_O17', 'X_O18', 'X_Ne20', 
             'blank1', 'blank2', 'blank3', 'blank4') 

def main(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                        help='FGONG file', type=str)
    parser.add_argument('-o', '--output', help='ASCII file', type=str)
    args = parser.parse_args(arguments)
    if args.output is None:
        args.output = args.input + '.dat'
    fgong2ascii(args.input, args.output)


def fgong2ascii(filename, output):
    fgong = load_fgong(filename)
    u = fgong['var'][:,3] / fgong['var'][:,4]
    np.savetxt(output, 
        np.vstack(( fgong['var'][:,0] / fgong['glob'][1], # x
                10**fgong['var'][:,1] * fgong['glob'][0], # m
                    fgong['var'][:,9] * u, # c2 = Gamma_1 * P / rho
                    u, # u = P / rho
                1 - fgong['var'][:,5] + fgong['var'][:,16], # Y = 1-X-Z
                    fgong['var'].T)).T, # rest of FGONG
        delimiter='\t', comments='',
        header='\t'.join((var_names[:add_vars+fgong['var'].shape[1]])))


def load_fgong(filename, N=16):
    '''Given an FGONG file, returns a Python dictionary containing
    NumPy arrays that correspond to the structures in the
    specification of the FGONG format:
    
    https://www.astro.up.pt/corot/ntools/docs/CoRoT_ESTA_Files.pdf
    
    That is, the dictionary has arrays indexed by 'nn', 'iconst',
    'ivar', 'ivers', 'glob' and 'var'.'''
    f = open(filename, 'r')
    
    fgong = {'header':[]}
    for i in range(4):
        fgong['header'].append(f.readline())
    
    tmp = [int(i) for i in f.readline().split()]
    fgong['nn'] = tmp[0]
    fgong['iconst'] = tmp[1]
    fgong['ivars'] = tmp[2]
    fgong['ivers'] = tmp[3]
    
    lines = f.readlines()
    tmp = []
    for line in lines:
        for i in range(len(line)//N):
            s = line[i*N:i*N+N]
            if s[-9:] == '-Infinity':
                s = '-Inf'
            elif s[-9:] == ' Infinity':
                s = 'Inf'
            elif s[-3:] == 'NaN':
                s = np.nan
            elif s[-4].lower() != 'e':
                s = s[:-4] + 'e' + s[-4:]
            
            tmp.append(float(s))
    
    fgong['glob'] = np.array(tmp[:fgong['iconst']])
    fgong['var'] = np.array(tmp[fgong['iconst']:]).reshape((-1,fgong['ivars']))
    
    f.close()
    
    return fgong


if __name__ == '__main__':
    exit(main(sys.argv[1:]))

