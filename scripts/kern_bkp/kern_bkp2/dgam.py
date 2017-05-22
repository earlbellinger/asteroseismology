#### Extracts Gamma_1 derivatives from an FGONG file 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

import numpy as np
import sys

def main(arguments):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',  default='ModelS_5b.FGONG',
                        help='FGONG file (default: ModelS_5b.FGONG)',
                        type=str)
    parser.add_argument('-o', '--output',  default="dgam.txt",
                        help='output file (default: dgam.txt)',
                        type=str)
    args = parser.parse_args(arguments)
    fgong = load_fgong(args.input)
    x = fgong['var'][:,0][::-1] / fgong['glob'][1]
    dgam_rho = fgong['var'][:,25][::-1]
    dgam_p = fgong['var'][:,26][::-1]
    dgam_Y = fgong['var'][:,27][::-1]
    #dgam = fgong['var'][:,25:28][::-1]
    #y = np.append(x.reshape(-1,1), dgam, 1)
    y = np.vstack((x, dgam_p, dgam_rho, dgam_Y)).T
    np.savetxt(args.output, y)


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
            elif s[-4].lower() != 'e':
                s = s[:-4] + 'e' + s[-4:]
            
            tmp.append(float(s))
    
    fgong['glob'] = np.array(tmp[:fgong['iconst']])
    fgong['var'] = np.array(tmp[fgong['iconst']:]).reshape((-1,fgong['ivars']))
    
    f.close()
    
    return fgong


if __name__ == '__main__':
    exit(main(sys.argv[1:]))

