#### Convert a sensible format back to FGONG 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

import numpy as np
import sys
import argparse
import fgong2ascii

def main(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                        help='ASCII file', type=str)
    parser.add_argument('-o', '--output', help='FGONG file', type=str)
    args = parser.parse_args(arguments)
    if args.output is None:
        args.output = args.input[-4:]
    ascii2fgong(args.input, args.output)


def ascii2fgong(filename, output):
    


if __name__ == '__main__':
    exit(main(sys.argv[1:]))

