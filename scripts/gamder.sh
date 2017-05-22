#!/bin/bash

#### Takes ASCII FGONG and calculates Gamma_1 derivatives using OPAL 2005 EOS 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### Parse command line tokens 
## -h and --help will display help instructions 
if [ -z ${1+x} ] || [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "Calculates Gamma_1 derivatives using OPAL 2005 EOS"
    echo "Usage: gamder.sh ASCIIFGONG"
    exit
fi

## Check that the first input (ASCII FGONG file) exists
if [ ! -e "$1" ]; then
    echo "Error: Cannot locate ASCII FGONG file $1"
    exit 1
fi
ASCII_FGONG=$1

## pull out the columns needed for gamder, i.e. 
# xr, t, p, rho, x, z, ggam1
R --slave -q -e "write.table(read.table('$ASCII_FGONG', header=1)[c('x', 'T', 'P', 'rho', 'X', 'Z', 'Gamma1')], 'gamder.in', row.names=F, col.names=F, quote=F)" 

# find metallicity in convective envelope and use that to pick an EOS table
Z=$(R --slave -q -e "Z <- read.table('$ASCII_FGONG', header=1, nrow=1)[['Z']]; DF <- read.table('$SCRIPTS_DIR/gamder/OPAL2005/eos.table', stringsAsFactors=F); cat(DF[,1][which.min(abs(Z-DF[,2]))])")

# pick core metallicity instead of env
#Z=$(R --slave -q -e "Z <- read.table('$ASCII_FGONG', header=1)[['Z']]; Z <- Z[length(Z)]; DF <- read.table('$SCRIPTS_DIR/gamder/OPAL2005/eos.table', stringsAsFactors=F); cat(DF[,1][which.min(abs(Z-DF[,2]))])")


ln -sf $SCRIPTS_DIR/gamder/OPAL2005/$Z eosfile #EOS5_02z8x

# call gamder
$SCRIPTS_DIR/gamder/gamder.x

