#!/bin/bash 
#### Make linear and "boundaries" simulations 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

## Run tracks at the boundaries of each parameter 
# e.g. a track with low mass and high mass, with all other vals fixed at solar 
#dir=simulations-boundaries

#maybe_sub.sh -p dispatch.sh -d $dir -M 0.7
#maybe_sub.sh -p dispatch.sh -d $dir -M 1.3

#maybe_sub.sh -p dispatch.sh -d $dir -Y 0.22
#maybe_sub.sh -p dispatch.sh -d $dir -Y 0.34

#maybe_sub.sh -p dispatch.sh -d $dir -Z 0.0001
#maybe_sub.sh -p dispatch.sh -d $dir -Z 0.04

#maybe_sub.sh -p dispatch.sh -d $dir -a 1.5
#maybe_sub.sh -p dispatch.sh -d $dir -a 2.5

## Vary each parameter linearly from min to max 
# keep all other values fixed at the solar value 
for i in $(seq 0.7 0.1 1.3)
  do maybe_sub.sh -p dispatch.sh -d simulations-linear -M $i
done

for i in $(seq 0.22 0.01 0.34)
  do maybe_sub.sh -p dispatch.sh -d simulations-linear -Y $i
done

for i in $(seq 1.5 0.1 2.5)
  do maybe_sub.sh -p dispatch.sh -d simulations-linear -a $i
done

maybe_sub.sh -p dispatch.sh -d simulations-linear -Z 0.0001
maybe_sub.sh -p dispatch.sh -d simulations-linear -Z 0.04
python3 sobol_dispatcher.py -d simulations-linear -N 10 \
    -M 1.0 1.0 -Y 0.018 0.018 -a 1.85 1.85 -s 10000
