#!/bin/bash

for mass in 1.064 1.08 1.096; do
    for logR in 0.07918125 0.08635983 0.09342169; do
        maybe_sub.sh -p 4 Rscript calibrate2.R $mass $logR 6900000000 0.1931246 0.096 
    done
done

