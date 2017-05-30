#!/bin/bash

for mass in 1.064 1.08 1.096; do
    for logR in 0.07918125 0.08635983 0.09342169; do
        maybe_sub.sh -p 4 Rscript calibrate2.R $mass $logR 6900000000 0.1931246 0.096 
    done
done


for mass in 1.015 1.03 1.045; do
    for logR in 0.04139269 0.04921802 0.05690485; do
        maybe_sub.sh -p 4 Rscript calibrate2.R $mass $logR 6800000000 0.1038037 0.052 
    done
done

