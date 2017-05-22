#!/bin/bash

for mass in 0.984 1 1.016; do
    for radius in 0.98 1 1.02; do
        maybe_sub.sh -p 2 Rscript vary_RM.R 0 $mass $radius
    done
done

for subdir in $( ls models ); do
    cd models/$subdir/LOGS_MS
    for ii in `seq 1 6`; do
        maybe_sub.sh kerexact.sh profile1.data.FGONG profile1-freqs $ii
    done
    cd -
done

