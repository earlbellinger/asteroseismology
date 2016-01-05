#!/bin/bash
#### Check if the queuing system is available. Otherwise, run locally. 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

maybe_sub() {
    cmd=$*
    
    ## Check if it should be run in parallel
    threads="request_cpus = 1"
    if [ "-p" == ${cmd:0:3} ]
      then
        threads="environment  = OMP_NUM_THREADS=$OMP_NUM_THREADS
request_cpus = $OMP_NUM_THREADS"
        cmd=${cmd:3}
    fi
    
    name=${cmd// /_}
    if command -v condor_submit >/dev/null 2>&1
      then
        ## Make directory for script
        dname=maybe_sub_logs
        mkdir -p $dname
        cd $dname
        
        ## Create a shell script
        echo "#!/usr/bin/sh
cd ..
$cmd
" > "condor-$name.sh"
        chmod +x "condor-$name.sh"
        
        ## Create and submit a condor job 
	    echo "Universe     = vanilla
getenv       = True
Executable   = condor-$name.sh
Output       = condor-$name.out
Error        = condor-$name.error
Log          = condor-$name.log
$threads

queue
" > "condor-$name.job"
        condor_submit "condor-$name.job"
        #condor_wait -status -wait 36000 "condor-$name.log"
        cd -
      else
        eval "$*"
    fi
}

if [ -z ${OMP_NUM_THREADS+x} ]; then OMP_NUM_THREADS=1; fi

maybe_sub $*

