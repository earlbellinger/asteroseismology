#!/bin/bash
#### Check if the queuing system is available. Otherwise, run locally. 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Updated: March 20, 2018

maybe_sub() {
    cmd=$*
    
    ## Check if job should be run in parallel
    environment="environment  = OMP_NUM_THREADS=$OMP_NUM_THREADS;PYTHONUNBUFFERED=1
request_cpus = $OMP_NUM_THREADS
"
    #threads=""
    #if [ "$OMP_NUM_THREADS" -gt 1 ]; then
    #threads=""
    #fi
    #environment="$environment;$threads"
    
    ## Set memory consumption
    image_size=""
    if [ "$MEMORY" -gt 0 ]; then
        image_size="image_size   = $MEMORY
"
    fi
    
    ## Check niceness control
    nice=""
    if [ $NICE -gt 0 ]; then
        nice="nice_user    = True
"
    fi
    
    requirements=""
    if [ $EXCLUDE -gt 0 ]; then
        requirements="Requirements = (Machine != \"seismo18.mps.mpg.de\" "\
"&& Machine != \"seismo19.mps.mpg.de\")
"
    fi
    
    ## Set machine
    machine=""
    if [ ! -z ${MACHINE+x} ]; then
        machine="Requirements = Machine==\"$MACHINE\"
"
    fi
    
    #name=${cmd// /_}
    name=${cmd//[-. \/]/_}
    #name=${name//[-._\/]/}
    #name=${name//./_}
    if command -v condor_submit >/dev/null 2>&1
      then
        ## Make directory for script
        dname="maybe_sub_logs/$name"
        mkdir -p $dname
        cd $dname
        
        ## Create a shell script
        echo "#!/usr/bin/sh
cd ../..
$cmd
" > "$name.sh"
        chmod +x "$name.sh"
        
        ## Create and submit a condor job 
	    echo "Universe     = vanilla
getenv       = True
Executable   = $name.sh
Output       = condor.out
Error        = condor.error
Log          = condor.log
$environment$image_size$nice$machine$requirements
queue
" > "condor.job"
        condor_submit "condor.job"
        #condor_wait -status -wait 36000 "condor-$name.log"
        cd -
      else
        eval "$*"
    fi
}

## Parse command line arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -h) HELP=1; break;;
    -n) NICE=1; shift 1;;
    -e) EXCLUDE=1; shift 1;;
    -p) OMP_NUM_THREADS="$2"; shift 2;;
    -m) MEMORY="$2"; shift 2;;
    -c) MACHINE="$2"; shift 2;;

     *) break;;
  esac
done

if [ -z ${HELP+x} ]; then HELP=0; fi
if [ -z ${OMP_NUM_THREADS+x} ]; then OMP_NUM_THREADS=1; fi
if [ -z ${NICE+x} ]; then NICE=0; fi
if [ -z ${MEMORY+x} ]; then MEMORY=0; fi
if [ -z ${EXCLUDE+x} ]; then EXCLUDE=0; fi

if [ $HELP -gt 0 ]; then
    echo "                        _                      _     ";
    echo "  _ __ ___   __ _ _   _| |__   ___   ___ _   _| |__  ";
    echo " | '_ \` _ \ / _\` | | | | '_ \ / _ \ / __| | | | '_ \ ";
    echo " | | | | | | (_| | |_| | |_) |  __/ \__ | |_| | |_) |";
    echo " |_| |_| |_|\__,_|\__, |_.__/ \___| |___/\__,_|_.__/ ";
    echo "                  |___/                              ";
    echo ""
    echo "  maybe_sub.sh: automatic job queuer for the condor"
    echo "                queueing system"
    echo
    echo "example:"
    echo "  maybe_sub.sh -n -p 16 echo hi"
    echo
    echo "flags:"
    echo "  -h   : show this helpful message and quit"
    echo "  -n   : run as nice job"
    echo "  -e   : exclude two cores (a.k.a., a 'really' nice job)"
    echo "  -p # : set the number of threads to run in parallel"
    echo "  -m # : set the memory consumption of the job"
    echo "  -c s : set the machine that this job will run on"
    echo
    exit
fi

maybe_sub $*

