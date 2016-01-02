#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

dualexp() {
    expname="M=$M""_""Y=$Y""_""Z=$Z""_""alpha=$alpha"
    dirname="$directory/$expname"
    
    mkdir -p "$dirname"
    cp -r $scriptdir/mesa_template/*\
        "$dirname"
    cd "$dirname"
    
    change 'initial_mass' '1.0' "$M"
    change 'initial_y' '-1' "$Y"
    change 'initial_z' '0.02' "$Z"
    change 'Zbase' '0.02' "$Z"
    change "new_Y" '-1' "$Y"
    change "new_Z" '-1' "$Z"
    change 'mixing_length_alpha' '2.1' "$alpha"
    
    mk
    #run "ZAMS"
    rn
    mv LOGS/history.data .
    
    change "stop_near_zams" ".true." ".false."
    change "create_pre_main_sequence_model" ".true." ".false."
    change "load_saved_model" ".false." ".true."
    change "write_profiles_flag" ".false." ".true."
    
    change "do_element_diffusion" ".false." ".true."
    
    #change "change_Y" ".false." ".true."
    #change "change_Z" ".false." ".true."
    
    rerun #"Hexh"
    mv history.data LOGS
    
    find "LOGS" -maxdepth 1 -type f -name "*.FGONG" | xargs -i \
        --max-procs=$OMP_NUM_THREADS bash -c \
        "echo start {}; fgong2freqs.sh {}; echo end {}"
    
    cd ..
    Rscript $scriptdir/seismology.R "$expname"
    #rm -rf "$expname"
}

change() { #param initval newval
    sed -i.bak "s/\!$1 = $2/$1 = $3/g" inlist_1.0
    sed -i.bak "s/$1 = $2/$1 = $3/g" inlist_1.0
}

#run() { # nameOfRun
#    #condor_submit mesa.job
#    #condor_wait condor.log
#    rn
#    mv final_profile.data "LOGS/profile_$1.data"
#    mv final_profile.data.FGONG "LOGS/profile_$1.data.FGONG"
#}

rerun() { # nameOfRun
    #run $1
    rn
    cp history.data "history.$1.data"
    tail -n+7 LOGS/history.data >> history.data
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;
    -d) directory="$2"; shift 2;;

    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done

if [ -z ${M+x} ]; then M=1; fi
if [ -z ${Y+x} ]; then Y=0.275; fi
if [ -z ${Z+x} ]; then Z=0.018; fi
if [ -z ${alpha+x} ]; then alpha=1.85; fi
if [ -z ${directory+x} ]; then directory=simulations; fi

dualexp

