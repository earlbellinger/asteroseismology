#!/bin/bash
#### Script for dispatching a stellar evolutionary track simulation 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

################################################################################
### GLOBAL VARIABLES ###########################################################
################################################################################
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
inlist='inlist_1pms'

################################################################################
### PARSE COMMAND LINE ARGUMENTS ###############################################
################################################################################
## Takes mass M, helium Y, metallicity Z, and mixing length parameter alpha
## -d is the directory where the simulations should be put 
while [ "$#" -gt 0 ]; do
  case "$1" in
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;
    -d) directory="$2"; shift 2;;
    -m) mesh_delta_coeff="$2"; shift 2;;
    -v) varcontrol_target="$2"; shift 2;;
    
    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done

## Set defaults if they weren't supplied
if [ -z ${M+x} ]; then M=1; fi
if [ -z ${Y+x} ]; then Y=0.275; fi
if [ -z ${Z+x} ]; then Z=0.018; fi
if [ -z ${alpha+x} ]; then alpha=1.81; fi
if [ -z ${directory+x} ]; then directory=simulations; fi
if [ -z ${mesh_delta_coeff+x} ]; then mesh_delta_coeff=1; fi
if [ -z ${varcontrol_target+x} ]; then varcontrol_target='1d-4'; fi

################################################################################
### HELPER FUNCTIONS ###########################################################
################################################################################
change_param() { 
    # Modifies a parameter in the current inlist. 
    # args: ($1) name of parameter 
    #       ($2) new value 
    #       ($3) filename of inlist where change should occur 
    # Additionally changes the 'inlist_0all' inlist. 
    # example command: change_param initial_mass 1.3 
    # example command: change_param log_directory 'LOGS_MS' 
    # example command: change_param do_element_diffusion .true. 
    param=$1 
    newval=$2 
    filename=$3 
    escapedParam=$(sed 's/[^^]/[&]/g; s/\^/\\^/g' <<< "$param")
    search="^\s*\!*\s*$escapedParam\s*=.+$" 
    replace="      $param = $newval" 
    sed -r -i.bak -e "s/$search/$replace/g" $filename 
    if [ ! "$filename" == 'inlist_0all' ]; then 
        change_param $param $newval 'inlist_0all'
    fi 
} 

set_inlist() { 
    # Changes to a different inlist by modifying where "inlist" file points 
    # args: ($1) filename of new inlist  
    # example command: change_inlists inlist_2ms 
    newinlist=$1 
    echo "Changing to $newinlist" 
    change_param "extra_star_job_inlist2_name" "'$newinlist'" "inlist" 
    change_param "extra_controls_inlist2_name" "'$newinlist'" "inlist" 
    inlist=$newinlist
}

set_params() { 
    change_param 'initial_mass' "$M" "$inlist"
    change_param 'initial_y' "$Y" "$inlist"
    change_param 'initial_z' "$Z" "$inlist"
    change_param 'new_Y' "$Y" "$inlist"
    change_param 'new_Z' "$Z" "$inlist"
    change_param 'Zbase' "$Z" "$inlist"
    change_param 'mixing_length_alpha' "$alpha" "$inlist"
    change_param 'mesh_delta_coeff' "$mesh_delta_coeff" "$inlist"
    change_param 'varcontrol_target' "$varcontrol_target" "$inlist"
}

fix_mod() { # final_mod
    final_mod=$1
    # MESA has a bug that makes it print numbers of the form *.***### 
    # (i.e. no letter for the exponent, and all the numbers are asterisks) 
    # This sed command replaces those with zeros. 
    sed -i.bak "s/\*\.\**-[0-9]*/0.0000000000000000D+00/g" $final_mod
}

################################################################################
### INITIALIZATION AND EVOLUTION ###############################################
################################################################################
## Make directory and copy over simulation files 
dirname="$directory"

if [ ! -d "$directory" ]; then
    mkdir -p "$dirname"
    cd "$dirname"
    cp -r $scriptdir/mesa/* .
else
    cd "$dirname"
    cp -r $scriptdir/mesa/inlist* .
fi

# pre-main sequence
set_params
./rn
fix_mod "zams.mod"

# main sequence
set_inlist "inlist_2ms"
set_params
./rn

# fin.

