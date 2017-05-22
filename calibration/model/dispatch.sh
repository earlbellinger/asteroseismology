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
## Takes mass M, helium Y, metallicity Z, mixing length parameter alpha,
##       overshoot o, and diffusion D
## -L means only hydrogen/helium diffusion (light elements)
## -d is the directory where the simulations should be put 
## -r means delete the calculations afterwards and leave only the product 
## -s means suppress frequency calculations 
while [ "$#" -gt 0 ]; do
  case "$1" in
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;
    -o) overshoot="$2"; shift 2;;
    -D) diffusion="$2"; shift 2;;
    -t) age="$2"; shift 2;;
    -d) directory="$2"; shift 2;;
    
    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done

## Set defaults if they weren't supplied
if [ -z ${HELP+x} ]; then HELP=0; fi
if [ -z ${M+x} ]; then M=1; fi
if [ -z ${Y+x} ]; then Y=0.275; fi
if [ -z ${Z+x} ]; then Z=0.018; fi
if [ -z ${alpha+x} ]; then alpha=1.81; fi
if [ -z ${overshoot+x} ]; then overshoot=0.09105378; fi
if [ -z ${diffusion+x} ]; then diffusion=1; fi
if [ -z ${age+x} ]; then age=20000000000; fi
if [ -z ${directory+x} ]; then directory=simulations; fi
if [ -z ${light+x} ]; then light=0; fi
if [ -z ${remove+x} ]; then remove=0; fi
if [ -z ${suppress+x} ]; then suppress=0; fi

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

## Deletes all of the calculations if the remove flag is set 
cleanup() {
    if [ $remove -eq 1 ]; then
        rm -rf "$dirname"
    fi
}

set_params() { 
    change_param 'initial_mass' "$M" "$inlist"
    change_param 'initial_y' "$Y" "$inlist"
    change_param 'initial_z' "$Z" "$inlist"
    change_param 'new_Y' "$Y" "$inlist"
    change_param 'new_Z' "$Z" "$inlist"
    change_param 'Zbase' "$Z" "$inlist"
    change_param 'mixing_length_alpha' "$alpha" "$inlist"
    change_param 'max_age' "$age" "$inlist"

    if (( $(echo "$overshoot > 0" | bc -l) )); then
        change_param 'step_overshoot_f_above_nonburn_core' "$overshoot" "$inlist"
        change_param 'step_overshoot_f_above_nonburn_shell' "$overshoot" "$inlist"
        change_param 'step_overshoot_f_below_nonburn_shell' "$overshoot" "$inlist"
        change_param 'step_overshoot_f_above_burn_h_core' "$overshoot" "$inlist"
        change_param 'step_overshoot_f_above_burn_h_shell' "$overshoot" "$inlist"
        change_param 'step_overshoot_f_below_burn_h_shell' "$overshoot" "$inlist"
        
        f0=$(echo "scale=8; $overshoot / 5" | bc -l)
        change_param 'overshoot_f0_above_nonburn_core' "$f0" "$inlist"
        change_param 'overshoot_f0_above_nonburn_shell' "$f0" "$inlist"
        change_param 'overshoot_f0_below_nonburn_shell' "$f0" "$inlist"
        change_param 'overshoot_f0_above_burn_h_core' "$f0" "$inlist"
        change_param 'overshoot_f0_above_burn_h_shell' "$f0" "$inlist"
        change_param 'overshoot_f0_below_burn_h_shell' "$f0" "$inlist"
    fi
}

set_diffusion() {
    if (( $(echo "$diffusion > 0" | bc -l) )); then
        change_param 'do_element_diffusion' '.true.' "$inlist"
        change_param 'diffusion_class_factor(1)' "$diffusion" "$inlist"
        change_param 'diffusion_class_factor(2)' "$diffusion" "$inlist"
        change_param 'diffusion_class_factor(3)' "$diffusion" "$inlist"
        change_param 'diffusion_class_factor(4)' "$diffusion" "$inlist"
        change_param 'diffusion_class_factor(5)' "$diffusion" "$inlist"
    fi
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
#change_inlists "inlist_2ms"
set_inlist "inlist_2ms"
set_diffusion
./rn

$SCRIPTS_DIR/fgong2freqs.sh LOGS_MS/profile1.data.FGONG 

# fin.

