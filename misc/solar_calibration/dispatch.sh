#!/bin/bash
#### Script for dispatching a stellar evolutionary track simulation 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

################################################################################
### GLOBAL VARIABLES ###########################################################
################################################################################
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
num_process=128 # try to get at least this many profile files 
init_mesh_delta_coeff=1 # the number in the inlist file 
mesh_delta_coeff=$init_mesh_delta_coeff
ms_mesh_delta=0.4 # the mesh spacing we will use on the main sequence
max_years_for_timestep=-1 # dt_max 
min_years_for_timestep=1000
profile_interval=1 # how often a profile should be output 
inlist='inlist_pms'

### Log files 
pmslog="pms.log" 
mslog="ms.log"
logfile="mesa.log"

### Success and failure conditions 
msuccess="termination code: max_age\|termination code: xa_central_lower_limit"
pmsuccess="termination code: Lnuc_div_L_zams_limit"
meshfail="mesh_plan problem"
convfail="terminated evolution: convergence problems"

################################################################################
### PARSE COMMAND LINE ARGUMENTS ###############################################
################################################################################
## Takes mass M, helium Y, metallicity Z, mixing length parameter alpha,
##       overshoot o, and diffusion D
## -L means only hydrogen/helium diffusion (light elements)
## -d is the directory where the simulations should be put 
## -r means delete the calculations afterwards and leave only the product 
while [ "$#" -gt 0 ]; do
  case "$1" in
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;
    -o) overshoot="$2"; shift 2;;
    -D) diffusion="$2"; shift 2;;
    -d) directory="$2"; shift 2;;
    -L) light=1; shift 1;;
    -r) remove=1; shift 1;;

    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done

## Set defaults if they weren't supplied
if [ -z ${M+x} ]; then M=1; fi
if [ -z ${Y+x} ]; then Y=0.266; fi
if [ -z ${Z+x} ]; then Z=0.018; fi
if [ -z ${alpha+x} ]; then alpha=1.81; fi
if [ -z ${overshoot+x} ]; then overshoot=0.07; fi
if [ -z ${diffusion+x} ]; then diffusion=1; fi
if [ -z ${directory+x} ]; then directory=simulations; fi
if [ -z ${light+x} ]; then light=0; fi
if [ -z ${remove+x} ]; then remove=0; fi

################################################################################
### HELPER FUNCTIONS ###########################################################
################################################################################
## Modifies the inlist. must specify what param you want to change, 
##     the current value it has, and the new value you want it to have
change_inlists() { #newinlist
    newinlist=$1
    echo "Changing to $newinlist"
    sed -i.bak "s/$inlist/$newinlist/g" "inlist"
    inlist=$newinlist
    set_params
}

change() { #param initval newval
    param=$1
    initval=$2
    newval=$3
    if grep -q "\!$param = $initval" "$inlist"; then
        sed -i.bak "s/\!$param = $initval/$param = $initval/g" "$inlist"
    fi
    sed -i.bak "s/$param = $initval/$param = $newval/g" "$inlist"
}

set_params() {
    change 'initial_mass' '1.0' "$M"
    change 'initial_y' '-1' "$Y"
    change 'initial_z' '0.02' "$Z"
    change "new_Y" '-1' "$Y"
    change "new_Z" '-1' "$Z"
    change "Zbase" '0.02' "$Z"
    change 'mixing_length_alpha' '2.1' "$alpha"

    if (( $(echo "$overshoot > 0" | bc -l) )); then
        change 'step_overshoot_f_above_nonburn_core' '0.005' "$overshoot"
        change 'step_overshoot_f_above_nonburn_shell' '0.005' "$overshoot"
        change 'step_overshoot_f_below_nonburn_shell' '0.005' "$overshoot"
        change 'step_overshoot_f_above_burn_h_core' '0.005' "$overshoot"
        change 'step_overshoot_f_above_burn_h_shell' '0.005' "$overshoot"
        change 'step_overshoot_f_below_burn_h_shell' '0.005' "$overshoot"
        
        f0=$(echo "scale=8; $overshoot / 5" | bc -l)
        change 'overshoot_f0_above_nonburn_core' '0.001' "$f0"
        change 'overshoot_f0_above_nonburn_shell' '0.001' "$f0"
        change 'overshoot_f0_below_nonburn_shell' '0.001' "$f0"
        change 'overshoot_f0_above_burn_h_core' '0.001' "$f0"
        change 'overshoot_f0_above_burn_h_shell' '0.001' "$f0"
        change 'overshoot_f0_below_burn_h_shell' '0.001' "$f0"
    fi
    
    if [ $diffusion -eq 1 ]; then
        change 'do_element_diffusion' '.false.' '.true.'
    fi
}

fix_mod() { # final_mod
    final_mod=$1
    # MESA r8118 has a bug that makes it print numbers of the form *.***### 
    # (i.e. no letter for the exponent, and all the numbers are asterisks) 
    # This sed command replaces those with zeros. 
    sed -i.bak "s/\*\.\**-[0-9]*/0.0000000000000000D+00/g" $final_mod
}

run() { # logs_dir  final_mod
    logs_dir=$1
    final_mod=$2
    ./rn
    fix_mod $final_mod
}

################################################################################
### INITIALIZATION AND EVOLUTION ###############################################
################################################################################
## Make directory and copy over simulation files 
#expname="M=$M""_""Y=$Y""_""Z=$Z""_""alpha=$alpha"\
#"_""overshoot=$overshoot""_""diffusion=$diffusion"
dirname="$directory" # /$expname"

if [ ! -d "$directory" ]; then
    mkdir -p "$dirname"
    cd "$dirname"
    cp -r $scriptdir/mesa/* .
else
    cd "$dirname"
    cp -r $scriptdir/mesa/inlist* .
fi
#rm -rf LOGS/*
#rm -rf LOGS_MS/*

# pre-main sequence
set_params
./rn
fix_mod "zams.mod"

# main sequence
change_inlists "inlist_ms"
run "LOGS_MS" "tams.mod" "history-ms.data"

#fgong2freqs.sh LOGS_MS/profile1.data.FGONG

# fin.

