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
inlist='inlist_1relax'
continue=1

################################################################################
### PARSE COMMAND LINE ARGUMENTS ###############################################
################################################################################
## Takes mass M, helium Y, metallicity Z, mixing length parameter alpha,
##       overshoot o, and diffusion D
## -L means only hydrogen/helium diffusion (light elements)
## -d is the directory where the simulations should be put 
## -r means delete the calculations afterwards and leave only the product 
## -s means suppress frequency calculations 
M=1
Y=0.27202387
Z=0.01830403
alpha=1.84663590
overshoot=0.09104194
age=4572000000
if (( $(echo "$M <= 1.2" | bc -l ) )); then 
    diffusion=1
elif (( $(echo "$M >= 1.3" | bc -l ) )); then 
    diffusion=0
else
    x=$(echo "scale=10; 1 - ( 1.3 - $M )*10" | bc -l)
    diffusion=$(sigmoid $x)
fi
directory=simulations
save=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;
    -t) age="$2"; shift 2;;
    -o) overshoot="$2"; shift 2;;
    -D) diffusion="$2"; shift 2;;
    -d) directory="$2"; shift 2;;
    -s) save=1; shift 1;;
    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done

################################################################################
### HELPER FUNCTIONS ###########################################################
################################################################################
change() { 
    # Modifies a parameter in the current inlist. 
    # args: ($1) name of parameter 
    #       ($2) new value 
    #       ($3) filename of inlist where change should occur 
    # Additionally changes the 'inlist_0all' inlist. 
    # example command: change initial_mass 1.3 
    # example command: change log_directory 'LOGS_MS' 
    # example command: change do_element_diffusion .true. 
    param=$1 
    newval=$2 
    filename=$3 
    escapedParam=$(sed 's/[^^]/[&]/g; s/\^/\\^/g' <<< "$param")
    search="^\s*\!*\s*$escapedParam\s*=.+$" 
    replace="      $param = $newval" 
    sed -r -i.bak -e "s/$search/$replace/g" $filename 
    if [ ! "$filename" == 'inlist_0all' ]; then 
        change $param $newval 'inlist_0all'
    fi 
} 

set_inlist() { 
    # Changes to a different inlist by modifying where "inlist" file points 
    # args: ($1) filename of new inlist  
    # example command: change_inlists inlist_2ms 
    newinlist=$1 
    echo "Changing to $newinlist" 
    change "extra_star_job_inlist2_name" "'$newinlist'" "inlist" 
    change "extra_controls_inlist2_name" "'$newinlist'" "inlist" 
    inlist=$newinlist
}

## Deletes all of the calculations if the remove flag is set 
cleanup() {
    if [ $remove -eq 1 ]; then
        rm -rf "$dirname"
    fi
}

set_params() {
    change 'initial_mass' "$M" "$inlist"
    change 'initial_y' "$Y" "$inlist"
    change 'initial_z' "$Z" "$inlist"
    change 'new_Y' "$Y" "$inlist"
    change 'new_Z' "$Z" "$inlist"
    change 'Zbase' "$Z" "$inlist"
    change 'max_age' "$age" "$inlist"
    change 'mixing_length_alpha' "$alpha" "$inlist"
}

set_diffusion() {
    if (( $(echo "$diffusion > 0" | bc -l) )); then
        change 'do_element_diffusion' '.true.' "$inlist"
        change 'diffusion_class_factor(1)' "$diffusion" "$inlist"
        change 'diffusion_class_factor(2)' "$diffusion" "$inlist"
        change 'diffusion_class_factor(3)' "$diffusion" "$inlist"
        change 'diffusion_class_factor(4)' "$diffusion" "$inlist"
        change 'diffusion_class_factor(5)' "$diffusion" "$inlist"
        if [[ taper -eq 1 ]]; then
            decrease=$(echo "scale=8; $diffusion / 100" | bc -l)
            change 'x_ctrl(2)' "$decrease" "$inlist"
        fi
    fi
}

set_overshoot() {
    if (( $(echo "$overshoot > 0" | bc -l) )); then
        change 'step_overshoot_f_above_nonburn_core' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_nonburn_shell' "$overshoot" "$inlist"
        change 'step_overshoot_f_below_nonburn_shell' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_burn_h_core' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_burn_h_shell' "$overshoot" "$inlist"
        change 'step_overshoot_f_below_burn_h_shell' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_burn_he_core' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_burn_he_shell' "$overshoot" "$inlist"
        change 'step_overshoot_f_below_burn_he_shell' "$overshoot" "$inlist"
        
        f0=$(echo "scale=8; $overshoot / 5" | bc -l)
        change 'overshoot_f0_above_nonburn_core' "$f0" "$inlist"
        change 'overshoot_f0_above_nonburn_shell' "$f0" "$inlist"
        change 'overshoot_f0_below_nonburn_shell' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_h_core' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_h_shell' "$f0" "$inlist"
        change 'overshoot_f0_below_burn_h_shell' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_he_core' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_he_shell' "$f0" "$inlist"
        change 'overshoot_f0_below_burn_he_shell' "$f0" "$inlist"
    fi
}

fix_mod() { # final_mod
    mod=$1
    check_mod $mod
    # MESA has a bug that makes it print numbers of the form *.***### 
    # (i.e. no letter for the exponent, and all the numbers are asterisks) 
    # This sed command replaces those with zeros. 
    if [[ continue -eq 1 ]]; then
        sed -i.bak "s/\*\.\**-[0-9]*/0.0000000000000000D+00/g" $mod
    fi
}

check_mod() {
    mod=$1
    # check that mod file got written
    if [[ ! -e $mod ]]; then 
        continue=0
    fi 
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
    cp -r $scriptdir/mesa/* . #inlist* .
fi

rm -rf LOGS LOGS_MS

# pre-main sequence
set_params
./rn
fix_mod "pms.mod"

if [[ continue -eq 0 ]]; then
    rm -rf zams.mod pms.mod
    exit 1
fi 

set_inlist "inlist_2pms"
set_diffusion
set_overshoot
./rn
fix_mod "zams.mod"

if [[ continue -eq 0 ]]; then
    rm -rf zams.mod pms.mod 
    exit 1
fi 

set_inlist "inlist_3ms"
set_params
if [ $save -eq 1 ]; then
    change 'write_profile_when_terminate' '.true.' "$inlist"
    change 'write_profiles_flag' '.true.' "$inlist"
    change 'history_interval' '1' "$inlist"
fi
./rn

rm -rf zams.mod pms.mod 

# fin.

