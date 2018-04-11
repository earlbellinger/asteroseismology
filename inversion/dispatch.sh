#!/bin/bash
#### Script for dispatching a stellar evolutionary track simulation 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

################################################################################
### GLOBAL VARIABLES ###########################################################
################################################################################
export MESA_DIR=/scratch/seismo/bellinger/MESA/mesa-r8118

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
inlist='inlist_1relax'

################################################################################
### PARSE COMMAND LINE ARGUMENTS ###############################################
################################################################################
M=1
Y=0.273766941652327
Z=0.0198435943135059
alpha=1.88474833300517
overshoot=0
f0=0
age=4570000000
calibrate=1
directory=simulations
save=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;
    -o) overshoot="$2"; shift 2;;
    -f) f0="$2"; shift 2;;
    -D) diffusion="$2"; shift 2;;
    -d) directory="$2"; shift 2;;
    -t) age="$2"; shift 2;;
    -s) save="$2"; shift 2;;
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

set_params() {
    change 'initial_mass' "$M" "$inlist"
    change 'initial_y' "$Y" "$inlist"
    change 'initial_z' "$Z" "$inlist"
    change 'new_Y' "$Y" "$inlist"
    change 'new_Z' "$Z" "$inlist"
    change 'Zbase' "$Z" "$inlist"
    change 'mixing_length_alpha' "$alpha" "$inlist"
    
    if (( $(echo "$calibrate > 0" | bc -l) )); then
        #change "write_profile_when_terminate" ".true." "$inlist"
        change "profile_interval" "99999" "$inlist"
        change "history_interval" "99999" "$inlist"
        change "max_num_profile_models" "1" "$inlist"
        change 'max_age' "$age" "$inlist"
    fi
    
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
    
    if (( $(echo "$save > 0" | bc -l) )); then
        change 'write_pulse_info_with_profile' ".true." "$inlist"
    fi
}

set_diffusion() {
    if (( $(echo "$diffusion > 0" | bc -l) )); then
        change 'do_element_diffusion' '.true.' "$inlist"
        #change 'diffusion_class_factor(1)' "$diffusion" "$inlist"
        #change 'diffusion_class_factor(2)' "$diffusion" "$inlist"
        #change 'diffusion_class_factor(3)' "$diffusion" "$inlist"
        #change 'diffusion_class_factor(4)' "$diffusion" "$inlist"
        #change 'diffusion_class_factor(5)' "$diffusion" "$inlist"
    fi
}

fix_mod() { # final_mod
    mod=$1
    # MESA has a bug that makes it print numbers of the form *.***### 
    # (i.e. no letter for the exponent, and all the numbers are asterisks) 
    # This sed command replaces those with zeros. 
    sed -i.bak "s/\*\.\**-[0-9]*/0.0000000000000000D+00/g" $mod
}

################################################################################
### INITIALIZATION AND EVOLUTION ###############################################
################################################################################
## Make directory and copy over simulation files 
expname=M="$M"_Y="$Y"_Z="$Z"_alpha="$alpha"\
_overshoot="$overshoot"_diffusion="$diffusion"
if (( $(echo "$calibrate > 0" | bc -l) )); then
    dirname="$directory"
else
    dirname="$directory"/"$expname"
fi

mkdir -p $dirname
cd $dirname
cp -r $scriptdir/mesa_template/* .
rm -rf LOGS/*

# pre-main sequence
set_params
./rn
fix_mod "pms.mod"

set_inlist "inlist_2pms"
set_params
./rn
fix_mod "zams.mod"

# main sequence
set_inlist "inlist_3ms"
set_params
#set_diffusion
./rn

# fin.

