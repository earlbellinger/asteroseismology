#!/bin/bash
#### Script for dispatching a stellar evolutionary track simulation 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

simulate() {
    expname="M=$M""_""Y=$Y""_""Z=$Z""_""alpha=$alpha"\
"_""diffusion=$diffusion""_""overshoot=$overshoot"
    dirname="$directory/$expname"
    
    mkdir -p "$dirname"
    cd "$dirname"
    cp -r $scriptdir/mesa_template/* .
    
    change 'initial_mass' '1.0' "$M"
    change 'initial_y' '-1' "$Y"
    change 'initial_z' '0.02' "$Z"
    change 'Zbase' '0.02' "$Z"
    change "new_Y" '-1' "$Y"
    change "new_Z" '-1' "$Z"
    change 'mixing_length_alpha' '2.1' "$alpha"
    
    change 'step_overshoot_f_above_nonburn_core' '0.005' "$overshoot"
    change 'step_overshoot_f_above_nonburn_shell' '0.005' "$overshoot"
    change 'step_overshoot_f_below_nonburn_shell' '0.005' "$overshoot"
    change 'step_overshoot_f_above_burn_h_core' '0.005' "$overshoot"
    change 'step_overshoot_f_above_burn_h_shell' '0.005' "$overshoot"
    change 'step_overshoot_f_below_burn_h_shell' '0.005' "$overshoot"
    
    mk
    rn
    mv LOGS/history.data .
    
    change "write_profiles_flag" ".false." ".true."
    change "create_pre_main_sequence_model" ".true." ".false."
    change "load_saved_model" ".false." ".true."
    change "save_model_when_terminate" ".true." ".false."
    change "stop_near_zams" ".true." ".false."
    change "Lnuc_div_L_zams_limit" "0.999d0" "-1"
    change 'relax_initial_Y' '.true.' '.false.'
    change 'relax_initial_Z' '.true.' '.false.'
    change 'which_atm_option' "'simple_photosphere'" "'Eddington_grey'"
    
    if (( $(echo "$M < 1" | bc -l) )); then
        change "profile_interval" "8" "12"
    fi
    if (( $(echo "$M < 0.9" | bc -l) )); then
        change "profile_interval" "12" "16"
    fi
    if (( $(echo "$M < 0.8" | bc -l) )); then
        change "profile_interval" "16" "20"
    fi
    
    if (( $(echo "$diffusion > 0" | bc -l) )); then
       change "do_element_diffusion" ".false." ".true."
       change 'diffusion_class_factor(:)' '1d0' "$diffusion"
    fi
    
    rn
    
    find "LOGS" -maxdepth 1 -type f -name "*.FGONG" | xargs -i \
        --max-procs=$OMP_NUM_THREADS bash -c \
        "echo start {}; fgong2freqs.sh {}; echo end {}"
    
    cd ../..
    
    sleep 1
    Rscript summarize.R "$dirname"
    rm -rf "$dirname"
}

change() { #param initval newval
    sed -i.bak "s/\!$1 = $2/$1 = $3/g" inlist_1.0
    sed -i.bak "s/$1 = $2/$1 = $3/g" inlist_1.0
}

## Parse command line arguments
# takes mass M, helium Y, metallicity Z, and mixing length parameter alpha
while [ "$#" -gt 0 ]; do
  case "$1" in
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;
    -D) diffusion="$2"; shift 2;;
    -f) overshoot="$2"; shift 2;;
    -d) directory="$2"; shift 2;;

    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done

# set defaults if they weren't supplied
if [ -z ${M+x} ]; then M=1; fi
if [ -z ${Y+x} ]; then Y=0.275; fi
if [ -z ${Z+x} ]; then Z=0.018; fi
if [ -z ${alpha+x} ]; then alpha=1.85; fi
if [ -z ${diffusion+x} ]; then diffusion=0; fi
if [ -z ${overshoot+x} ]; then overshoot=0; fi
if [ -z ${directory+x} ]; then directory=simulations; fi

simulate

