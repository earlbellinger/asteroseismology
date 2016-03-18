#!/bin/bash
#### Script for dispatching a stellar evolutionary track simulation 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

################################################################################
### GLOBAL VARIABLES ###########################################################
################################################################################
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
num_process=300 # try to get at least this many fgong files 
init_mesh_delta_coeff=1 # the number in the inlist file 
mesh_delta_coeff=$init_mesh_delta_coeff
ms_mesh_delta=0.4 # the mesh spacing we will use on the main sequence
mesh_delta_limit=0.2 # a lower bound on mesh spacing 
mesh_delta_upper=3 # an upper bound on mesh spacing
max_years_limit=1000000 # a lower bound on time stepping 
max_bounces=20 # quit if we require more than this many mesh adjustments 
n_bounces=0 # the amount of adjustments we've already performed 
max_years_for_timestep=-1
#min_process=100 # minimum number of fgong files 

### Log files 
pmslog="pms.log" 
logfile="mesa.log"

### Success and failure conditions 
msuccess="termination code: max_age\|termination code: xa_central_lower_limit"
pmsuccess="termination code: Lnuc_div_L_zams_limit"
meshfail="mesh_plan problem"
convfail="terminated evolution: convergence problems"

################################################################################
### HELPER FUNCTIONS ###########################################################
################################################################################
## Modifies the inlist. must specify what param you want to change, 
##     the current value it has, and the new value you want it to have
change() { #param initval newval
    if grep -q "!$1 = $2" inlist_1.0; then
        sed -i.bak "s/\!$1 = $2/$1 = $2/g" inlist_1.0
    fi
    sed -i.bak "s/$1 = $2/$1 = $3/g" inlist_1.0
}

## Deletes all of the calculations if the remove flag is set 
cleanup() {
    if [ $remove -eq 1 ]; then
        rm -rf "$dirname"
    fi
}

## Increase global bounce counter and check that it hasn't exceeded max bounces
bounce() {
    n_bounces=$((n_bounces+1))
    if [ $n_bounces -gt $max_bounces ]; then
        echo "Error: Bounced too much"
        cleanup
        exit 1
    fi
}

## By default, increase the number of points
## but if there are convergence problems, decrease instead 
adjust_mesh() {
    new_mesh_delta_coeff=$(echo "scale=3; $mesh_delta_coeff * 0.9" | bc -l)
    if grep -q "$convfail" "$logfile"; then
        new_mesh_delta_coeff=$(echo "scale=3; $mesh_delta_coeff * 1.1" |
            bc -l)
    fi
}

## See that the new mesh spacing is within bounds 
check_mesh_limits() {
    if (( $(echo "$new_mesh_delta_coeff < $mesh_delta_limit" | bc -l) )) ||
       (( $(echo "$new_mesh_delta_coeff > $mesh_delta_upper" | bc -l) ));
      then
        echo "Error: Couldn't achieve MS convergence (mesh limit)" | tee -a "$logfile"
        cleanup
        exit 1
    fi
}

## Actually change the meshing in the inlist and delete all the logs in there
change_meshing() {
    mv LOGS/history.data "history_""$mesh_delta_coeff""_"\
"$max_years_for_timestep"".data"
    rm -f LOGS/*
    change "mesh_delta_coeff" "$mesh_delta_coeff" "$new_mesh_delta_coeff"
    mesh_delta_coeff=$new_mesh_delta_coeff
    echo "Retrying with mesh_delta_coeff = $mesh_delta_coeff" |
        tee "$logfile"
}

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
if [ -z ${Y+x} ]; then Y=0.272; fi
if [ -z ${Z+x} ]; then Z=0.0187; fi
if [ -z ${alpha+x} ]; then alpha=1.936; fi
if [ -z ${overshoot+x} ]; then overshoot=0.05; fi
if [ -z ${diffusion+x} ]; then diffusion=1; fi
if [ -z ${directory+x} ]; then directory=simulations; fi
if [ -z ${light+x} ]; then light=0; fi
if [ -z ${remove+x} ]; then remove=0; fi

################################################################################
### INITIALIZATION AND PRE-MAIN SEQUENCE #######################################
################################################################################
## Make directory and copy over simulation files 
expname="M=$M""_""Y=$Y""_""Z=$Z""_""alpha=$alpha"\
"_""overshoot=$overshoot""_""diffusion=$diffusion"
dirname="$directory/$expname"

mkdir -p "$dirname"
cd "$dirname"
cp -r $scriptdir/mesa_template/* .

### Set up initial parameters 
change 'initial_mass' '1.0' "$M"
change 'initial_y' '-1' "$Y"
change 'initial_z' '0.02' "$Z"
change "new_Y" '-1' "$Y"
change "new_Z" '-1' "$Z"
change 'mixing_length_alpha' '2.1' "$alpha"
#change 'Zbase' '0.02' "$Z"

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

## Obtain pre-main sequence model 
./rn | tee "$pmslog"

## Check that the PMS converged. If not, reduce mesh spacing 
while ! grep -q "$pmsuccess" "$pmslog"; do
    new_mesh_delta_coeff=$(echo "scale=3; $mesh_delta_coeff * 0.9" | bc -l)
    if (( $(echo "$new_mesh_delta_coeff < $mesh_delta_limit" | bc -l) )); 
      then
        echo "Error: Couldn't achieve PMS convergence" | tee -a "$pmslog"
        cleanup
        exit 1
    fi
    echo "Retrying PMS with mesh_delta_coeff = $new_mesh_delta_coeff" | 
        tee "$pmslog"
    change "mesh_delta_coeff" "$mesh_delta_coeff" "$new_mesh_delta_coeff"
    mesh_delta_coeff=$new_mesh_delta_coeff
    ./rn | tee -a "$pmslog"
done
mv LOGS/history.data history-pms.data

################################################################################
### MAIN SEQUENCE ##############################################################
################################################################################
## Turn off PMS conditions
change "create_pre_main_sequence_model" ".true." ".false."
change "load_saved_model" ".false." ".true."
change "save_model_when_terminate" ".true." ".false."
change "stop_near_zams" ".true." ".false."
change "Lnuc_div_L_zams_limit" "0.999d0" "-1"
change 'relax_initial_Y' '.true.' '.false.'
change 'relax_initial_Z' '.true.' '.false.'

## Set up main-sequence conditions like diffusion 
if (( $(echo "$diffusion > 0" | bc -l) )); then
   change "do_element_diffusion" ".false." ".true."
   change 'diffusion_class_factor(1)' '1' "$diffusion"
   change 'diffusion_class_factor(2)' '1' "$diffusion"
   change 'diffusion_class_factor(3)' '1' "$diffusion"
   change 'diffusion_class_factor(4)' '1' "$diffusion"
   change 'diffusion_class_factor(5)' '1' "$diffusion"
   if [ $light -eq 1 ]; then
       change 'diffusion_class_representative(4)' "'o16'" "'he4'"
       change 'diffusion_class_representative(5)' "'fe56'" "'he4'"
   fi
fi
if (( $(echo "$mesh_delta_coeff > $ms_mesh_delta" | bc -l) )); then
    new_mesh_delta_coeff=$ms_mesh_delta
    change 'mesh_delta_coeff' "$mesh_delta_coeff" "$ms_mesh_delta"
fi
change "min_timestep_limit" "1d-12" "1d10"
change 'which_atm_option' "'simple_photosphere'" "'Eddington_grey'"

## Run main sequence without any time limit to obtain max age 
timeout 12h ./rn | tee "$logfile"

## Check that it succeeded. If not, adjust mesh until it does. 
while ! grep -q "$msuccess" "$logfile"; do
    bounce
    adjust_mesh
    check_mesh_limits
    change_meshing
    timeout 12h ./rn | tee "$logfile"
done

## Obtain max age and rerun with enough points 
max_age=$(R --slave -q -e "options(scipen = 999);cat("\
"max(read.table('LOGS/history.data',header=1, skip=5)[['star_age']]))")
timestep=$(echo "scale=0; $max_age / $num_process" | bc -l)
if (( $(echo "$timestep < $max_years_limit" | bc -l) )); then
    max_years_limit=$(echo "scale=0; timestep / 2" | bc -l)
fi
change "max_years_for_timestep" "-1" "$timestep"
change "write_profiles_flag" ".false." ".true."
timeout 24h ./rn | tee "$logfile"

## Check that it worked and that there are no discontinuities 
while ! grep -q "$msuccess" "$logfile" ||
        [ $(Rscript ../../discontinuity.R) == 1 ]; do
    bounce
    adjust_mesh
    
    ## Check that we're still within time bounds 
    if (( $(echo "$max_years_for_timestep < $max_years_limit" | bc -l) )); 
      then
        echo "Error: Couldn't achieve MS convergence (timestep limit)" | tee -a "$logfile"
        cleanup
        exit 1
    fi
    
    ## If the meshing failed, decrease timestep and reset meshing 
    if (( $(echo "$new_mesh_delta_coeff < $mesh_delta_limit" | bc -l) )) ||
            grep -q "$meshfail" "$logfile"; then
        new_max_years_for_timestep=$(echo "scale=0;
            $max_years_for_timestep/2" | bc -l)
        new_mesh_delta_coeff=$init_mesh_delta_coeff
        change "max_years_for_timestep" "$max_years_for_timestep" \
            "$new_max_years_for_timestep" 
        max_years_for_timestep=$new_max_years_for_timestep
        echo "Retrying with max_years_for_timestep = $max_years_for_timestep" |
            tee -a "$logfile"
    fi
    
    check_mesh_limits
    change_meshing
    timeout 24h ./rn | tee -a "$logfile"
done

## Copy the working history file 
cp LOGS/history.data history-ms.data

## Process FGONG files with ADIPLS 
ls LOGS/*.FGONG | xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
    "echo start {}; fgong2freqs.sh {}; echo end {}"

## Summarize what we've learned and finish 
cd ../..
Rscript summarize.R "$dirname"
cleanup

