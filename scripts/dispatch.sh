#!/bin/bash
#### Script for dispatching a stellar evolutionary track simulation 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
min_process=10
num_process=200
init_mesh_delta_coeff=1
mesh_delta_coeff=$init_mesh_delta_coeff
mesh_delta_limit=0.3
mesh_delta_upper=3
max_years_for_timestep=1000000
max_years_limit=1000000
max_bounces=10
n_bounces=0

pmslog="pms.log"
logfile="mesa.log"
msuccess="termination code: max_age\|termination code: xa_central_lower_limit"
pmsuccess="termination code: Lnuc_div_L_zams_limit"

meshfail="mesh_plan problem"
convfail="terminated evolution: convergence problems"

simulate() {
    expname="M=$M""_""Y=$Y""_""Z=$Z""_""alpha=$alpha"\
"_""overshoot=$overshoot""_""diffusion=$diffusion"
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
    
    if (( $(echo "$overshoot > 0" | bc -l) )); then
        change 'step_overshoot_f_above_nonburn_core' '0.005' "$overshoot"
        change 'step_overshoot_f_above_nonburn_shell' '0.005' "$overshoot"
        change 'step_overshoot_f_below_nonburn_shell' '0.005' "$overshoot"
        change 'step_overshoot_f_above_burn_h_core' '0.005' "$overshoot"
        change 'step_overshoot_f_above_burn_h_shell' '0.005' "$overshoot"
        change 'step_overshoot_f_below_burn_h_shell' '0.005' "$overshoot"
    fi
    
    ./rn | tee "$pmslog"
    
    while ! grep -q "$pmsuccess" "$pmslog"; do
        new_mesh_delta_coeff=$(echo "scale=2; $mesh_delta_coeff * 0.9" | bc -l)
        if (( $(echo "$new_mesh_delta_coeff < $mesh_delta_limit" | bc -l) )); 
          then
            echo "Couldn't achieve convergence" | tee -a "$pmslog"
            cleanup
            exit 1
        fi
        echo "Retrying PMS with mesh_delta_coeff = $new_mesh_delta_coeff" | 
            tee "$pmslog"
        change "mesh_delta_coeff" "$mesh_delta_coeff" "$new_mesh_delta_coeff"
        mesh_delta_coeff=$new_mesh_delta_coeff
        ./rn | tee -a "$pmslog"
    done
    
    mv LOGS/history.data .
    
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
    
    change "max_years_for_timestep" "-1" "$max_years_for_timestep"
    change "min_timestep_limit" "1d-12" "1d9"
    change "create_pre_main_sequence_model" ".true." ".false."
    change "load_saved_model" ".false." ".true."
    change "save_model_when_terminate" ".true." ".false."
    change "stop_near_zams" ".true." ".false."
    change "Lnuc_div_L_zams_limit" "0.999d0" "-1"
    change 'relax_initial_Y' '.true.' '.false.'
    change 'relax_initial_Z' '.true.' '.false.'
    change 'which_atm_option' "'simple_photosphere'" "'Eddington_grey'"
    
    ./rn | tee "$logfile"
    
    while ! grep -q "$msuccess" "$logfile" ||
            [ $(Rscript ../../discontinuity.R) == 1 ]; do
        
        tmp=$(Rscript ../../profile_interval.R)
        tmp=(${tmp//;/ })
        num_lines=${tmp[0]}
        max_age=${tmp[1]}
        if [ $num_lines -lt $min_process ] && 
           [ $(cat LOGS/history.data | wc -l) -gt $num_process ]; then
            echo "Fully radiative track"
            cleanup
            exit 1
        fi
        
        n_bounces=$n_bounces+1
        if [ $n_bounces -gt $max_bounces ]; then
            echo "Bounced too much"
            cleanup
            exit 1
        fi
        
        mv LOGS/history.data "history_""$mesh_delta_coeff""_"\
"$max_years_for_timestep"".data"
        
        # decrease the mesh spacing
        new_mesh_delta_coeff=$(echo "scale=2; $mesh_delta_coeff * 0.9" | bc -l)
        
        # check that we're still within bounds
        if (( $(echo "$new_mesh_delta_coeff < $mesh_delta_limit" | bc -l) &&
              $(echo "$max_years_for_timestep < $max_years_limit" | bc -l) )) ||
           (( $(echo "$new_mesh_delta_coeff > $mesh_delta_upper" | bc -l) ));
          then
            echo "Couldn't achieve convergence" | tee -a "$logfile"
            cleanup
            exit 1
        fi
        
        # if the meshing failed, decrease timestep and reset meshing
        if (( $(echo "$new_mesh_delta_coeff < $mesh_delta_limit" | bc -l) )) ||
                grep -q "$meshfail" "$logfile"; then
            new_max_years_for_timestep=$(echo "scale=0;
                $max_years_for_timestep/2" | bc -l)
            new_mesh_delta_coeff=$init_mesh_delta_coeff
            change "max_years_for_timestep" "$max_years_for_timestep" \
                "$new_max_years_for_timestep" 
            max_years_for_timestep=$new_max_years_for_timestep
        fi
        
        # if there are convergence problems, decrease the number of points used
        if grep -q "$convfail" "$logfile"; then
            new_mesh_delta_coeff=$(echo "scale=2; $mesh_delta_coeff * 1.1" |
                bc -l)
        fi
        
        change "mesh_delta_coeff" "$mesh_delta_coeff" "$new_mesh_delta_coeff"
        mesh_delta_coeff=$new_mesh_delta_coeff
        echo "Retrying with mesh_delta_coeff = $mesh_delta_coeff" |
            tee "$logfile"
        echo "Retrying with max_years_for_timestep = $max_years_for_timestep" |
            tee -a "$logfile"
        ./rn | tee -a "$logfile"
    done
    
    # enable profile writing and rerun track with good settings 
    change "write_profiles_flag" ".false." ".true."
    
    tmp=$(Rscript ../../profile_interval.R)
    tmp=(${tmp//;/ })
    num_lines=${tmp[0]}
    max_age=${tmp[1]}
    if [ $num_lines -lt 10 ]; then
        echo "Fully radiative track"
        cleanup
        exit 1
    fi
    if [ $num_lines -gt $num_process ]; then
        new_profile_interval=$(echo "scale=0; 
            $num_lines / $num_process" | bc -l)
        if [ $new_profile_interval -gt 1 ]; then
            change "profile_interval" "1" "$new_profile_interval"
        fi
    fi
    change "max_age" "16e9" "$max_age"
    ./rn | tee -a "$logfile"
    
    # only process some of the adipls files 
    num_files="$(find 'LOGS' -maxdepth 1 -type f -name '*.FGONG' | wc -l)"
    if [ $num_files -lt $min_process ]; then
        echo "Small number of logs generated ($num_files < $min_process)" | 
            tee -a "$logfile"
        #exit 1 # maybe some day I will change this to rerun with more profiles
    fi
    Rscript ../../fgong_enumerate.R "$num_process" | 
        xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
        "echo start {}; fgong2freqs.sh {}; echo end {}"
    
    cd ../..
    Rscript summarize.R "$dirname"
    cleanup
}

change() { #param initval newval
    if grep -q "!$1 = $2" inlist_1.0; then
        sed -i.bak "s/\!$1 = $2/$1 = $2/g" inlist_1.0
    fi
    sed -i.bak "s/$1 = $2/$1 = $3/g" inlist_1.0
    sleep 1
}

cleanup() {
    if [ $remove -eq 1 ]; then
        rm -rf "$dirname"
    fi
}

## Parse command line arguments
# takes mass M, helium Y, metallicity Z, and mixing length parameter alpha
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

# set defaults if they weren't supplied
if [ -z ${M+x} ]; then M=1; fi
if [ -z ${Y+x} ]; then Y=0.27; fi
if [ -z ${Z+x} ]; then Z=0.018; fi
if [ -z ${alpha+x} ]; then alpha=1.9; fi
if [ -z ${overshoot+x} ]; then overshoot=0.2; fi
if [ -z ${diffusion+x} ]; then diffusion=1; fi
if [ -z ${directory+x} ]; then directory=simulations; fi
if [ -z ${light+x} ]; then light=0; fi
if [ -z ${remove+x} ]; then remove=0; fi

simulate
