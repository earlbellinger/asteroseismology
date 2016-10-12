#!/bin/bash
#### Script for dispatching a stellar evolutionary track simulation 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

################################################################################
### GLOBAL VARIABLES ###########################################################
################################################################################
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
num_process=70 # try to get at least this many fgong files 
init_mesh_delta_coeff=1 # the number in the inlist file 
mesh_delta_coeff=$init_mesh_delta_coeff
ms_mesh_delta=0.4 # the mesh spacing we will use on the main sequence
max_years_for_timestep=-1 # dt_max 
profile_interval=1 # how often a profile should be output 

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
change() { #param initval newval inlist_name
    sed -i.bak "s/$1 = $2/$1 = $3/g" "$4"
    sed -i.bak "s/\!$1/$1/g" "$4"
}

## Deletes all of the calculations if the remove flag is set 
cleanup() {
    if [ $remove -eq 1 ]; then
        rm -rf "$dirname"
    fi
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
if [ -z ${Y+x} ]; then Y=0.268; fi
if [ -z ${Z+x} ]; then Z=0.0198; fi
if [ -z ${alpha+x} ]; then alpha=1.86; fi
if [ -z ${overshoot+x} ]; then overshoot=0.05; fi
if [ -z ${diffusion+x} ]; then diffusion=1; fi
if [ -z ${directory+x} ]; then directory=simulations; fi
if [ -z ${light+x} ]; then light=0; fi
if [ -z ${remove+x} ]; then remove=0; fi

################################################################################
### INITIALIZATION #############################################################
################################################################################
## Make directory and copy over simulation files 
expname="M=$M""_""Y=$Y""_""Z=$Z""_""alpha=$alpha"\
"_""overshoot=$overshoot""_""diffusion=$diffusion"
dirname="$directory/$expname"

mkdir -p "$dirname"
cd "$dirname"
cp -r $scriptdir/mesa/* .
rm -rf LOGS/*

### Set up initial parameters 
change 'initial_mass' '1.0' "$M" "inlist_pms"
change 'initial_y' '-1' "$Y" "inlist_pms"
change 'initial_z' '0.02' "$Z" "inlist_pms"
change "new_Y" '-1' "$Y" "inlist_pms"
change "new_Z" '-1' "$Z" "inlist_pms"

for inl in $( ls | grep "inlist_" ); do
    change "Zbase" '0.02' "$Z" "$inl"
    change 'mixing_length_alpha' '2.1' "$alpha" "$inl"
    
    if (( $(echo "$overshoot > 0" | bc -l) )); then
       change 'step_overshoot_f_above_nonburn_core'  '0.005' "$overshoot" "$inl"
       change 'step_overshoot_f_above_nonburn_shell' '0.005' "$overshoot" "$inl"
       change 'step_overshoot_f_below_nonburn_shell' '0.005' "$overshoot" "$inl"
       change 'step_overshoot_f_above_burn_h_core'   '0.005' "$overshoot" "$inl"
       change 'step_overshoot_f_above_burn_h_shell'  '0.005' "$overshoot" "$inl"
       change 'step_overshoot_f_below_burn_h_shell'  '0.005' "$overshoot" "$inl"
       
       f0=$(echo "scale=8; $overshoot / 5" | bc -l)
       change 'overshoot_f0_above_nonburn_core'  '0.001' "$f0" "$inl"
       change 'overshoot_f0_above_nonburn_shell' '0.001' "$f0" "$inl"
       change 'overshoot_f0_below_nonburn_shell' '0.001' "$f0" "$inl"
       change 'overshoot_f0_above_burn_h_core'   '0.001' "$f0" "$inl"
       change 'overshoot_f0_above_burn_h_shell'  '0.001' "$f0" "$inl"
       change 'overshoot_f0_below_burn_h_shell'  '0.001' "$f0" "$inl"
    fi
done

if (( $(echo "$diffusion > 0" | bc -l) )); then
   change "do_element_diffusion" ".false." ".true." "inlist_ms"
   change 'diffusion_class_factor(1)' '1' "$diffusion" "inlist_ms"
   change 'diffusion_class_factor(2)' '1' "$diffusion" "inlist_ms"
   change 'diffusion_class_factor(3)' '1' "$diffusion" "inlist_ms"
   change 'diffusion_class_factor(4)' '1' "$diffusion" "inlist_ms"
   change 'diffusion_class_factor(5)' '1' "$diffusion" "inlist_ms"
   if [ $light -eq 1 ]; then
       change 'diffusion_class_representative(4)' "'o16'" "'he4'" "inlist_ms"
       change 'diffusion_class_representative(5)' "'fe56'" "'he4'" "inlist_ms"
   fi
fi

################################################################################
## RUN #########################################################################
################################################################################
./rn | tee "$pmslog"
sed -i.bak "s/\*\.\**-[0-9]*/0.0000000000000000D+00/g" zams.mod 
mv LOGS/history.data history-pms.data

sed -i.bak "s/'inlist_pms'/'inlist_ms'/g" inlist
./rn
sed -i.bak "s/\*\.\**-[0-9]*/0.0000000000000000D+00/g" tams.mod 
cp LOGS/history.data history-ms.data

sed -i.bak "s/'inlist_ms'/'inlist_sg'/g" inlist
./rn
sed -i.bak "s/\*\.\**-[0-9]*/0.0000000000000000D+00/g" brgb.mod 
cp LOGS/history.data history-sg.data

cd LOGS
ls *.FGONG | xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
    "echo start {}; fgong2freqs-giants.sh {}; echo end {}"

