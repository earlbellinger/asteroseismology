#!/bin/bash
#### Script for dispatching a stellar evolutionary track simulation 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

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
    -h) HELP=1; break;;
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
if [ -z ${HELP+x} ]; then HELP=0; fi
if [ -z ${M+x} ]; then M=1; fi
if [ -z ${Y+x} ]; then Y=0.266; fi
if [ -z ${Z+x} ]; then Z=0.018; fi
if [ -z ${alpha+x} ]; then alpha=1.81; fi
if [ -z ${overshoot+x} ]; then overshoot=0.07; fi
if [ -z ${diffusion+x} ]; then diffusion=1; fi
if [ -z ${directory+x} ]; then directory=simulations; fi
if [ -z ${light+x} ]; then light=0; fi
if [ -z ${remove+x} ]; then remove=0; fi

if [ $HELP -gt 0 ]; then
    echo
    echo "   _____ _             _         _____  _                 _       _     ";
    echo "  / ____(_)           | |       |  __ \(_)               | |     | |    ";
    echo " | |  __ _  __ _ _ __ | |_ ___  | |  | |_ ___ _ __   __ _| |_ ___| |__  ";
    echo " | | |_ | |/ _\` | '_ \| __/ __| | |  | | / __| '_ \ / _\` | __/ __| '_ \ ";
    echo " | |__| | | (_| | | | | |_\__ \ | |__| | \__ | |_) | (_| | || (__| | | |";
    echo "  \_____|_|\__,_|_| |_|\__|___/ |_____/|_|___| .__/ \__,_|\__\___|_| |_|";
    echo "                                             | |                        ";
    echo "                                             |_|                        ";
    echo
    echo "  dispatch.sh: evolve an asteroseismic evolutionary track until the "
    echo "               base of the RGB with MESA & GYRE"
    echo
    echo "example:"
    echo "  dispatch.sh -d my_directory -M 1.2 -Z 0.001"
    echo
    echo "flags:"
    echo "  -h   : show this helpful message and quit"
    echo "  -r   : remove the profile files afterwards to save space"
    echo "                                               [default: false      ]"
    echo "  -d s : set the output directory to 's'       [default: simulations]"
    echo "  -L   : only diffuse light elements (X,Y)     [default: false      ]"
    echo
    echo "controls:                             (qty)         (solar defaults)"
    echo "  -M # : initial stellar mass        M/M_sun       [default:  1     ]"
    echo "  -Y # : initial helium abundance    Y_0           [default:  0.266 ]"
    echo "  -Z # : initial metallicity         Z_0           [default:  0.018 ]"
    echo "  -a # : mixing length parameter     \\alpha_MLT    [default:  1.81  ]"
    echo "  -o # : overshooting parameter      \\alpha_ov     [default:  0.07  ]"
    echo "  -D # : diffusion factor            D             [default:  1     ]"
    echo
    exit
fi

################################################################################
### GLOBAL VARIABLES ###########################################################
################################################################################
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
num_process=75 # try to get at least this many profile files 
init_mesh_delta_coeff=1 # the number in the inlist file 
mesh_delta_coeff=$init_mesh_delta_coeff
ms_mesh_delta=0.4 # the mesh spacing we will use on the main sequence
max_years_for_timestep=-1 # dt_max 
min_years_for_timestep=100000
profile_interval=1 # how often a profile should be output 
inlist='inlist_pms'

### Log files 
pmslog="pms.log" 
mslog="ms.log"
sglog="sg.log"
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

## Deletes all of the calculations if the remove flag is set 
cleanup() {
    if [ $remove -eq 1 ]; then
        rm -rf "$dirname"
    fi
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
    num_logs=$(cat $logs_dir/history.data | wc -l)
    if [ $num_logs -lt $num_process ]; then
        age_range=$(R --slave -q -e "options(scipen = 999);"\
"DF <- read.table('$logs_dir/history.data', header=1, skip=5);"\
"cat(max(DF[['star_age']])-min(DF[['star_age']]))")
        max_years_for_timestep=$(echo "$age_range / $num_process" | bc)
        if [ $max_years_for_timestep -lt $min_years_for_timestep ]; then
            echo "Too small timesteps"
            exit 1
        fi
        change 'max_years_for_timestep' '-1' $max_years_for_timestep
    fi
    change 'write_profiles_flag' '.false.' '.true.'
    ./rn
    fix_mod $final_mod
    ls $logs_dir/*.GYRE | xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
        "echo start {}; fgong2freqs-gyre.sh {}; echo end {}"
}

################################################################################
### INITIALIZATION AND PRE-MAIN SEQUENCE #######################################
################################################################################
## Make directory and copy over simulation files 
expname="M=$M""_""Y=$Y""_""Z=$Z""_""alpha=$alpha"\
"_""overshoot=$overshoot""_""diffusion=$diffusion"
dirname="$directory/$expname"

mkdir -p "$dirname"
cd "$dirname"
cp -r $scriptdir/mesa/* .
rm -rf LOGS/*

set_params
./rn
fix_mod "zams.mod"

################################################################################
### MAIN SEQUENCE ##############################################################
################################################################################
change_inlists "inlist_ms"
run "LOGS_MS" "tams.mod" "history-ms.data"

################################################################################
### SUB-GIANT ##################################################################
################################################################################
change_inlists "inlist_sg"
run "LOGS_SG"

## Process pulsation files with GYRE 
#ls LOGS/*.FGONG | xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
#    "echo start {}; fgong2freqs.sh {}; echo end {}"

