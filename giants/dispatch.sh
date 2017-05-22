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
continue=1

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
if [ -z ${Y+x} ]; then Y=0.26929755; fi
if [ -z ${Z+x} ]; then Z=0.01758726; fi
if [ -z ${alpha+x} ]; then alpha=1.80003929; fi
if [ -z ${overshoot+x} ]; then overshoot=0.08212582; fi
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
    elif grep -q "\! $param = $initval" "$inlist"; then
        sed -i.bak "s/\! $param = $initval/$param = $initval/g" "$inlist"
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
        change 'step_overshoot_f_above_burn_he_core' '0.005' "$overshoot"
        
        f0=$(echo "scale=8; $overshoot / 5" | bc -l)
        change 'overshoot_f0_above_nonburn_core' '0.001' "$f0"
        change 'overshoot_f0_above_nonburn_shell' '0.001' "$f0"
        change 'overshoot_f0_below_nonburn_shell' '0.001' "$f0"
        change 'overshoot_f0_above_burn_h_core' '0.001' "$f0"
        change 'overshoot_f0_above_burn_h_shell' '0.001' "$f0"
        change 'overshoot_f0_below_burn_h_shell' '0.001' "$f0"
        change 'overshoot_f0_above_burn_he_core' '0.001' "$f0"
    fi
}

set_diffusion() {
    if (( $(echo "$diffusion > 0" | bc -l) )); then
        change "do_element_diffusion" ".false." ".true." "$inlist"
        change 'diffusion_class_factor(1)' '1' "$diffusion" "$inlist"
        change 'diffusion_class_factor(2)' '1' "$diffusion" "$inlist"
        change 'diffusion_class_factor(3)' '1' "$diffusion" "$inlist"
        change 'diffusion_class_factor(4)' '1' "$diffusion" "$inlist"
        change 'diffusion_class_factor(5)' '1' "$diffusion" "$inlist"
        decrease=$(echo "scale=8; $diffusion / 100" | bc -l)
        change 'x_ctrl(2)' '0.01' "$decrease" "$inlist"
    fi
}

fix_mod() { # final_mod
    final_mod=$1
    # MESA has a bug that makes it print numbers of the form *.***### 
    # (i.e. no letter for the exponent, and all the numbers are asterisks) 
    # This sed command replaces those with zeros. 
    sed -i.bak "s/\*\.\**-[0-9]*/0.0000000000000000D+00/g" $final_mod
}

run() { # logs_dir  final_mod
    if ((!continue)); then return 0; fi
    
    logs_dir=$1
    final_mod=$2
    
    ./rn | tee output
    num_logs=$(cat $logs_dir/history.data | wc -l)
    while [ $num_logs -lt $num_process ]; do
        age_range=$(R --slave -q -e "options(scipen = 999);"\
"DF <- read.table('"$logs_dir"/history.data', header=1, skip=5);"\
"cat(max(DF[['star_age']])-min(DF[['star_age']]))")
        max_years_for_timestep=$(echo "$age_range / $num_process" | bc)
        if [ $max_years_for_timestep -lt $min_years_for_timestep ]; then
            echo "Too small timesteps"
            exit 1
        fi
        change 'max_years_for_timestep' '.\+' $max_years_for_timestep
        rm -rf "$logs_dir"
        ./rn
        num_logs=$(cat $logs_dir/history.data | wc -l)
    done
    #change 'write_profiles_flag' '.false.' '.true.'
    #./rn
    fix_mod $final_mod
    Rscript $scriptdir/model_select.R $num_process $logs_dir | 
        xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
            "echo start {}; fgong2freqs-gyre.sh {}; echo end {}"
    
    if grep -q "termination code: max_age" output; then continue=0; fi
}

################################################################################
### INITIALIZATION AND EVOLUTION ###############################################
################################################################################
## Make directory and copy over simulation files 
expname="M=$M""_""Y=$Y""_""Z=$Z""_""alpha=$alpha"\
"_""overshoot=$overshoot""_""diffusion=$diffusion"
dirname="$directory/$expname"

mkdir -p "$dirname"
cd "$dirname"
cp -r $scriptdir/mesa/* .
rm -rf LOGS/*

# pre-main sequence
set_params
./rn
fix_mod "zams.mod"

# main sequence
change_inlists "inlist_ms"
set_diffusion
run "LOGS_MS" "tams.mod"

# sub-giant
change_inlists "inlist_sg"
set_diffusion
run "LOGS_SG" "brgb.mod"

# summarize and be done 
cd ../..
Rscript summarize.R "$dirname"
cleanup

# fin.
