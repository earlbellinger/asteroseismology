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
num_process=128 # try to get at least this many profile files 
min_years_for_timestep=0.0000001
inlist='inlist_1pms'
continue=1

M0=1
Y0=0.27202387
Z0=0.01830403
solar_alpha=1.84663590
solar_overshoot=0.09104194

sigmoida=$(echo "scale=10; (1+e(-6))/(1-e(-6))" | bc -l)
sigmoidb=$(echo "scale=10; 1/(e(6)-1)" | bc -l)
sigmoid() {
    echo "scale=10; $sigmoida / (1+e(6*(2*$1 - 1))) - $sigmoidb" | bc -l
}


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
    -s) suppress=1; shift 1;;
    
    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done

## Set defaults if they weren't supplied
if [ -z ${HELP+x} ]; then HELP=0; fi
if [ -z ${M+x} ]; then M=$M0; fi
if [ -z ${Y+x} ]; then Y=$Y0; fi
if [ -z ${Z+x} ]; then Z=$Z0; fi
if [ -z ${alpha+x} ]; then alpha=$solar_alpha; fi
if [ -z ${o+x} ]; then overshoot=$solar_overshoot; fi
if [ -z ${directory+x} ]; then directory=simulations; fi
if [ -z ${light+x} ]; then light=0; fi
if [ -z ${remove+x} ]; then remove=0; fi
if [ -z ${suppress+x} ]; then suppress=0; fi

if [ -z ${diffusion} ]; then
    if (( $(echo "$M <= 1.2" | bc -l ) )); then 
        diffusion=1
    elif (( $(echo "$M >= 1.3" | bc -l ) )); then 
        diffusion=0
    else
        x=$(echo "scale=10; 1 - ( 1.3 - $M )*10" | bc -l)
        diffusion=$(sigmoid $x)
    fi
fi 

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
    echo "  -Y # : initial helium abundance    Y_0           [default:  0.27  ]"
    echo "  -Z # : initial metallicity         Z_0           [default:  0.018 ]"
    echo "  -D # : diffusion factor            D             [default:  1     ]"
    echo
    exit
fi

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
    change 'mixing_length_alpha' "$alpha" "$inlist"
    
    if [[ suppress -eq 1 ]]; then
        change 'write_profiles_flag' '.false.' "$inlist"
    fi
    
    set_overshoot
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
        
        f0=$(echo "scale=8; $overshoot / 5" | bc -l)
        change 'overshoot_f0_above_nonburn_core' "$f0" "$inlist"
        change 'overshoot_f0_above_nonburn_shell' "$f0" "$inlist"
        change 'overshoot_f0_below_nonburn_shell' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_h_core' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_h_shell' "$f0" "$inlist"
        change 'overshoot_f0_below_burn_h_shell' "$f0" "$inlist"
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
    if [[ continue -eq 0 ]]; then 
        return 0
    fi
    
    logs_dir=$1
    final_mod=$2
    
    ./rn | tee output
    num_logs=$(cat $logs_dir/history.data | wc -l)
    num_logs=$(echo "$num_logs - 7" | bc -l)
    while [ $num_logs -lt $num_process ]; do
        
        if ! egrep -q "history_interval = 1\s*$" $inlist; then
            hist_int=$(grep "history_interval" $inlist |
                sed -e 's/[^0-9]//g')
            interval=$(echo "scale=0;"\
                "$hist_int * $num_logs / $num_process" | bc -l)
            if [ $interval -lt 1 ]; then 
                interval=1
            fi
            change 'profile_interval' $interval $inlist
            change 'history_interval' $interval $inlist
        fi
        
        if egrep -q "history_interval = 1\s*$" $inlist; then
            timestep=$(R --slave -q -e "options(scipen = 999);"\
"DF <- read.table('"$logs_dir"/history.data', header=1, skip=5);"\
"res <- (max(DF[['star_age']])-min(DF[['star_age']]))/"$num_process";"\
"res.str <- sub('\\\..+', '', paste(res));"\
"cat(sub('e\\\+0+|e0+', 'd', res.str))")
            if (( $(echo "$timestep < $min_years_for_timestep" | bc -l) )); then
                echo "Too small timesteps"
                exit 1
            fi
            change 'max_years_for_timestep' "$timestep" $inlist
        fi
        
        rm -rf "$logs_dir"
        ./rn | tee output
        num_logs=$(cat $logs_dir/history.data | wc -l)
        num_logs=$(echo "$num_logs - 7" | bc -l)
    done
    
    fix_mod $final_mod
    
    if [[ suppress -eq 0 ]]; then
        Rscript $scriptdir/model_select.R $num_process $logs_dir | 
            xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
                "echo start {}; gyre-l0.sh {}; echo end {}"
                #"echo start {}; gyre-Dnu.sh {}; echo end {}"
    fi
    
    if grep -q "termination code: max_age" output; then 
        continue=0
    fi
}

################################################################################
### INITIALIZATION AND EVOLUTION ###############################################
################################################################################
## Make directory and copy over simulation files 
expname="M=$M""_""Y=$Y""_""Z=$Z"
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
set_inlist "inlist_2ms"
set_diffusion
run "LOGS_MS" "tams.mod"

# sub-giant
set_inlist "inlist_3sg"
set_diffusion
run "LOGS_SG" "brgb.mod"

# red giant
set_inlist "inlist_4rgb"
run "LOGS_RGB" "bump.mod"

# bump to tip
set_inlist "inlist_5bump"
run "LOGS_BUMP" "flash.mod"

# helium burning
set_inlist "inlist_6heb"
run "LOGS_HEB" "heb.mod"

# summarize and be done 
cd ../..
Rscript summarize.R "$dirname"
cleanup

# fin.

