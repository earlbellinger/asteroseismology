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
num_process=16 # try to get at least this many profile files 
min_years_for_timestep=0.0000001
inlist='inlist_1pms'
continue=1

M0=1
Y0=0.27202387
Z0=0.01830403
alpha0=1.81
overshoot0=0.3

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
if [ -z ${M+x} ]; then M=$M0; fi
if [ -z ${Y+x} ]; then Y=$Y0; fi
if [ -z ${Z+x} ]; then Z=$Z0; fi
if [ -z ${alpha+x} ]; then alpha=$alpha0; fi
if [ -z ${overshoot+x} ]; then overshoot=$overshoot0; fi
if [ -z ${directory+x} ]; then directory=simulations; fi
if [ -z ${light+x} ]; then light=0; fi
if [ -z ${remove+x} ]; then remove=0; fi
if [ -z ${suppress+x} ]; then suppress=0; fi

if [ -z ${diffusion+x} ]; then 
    if (( $(echo "$M <= 1.2" | bc -l ) )); then 
        diffusion=1
    elif (( $(echo "$M >= 1.3" | bc -l ) )); then 
        diffusion=0
    else
        #diffusion=$(echo "( 1.3 - $M )*10" | bc -l)
        x=$(echo "scale=10; 1 - ( 1.3 - $M )*10" | bc -l)
        diffusion=$(sigmoid $x)
    fi
fi


################################################################################
### HELPER FUNCTIONS ###########################################################
################################################################################
change_param() { 
    # Modifies a parameter in the current inlist. 
    # args: ($1) name of parameter 
    #       ($2) new value 
    #       ($3) filename of inlist where change should occur 
    # Additionally changes the 'inlist_0all' inlist. 
    # example command: change_param initial_mass 1.3 
    # example command: change_param log_directory 'LOGS_MS' 
    # example command: change_param do_element_diffusion .true. 
    param=$1 
    newval=$2 
    filename=$3 
    escapedParam=$(sed 's/[^^]/[&]/g; s/\^/\\^/g' <<< "$param")
    search="^\s*\!*\s*$escapedParam\s*=.+$" 
    replace="      $param = $newval" 
    sed -r -i.bak -e "s/$search/$replace/g" $filename 
    if [ ! "$filename" == 'inlist_0all' ]; then 
        change_param $param $newval 'inlist_0all'
    fi 
} 

set_inlist() { 
    # Changes to a different inlist by modifying where "inlist" file points 
    # args: ($1) filename of new inlist  
    # example command: change_inlists inlist_2ms 
    newinlist=$1 
    echo "Changing to $newinlist" 
    change_param "extra_star_job_inlist2_name" "'$newinlist'" "inlist" 
    change_param "extra_controls_inlist2_name" "'$newinlist'" "inlist" 
    inlist=$newinlist
}

## Deletes all of the calculations if the remove flag is set 
cleanup() {
    if [ $remove -eq 1 ]; then
        rm -rf "$dirname"
    fi
}

set_params() {
    change_param 'initial_mass' "$M" "$inlist"
    change_param 'initial_y' "$Y" "$inlist"
    change_param 'initial_z' "$Z" "$inlist"
    change_param 'mixing_length_alpha' "$alpha" "$inlist"
    change_param 'new_Y' "$Y" "$inlist"
    change_param 'new_Z' "$Z" "$inlist"
    change_param 'Zbase' "$Z" "$inlist"
    
    change_param 'step_overshoot_f_above_nonburn_core' "$overshoot" "$inlist"
    change_param 'step_overshoot_f_above_nonburn_shell' "$overshoot" "$inlist"
    change_param 'step_overshoot_f_below_nonburn_shell' "$overshoot" "$inlist"
    change_param 'step_overshoot_f_above_burn_h_core' "$overshoot" "$inlist"
    change_param 'step_overshoot_f_above_burn_h_shell' "$overshoot" "$inlist"
    change_param 'step_overshoot_f_below_burn_h_shell' "$overshoot" "$inlist"
    change_param 'step_overshoot_f_above_burn_he_core' "$overshoot" "$inlist"
    change_param 'step_overshoot_f_above_burn_he_shell' "$overshoot" "$inlist"
    change_param 'step_overshoot_f_below_burn_he_shell' "$overshoot" "$inlist"
    f0=$(echo "scale=8; $overshoot / 5" | bc -l)
    change_param 'overshoot_f0_above_nonburn_core' "$f0" "$inlist"
    change_param 'overshoot_f0_above_nonburn_shell' "$f0" "$inlist"
    change_param 'overshoot_f0_below_nonburn_shell' "$f0" "$inlist"
    change_param 'overshoot_f0_above_burn_h_core' "$f0" "$inlist"
    change_param 'overshoot_f0_above_burn_h_shell' "$f0" "$inlist"
    change_param 'overshoot_f0_below_burn_h_shell' "$f0" "$inlist"
    change_param 'overshoot_f0_above_burn_he_core' "$f0" "$inlist"
    change_param 'overshoot_f0_above_burn_he_shell' "$f0" "$inlist"
    change_param 'overshoot_f0_below_burn_he_shell' "$f0" "$inlist"
    
    if [[ suppress -eq 1 ]]; then
        change_param 'write_profiles_flag' '.false.' "$inlist"
    fi
}

set_diffusion() {
    if (( $(echo "$diffusion > 0" | bc -l) )); then
        change_param 'do_element_diffusion' '.true.' "$inlist"
        change_param 'diffusion_class_factor(1)' "$diffusion" "$inlist"
        change_param 'diffusion_class_factor(2)' "$diffusion" "$inlist"
        change_param 'diffusion_class_factor(3)' "$diffusion" "$inlist"
        change_param 'diffusion_class_factor(4)' "$diffusion" "$inlist"
        change_param 'diffusion_class_factor(5)' "$diffusion" "$inlist"
        decrease=$(echo "scale=8; $diffusion / 100" | bc -l)
        change_param 'x_ctrl(2)' "$decrease" "$inlist"
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
    if grep -q "cannot find acceptable model" output; then
        rm -rf "LOGS_MS" 
        exit 1 
    fi
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
            change_param 'profile_interval' $interval $inlist
            change_param 'history_interval' $interval $inlist
        fi
        
        if egrep -q "history_interval = 1\s*$" $inlist; then
            timestep=$(R --slave -q -e "options(scipen = 999);"\
"DF <- read.table('"$logs_dir"/history.data', header=1, skip=5);"\
"res <- (max(DF[['star_age']])-min(DF[['star_age']]))/"$num_process";"\
"res.str <- sub('\\\..+', '', paste(res));"\
"cat(sub('e\\\+0+|e0+', 'd', res.str))")
            if (( $(echo "$timestep < $min_years_for_timestep" | bc -l) )); then
                echo "Too small timesteps"
                rm -rf "LOGS_MS"
                exit 1
            fi
            change_param 'max_years_for_timestep' "$timestep" $inlist
        fi
        
        rm -rf "$logs_dir"
        ./rn | tee output
        num_logs=$(cat $logs_dir/history.data | wc -l)
        num_logs=$(echo "$num_logs - 7" | bc -l)
        
        if grep -q "cannot find acceptable model" output; then
            rm -rf "LOGS_MS" 
            exit 1 
        fi
    done
    
    fix_mod $final_mod
    
    if [[ suppress -eq 0 ]]; then
        Rscript $scriptdir/model_select.R $num_process $logs_dir | 
            xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
                "echo start {}; 
                 fgong2freqs.sh {}; 
                 python3 $SCRIPTS_DIR/fgong2ascii.py -i {};
                 echo end {}"
        cd $logs_dir
        find . -d -type d -exec rm -rf '{}' \;
        cd -
    fi
    
    if grep -q "termination code: max_age" output; then 
        continue=0
    fi
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
set_inlist "inlist_2ms"
set_diffusion
run "LOGS_MS" "tams.mod"

# sub-giant
#set_inlist "inlist_3sg"
#set_diffusion
#run "LOGS_SG" "brgb.mod"

# summarize and be done 
cd ../..
#Rscript summarize.R "$dirname"

# fin.
