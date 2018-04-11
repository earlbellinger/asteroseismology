#!/bin/bash
#### Script for dispatching a stellar evolutionary track simulation 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

################################################################################
### GLOBAL VARIABLES ###########################################################
################################################################################
export MESA_DIR=/scratch/seismo/bellinger/MESA/mesa-r10108

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
min_years_for_timestep=0.0000001
inlist='inlist_1relax'
proceed=1

#sigmoida=$(echo "scale=10; (1+e(-6))/(1-e(-6))" | bc -l)
#sigmoidb=$(echo "scale=10; 1/(e(6)-1)" | bc -l)
#sigmoid() {
#    echo "scale=10; $sigmoida / (1+e(6*(2*$1 - 1))) - $sigmoidb" | bc -l
#}


################################################################################
### PARSE COMMAND LINE ARGUMENTS ###############################################
################################################################################
M=1
Y=0.272957104887671
Z=0.0185827894799524
alpha=1.8367593860737
overshoot=0
undershoot=0
overexp=0
underexp=0
diffusion=1
settling=1
eta=0
HELP=0
directory=simulations
taper=0
remove=0
suppress=0
calibrate=0
mainseq=0
subgiant=0
bump=0
nopulse=0
FGONG=0
f0=0.01
num_process=128 # try to get at least this many profile files 
net="'pp_cno_extras_o18_ne22.net'"

while [ "$#" -gt 0 ]; do
  case "$1" in
    -h) HELP=1; break;;
    -n) expname="$2"; shift 2;;
    -N) num_process="$2"; shift 2;;
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;
    -o) overshoot="$2"; shift 2;;
   -oe) overexp="$2"; shift 2;;
    -u) undershoot="$2"; shift 2;;
   -ue) underexp="$2"; shift 2;;
   -f0) f0="$2"; shift 2;;
    -D) diffusion="$2"; shift 2;;
    -g) settling="$2"; shift 2;;
    -e) eta="$2"; shift 2;;
    -d) directory="$2"; shift 2;;
    -c) calibrate="$2"; shift 2;;
  -net) net="$2"; shift 2;;
    -t) taper=1; shift 1;;
    -r) remove=1; shift 1;;
   -MS) mainseq=1; shift 1;;
    -S) subgiant=1; shift 1;;
    -B) bump=1; shift 1;;
    -s) suppress=1; shift 1;;
    -p) nopulse=1; shift 1;;
    -f) FGONG=1; shift 1;;
    
    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done

if (( $(echo "$taper > 0" | bc -l) )); then
    if (( $(echo "$M <= 1.25" | bc -l) )); then
        diffusion=1
        settling=1
    else
        #val=$(echo "scale=10; (1+e(-6))/(1-e(-6))" | bc -l)
        val=$(echo "scale=10; e(-($M-1.25)*($M-1.25) / 0.01445)" | bc -l)
        diffusion=$val
        settling=$val
    fi
fi

#if (( $(echo "$taper > 0" | bc -l) )); then
#    if (( $(echo "$M <= 1.2" | bc -l ) )); then 
#        diffusion=1
#        settling=1
#    elif (( $(echo "$M >= 1.3" | bc -l ) )); then 
#        diffusion=0
#        settling=0
#    else
#        x=$(echo "scale=10; 1 - ( 1.3 - $M )*10" | bc -l)
#        diffusion=$(sigmoid $x)
#        settling=$(sigmoid $x)
#    fi
#fi

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
    echo "  dispatch.sh: evolve an asteroseismic evolutionary track through to "
    echo "               the red clump with MESA & GYRE"
    echo
    echo "example:"
    echo "  dispatch.sh -d my_directory -M 1.2 -Z 0.001"
    echo
    echo "flags:"
    echo "  -h   : show this helpful message and quit"
    echo "  -t   : turn off diffusion as mass increases"
    echo "  -S   : stop at the base of the RGB"
    echo " -MS   : stop at core H exhaustion" 
    echo "  -B   : stop at the RGB bump"
    echo "  -s   : suppress the writing of profile files"
    echo "  -p   : suppress pulsation calculations"
    echo "  -r   : delete profile files"
    echo "  -f   : use FGONG file format (otherwise GYRE)"
    echo "  -c # : calibrate to age # Gyr"
    echo "  -n s : name the track 's'                    [default: from params]"
    echo "  -N # : output # profile files                [default: 128        ]"
    echo "  -d s : set the output directory to 's'       [default: simulations]"
    echo "-net s : use the 's' isotope network           [default:            ]"
    echo
    echo "controls:                             (qty)         (solar defaults)"
    echo "  -M # : initial stellar mass        M/M_sun       [default:  1     ]"
    echo "  -Y # : initial helium abundance    Y_0           [default:  0.2729]"
    echo "  -Z # : initial metallicity         Z_0           [default:  0.0185]"
    echo "  -a # : mixing length parameter     alpha_MLT     [default:  1.8367]"
    echo "  -D # : diffusion factor            D             [default:  1     ]"
    echo "  -g # : gravitational settling      g             [default:  1     ]"
    echo "  -e # : Reimer's mass loss          eta           [default:  0     ]"
    echo "  -o # : overshooting                alpha_ov      [default:  0     ]"
    echo " -oe # : exponential overshooting                  [default:  0     ]"
    echo "  -u # : undershooting                             [default:  0     ]"
    echo " -ue # : exponential undershooting                 [default:  0     ]"
    echo " -f0 # : location from which to get velocity       [default:  0.05  ]"
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
    cd - 
    Rscript summarize.R "$dirname" $num_process #$(echo "$num_process / 2" | bc)
    if [ $remove -eq 1 ]; then
        rm -rf "$dirname"
    fi
    exit
}

set_params() {
    change 'initial_mass' "$M" "$inlist"
    change 'initial_y' "$Y" "$inlist"
    change 'initial_z' "$Z" "$inlist"
    change 'new_Y' "$Y" "$inlist"
    change 'new_Z' "$Z" "$inlist"
    change 'Zbase' "$Z" "$inlist"
    change 'mixing_length_alpha' "$alpha" "$inlist"
    change 'new_net_name' "$net" "$inlist"
    
    if [[ FGONG -eq 1 ]]; then
        change 'pulse_data_format' "'FGONG'" "$inlist"
    fi
    
    if [[ suppress -eq 1 ]]; then
        change 'write_profiles_flag' '.false.' "$inlist"
    fi
    if [[ nopulse -eq 1 ]]; then
        change 'write_pulse_data_with_profile' '.false.' "$inlist"
    fi
}

set_diffusion() {
    if (( $(echo "$diffusion > 0" | bc -l) )); then
        change 'do_element_diffusion' '.true.' "$inlist"
        change 'diffusion_SIG_factor' "$diffusion" "$inlist"
        change 'diffusion_GT_factor' "$settling" "$inlist" 
        #change 'diffusion_class_factor(:)' "$diffusion" "$inlist"
        #change 'diffusion_class_factor(1)' "$diffusion" "$inlist"
        #change 'diffusion_class_factor(2)' "$diffusion" "$inlist"
        #change 'diffusion_class_factor(3)' "$diffusion" "$inlist"
        #change 'diffusion_class_factor(4)' "$diffusion" "$inlist"
        #change 'diffusion_class_factor(5)' "$diffusion" "$inlist"
    fi
    if (( $(echo "$settling > 0" | bc -l) )); then
        change 'do_element_diffusion' '.true.' "$inlist"
        change 'diffusion_SIG_factor' "$diffusion" "$inlist"
        change 'diffusion_GT_factor' "$settling" "$inlist" 
    fi
}

set_overshoot() {
    if (( $(echo "$overshoot > 0" | bc -l) )); then
        change 'step_overshoot_f_above_nonburn_core' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_nonburn_shell' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_burn_h_core' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_burn_h_shell' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_burn_he_core' "$overshoot" "$inlist"
        change 'step_overshoot_f_above_burn_he_shell' "$overshoot" "$inlist"
        
        #f0=$(echo "scale=8; $overshoot / 5" | bc -l)
        #change 'overshoot_f0_above_nonburn_core' "$f0" "$inlist"
        #change 'overshoot_f0_above_nonburn_shell' "$f0" "$inlist"
        #change 'overshoot_f0_above_burn_h_core' "$f0" "$inlist"
        #change 'overshoot_f0_above_burn_h_shell' "$f0" "$inlist"
        #change 'overshoot_f0_above_burn_he_core' "$f0" "$inlist"
        #change 'overshoot_f0_above_burn_he_shell' "$f0" "$inlist"
        #change 'overshoot_f0_below_nonburn_shell' "$f0" "$inlist"
        #change 'overshoot_f0_below_burn_h_shell' "$f0" "$inlist"
        #change 'overshoot_f0_below_burn_he_shell' "$f0" "$inlist"
    fi
    
    if (( $(echo "$undershoot > 0" | bc -l) )); then
        change 'step_overshoot_f_below_nonburn_shell' "$undershoot" "$inlist"
        change 'step_overshoot_f_below_burn_h_shell' "$undershoot" "$inlist"
        change 'step_overshoot_f_below_burn_he_shell' "$undershoot" "$inlist"
    fi
    
    if (( $(echo "$overexp > 0" | bc -l) )); then
        change 'overshoot_f_above_nonburn_core' "$overexp" "$inlist"
        change 'overshoot_f_above_nonburn_shell' "$overexp" "$inlist"
        change 'overshoot_f_above_burn_h_core' "$overexp" "$inlist"
        change 'overshoot_f_above_burn_h_shell' "$overexp" "$inlist"
        change 'overshoot_f_above_burn_he_core' "$overexp" "$inlist"
        change 'overshoot_f_above_burn_he_shell' "$overexp" "$inlist"
    fi
    
    if (( $(echo "$underexp > 0" | bc -l) )); then
        change 'overshoot_f_below_nonburn_shell' "$underexp" "$inlist"
        change 'overshoot_f_below_burn_h_shell' "$underexp" "$inlist"
        change 'overshoot_f_below_burn_he_shell' "$underexp" "$inlist"
    fi
    
    if (( $(echo "$f0 > 0" | bc -l) )); then
        change 'overshoot_f0_above_nonburn_core' "$f0" "$inlist"
        change 'overshoot_f0_above_nonburn_shell' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_h_core' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_h_shell' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_he_core' "$f0" "$inlist"
        change 'overshoot_f0_above_burn_he_shell' "$f0" "$inlist"
        change 'overshoot_f0_below_nonburn_shell' "$f0" "$inlist"
        change 'overshoot_f0_below_burn_h_shell' "$f0" "$inlist"
        change 'overshoot_f0_below_burn_he_shell' "$f0" "$inlist"
    fi
}

fix_mod() { # final_mod
    mod=$1
    check_mod "$mod"
    # MESA has a bug that makes it print numbers of the form *.***### 
    # (i.e. no letter for the exponent, and all the numbers are asterisks) 
    # This sed command replaces those with zeros. 
    if [[ $proceed -eq 1 ]]; then
        sed -i.bak "s/\*\.\**-[0-9]*/0.0000000000000000D+00/g" $mod
    fi
}

check_mod() {
    mod=$1
    # check that mod file got written
    if [[ ! -e $mod ]]; then 
        proceed=0
    fi
}

run() { # logs_dir  final_mod
    logs_dir=$1
    final_mod=$2
    
    if [[ $proceed -eq 0 ]]; then 
        return 0
    fi
    
    ./rn | tee output
    check_mod $final_mod
    num_logs=$(cat $logs_dir/history.data | wc -l)
    num_logs=$(echo "$num_logs - 7" | bc -l)
    while [[ $num_logs -lt $num_process && $proceed -eq 1 ]]; do
        
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
"if (nrow(DF) <= 1) { cat(0); quit('no'); } ;"\
"res <- (max(DF[['star_age']])-min(DF[['star_age']]))/"$num_process";"\
"res.str <- sub('\\\..+', '', paste(res));"\
"cat(sub('e\\\+0+|e0+', 'd', res.str))")
            if (( $(echo "$timestep < $min_years_for_timestep" | bc -l) )); then
                echo "Too small timesteps"
                proceed=0
            fi
            change 'max_years_for_timestep' "$timestep" $inlist
        fi
        
        rm -rf "$logs_dir"
        ./rn | tee output
        check_mod $final_mod
        num_logs=$(cat $logs_dir/history.data | wc -l)
        num_logs=$(echo "$num_logs - 7" | bc -l)
    done
    
    #if [[ $proceed -eq 1 ]]; then
    if [[ $proceed -eq 1 ]] && (( $(echo "$calibrate <= 0" | bc -l) )); then
        fix_mod "$final_mod"
        if [[ suppress -eq 0 ]] && [[ nopulse -eq 0 ]]; then
            if [[ "$inlist" = "inlist_3ms" ]] || [[ "$inlist" = "inlist_4sg" ]]
              then 
                Rscript $scriptdir/model_select.R $num_process $logs_dir | 
                    xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
                        "echo start {}; gyre2freqs.sh {}; echo end {}"
                        #"echo start {}; gyre-Dnu.sh {}; echo end {}"
            else
                Rscript $scriptdir/model_select.R $num_process $logs_dir | 
                    xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
                        "echo start {}; gyre-l0.sh {}; echo end {}"
                        #"echo start {}; gyre-Dnu.sh {}; echo end {}"
            fi
        fi
    fi
    
    if grep -q "termination code: max_age" output; then 
        proceed=0
    fi
}

################################################################################
### INITIALIZATION AND EVOLUTION ###############################################
################################################################################
## Make directory and copy over simulation files 
if [ -z ${expname+x} ]; then
    expname=M="$M"_Y="$Y"_Z="$Z"_alpha="$alpha"\
_diffusion="$diffusion"_settling="$settling"_eta="$eta"\
_overshoot="$overshoot"
fi
dirname=$directory/$expname
#if (( $(echo "$calibrate > 0" | bc -l) )); then
#    dirname=$directory
#fi

mkdir -p $dirname
cd $dirname
cp -r $scriptdir/mesa/* .
rm -rf LOGS/*
rm -f track
echo "id M Y Z alpha diffusion settling eta overshoot undershoot overexp underexp
$expname $M $Y $Z $alpha $diffusion $settling $eta $overshoot $undershoot $overexp $underexp" >> track
#echo M="$M"_Y="$Y"_Z="$Z"_alpha="$alpha"\
#_diffusion="$diffusion"_settling="$settling"_eta="$eta"\
#_overshoot="$overshoot"_undershoot="$undershoot"\
#_overexp="$overexp"_underexp="$underexp" >>> track

# pre-main sequence
set_params
./rn
fix_mod "pms.mod"

set_inlist "inlist_2pms"
set_overshoot
./rn
fix_mod "zams.mod"

# main sequence
set_inlist "inlist_3ms"
set_diffusion
if (( $(echo "$eta > 0" | bc -l) )); then
    change 'Reimers_scaling_factor' "$eta" "$inlist"
fi
#if (( $(echo "$calibrate > 0" | bc -l) )); then
if [ ! $calibrate = 0 ]; then
    change 'max_age' "$calibrate" "$inlist"
    change 'profile_interval' '99999' "$inlist"
    change 'history_interval' '99999' "$inlist"
    change 'max_num_profile_models' '1' "$inlist"
    num_process=1
fi
run "LOGS_MS" "tams.mod"

if [ ! $calibrate = 0 ]; then
    exit
fi

#if (( $(echo "$mainseq > 0" | bc -l) )); then
if [ ! $mainseq = 0 ]; then
    cleanup
fi

# sub-giant
set_inlist "inlist_4sg"
run "LOGS_SG" "brgb.mod"

#if (( $(echo "$subgiant > 0" | bc -l) )); then
if [ ! $subgiant = 0 ]; then 
    cleanup
fi

# red giant
set_inlist "inlist_5rgb"
run "LOGS_RGB" "bump.mod"

if [ ! $bump = 0 ]; then 
    cleanup
fi

if ! grep -q "stopping: phase_of_evolution >= x_integer_ctrl(1)" output; then 
    # bump to tip 
    set_inlist "inlist_6bump"
    run "LOGS_BUMP" "flash.mod"
    
    # helium burning 
    set_inlist "inlist_7heb"
    run "LOGS_HEB" "hexh.mod"

else
    # no bump: straight to helium burning 
    set_inlist "inlist_7heb"
    change "saved_model_name" "bump.mod" "$inlist"
    run "LOGS_HEB" "hexh.mod"
fi


# summarize and be done 
#cd ../..
#Rscript summarize.R "$dirname"
cleanup

# fin.

