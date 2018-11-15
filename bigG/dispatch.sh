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
Y=0.272804507715114
Z=0.0185654918074489
alpha=1.83454527401853
diffusion=1
age='4.572d9'
beta=0
HELP=0
remove=0
directory=simulations

while [ "$#" -gt 0 ]; do
  case "$1" in
    -h) HELP=1; break;;
    -r) remove=1; shift 1;;
    -n) expname="$2"; shift 2;;
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;
    -D) diffusion="$2"; shift 2;;
	-t) age="$2"; shift 2;;
	-b) beta="$2"; shift 2;;
    -d) directory="$2"; shift 2;;
    
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
    change 'x_ctrl(1)' "$beta" "$inlist"
    change 'max_age' "$age" "$inlist"
}

set_diffusion() {
    if (( $(echo "$diffusion > 0" | bc -l) )); then
        change 'do_element_diffusion' '.true.' "$inlist"
        change 'diffusion_SIG_factor' "$diffusion" "$inlist"
        change 'diffusion_GT_factor' "$diffusion" "$inlist" 
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

################################################################################
### INITIALIZATION AND EVOLUTION ###############################################
################################################################################
## Make directory and copy over simulation files 
if [ -z ${expname+x} ]; then
    expname=M="$M"_Y="$Y"_Z="$Z"_alpha="$alpha"\
_diffusion="$diffusion"_age="$age"_beta="$beta"
fi
dirname=$directory/$expname

mkdir -p $dirname
cd $dirname
cp -r $scriptdir/mesa/* .
rm -rf LOGS/*
rm -f track
echo "id M Y Z alpha diffusion age beta
$expname $M $Y $Z $alpha $diffusion $age $beta" >> track

# pre-main sequence
set_params
./rn
fix_mod "pms.mod"

set_inlist "inlist_2pms"
set_diffusion
./rn
fix_mod "zams.mod"

# main sequence
# set_inlist "inlist_3ms"
# ./rn

gyre2freqs.sh -i LOGS/profile1.data.GYRE -S

cd - 
Rscript summarize.R "$dirname"
if [ $remove -eq 1 ]; then
    rm -rf "$dirname"
fi
# fin.
