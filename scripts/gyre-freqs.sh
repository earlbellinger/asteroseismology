#!/bin/bash

#### Converter for .GYRE file to oscillation mode frequencies with GYRE
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

### Parse command line tokens 
## -h and --help will display help instructions 
if [ -z ${1+x} ] || [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "Converter for .GYRE to oscillation mode frequencies."
    echo "Usage: ./gyre-l02.sh .GYRE OUTPUT"
    echo
    echo "Will make a directory called OUTPUT and place OUTPUT.dat in there."
    echo "In absence of OUTPUT, will be created from the GYRE filename."
    exit
fi

## Check that the first input (GYRE file) exists 
if [ ! -e "$1" ]; then
    echo "Error: Cannot locate file $1"
    exit 1
fi

## Pull out the name of the GYRE file 
bname="$(basename $1)"
fname="${bname%%.*}-freqs"
pname="${bname::-5}"


## If the second (OUTPUT) argument doesn't exist, create one from the first 
if [ -z ${2+x} ]; then
    path=$(dirname "$1")/"$fname"
  else
    path="$2"
fi

## Create a directory for the results and go there 
mkdir -p "$path" 
#cp "$1" "$path" 
cd "$path" 

logfile="gyre-l02.log"
#exec > $logfile 2>&1

## Create a gyre.in file to find the large frequency separation 
echo "
&model
    model_type = 'EVOL'
    file = '../$bname'
    file_format = 'MESA'
    ! file_format = 'FGONG'
    ! data_format = '(1P5E16.9,x)'
/

&constants
    G_GRAVITY = 6.67428d-8 ! 6.67408d-8 !
/

&mode
    l = 0          ! Harmonic degree
/

&mode
    l = 2          ! Harmonic degree
/

&osc
    outer_bound = 'JCD'
    variables_set = 'JCD'
    inertia_norm = 'BOTH'
    !reduce_order = .FALSE.
    nonadiabatic = .TRUE.
/

&num
    !ivp_solver = 'MAGNUS_GL2'
/

&scan
    grid_type = 'LINEAR'
    freq_units = 'ACOUSTIC_CUTOFF'
    freq_min = 0.01
    freq_max = 1
    n_freq = 1000
/

&grid
/

&shoot_grid
/

&shoot_grid
    op_type = 'RESAMP_DISPERSION'
    alpha_osc = 50
    alpha_exp = 10
/

&recon_grid
/

&ad_output
/

&nad_output
    summary_file = '$fname.dat'
    summary_file_format = 'TXT'
    summary_item_list = 'l,n_pg,n_p,n_g,freq,E_norm'
    freq_units = 'UHZ'
/

" >| "gyre.in"

## Run GYRE
timeout 3600 $GYRE_DIR/bin/gyre gyre.in &>gyre.out

RETVAL=$?
if [ $RETVAL -eq 124 ]; then
    echo "search failed"
    exit 1
fi


### Hooray!
cp "$fname.dat" ..
#rm -rf *
currdir=$(pwd)
cd ..
#rm -rf "$currdir"
echo "Conversion complete. Results can be found in $fname.dat"
exit 0
