#!/bin/bash

#### Converter for .GYRE file to oscillation mode frequencies with GYRE
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

### Parse command line tokens 
## -h and --help will display help instructions 
input=$1
output=$2
if [ -z ${1+x} ] || [ $input == '-h' ] || [ $input == '--help' ]; then
    echo "Converter for .GYRE to oscillation mode frequencies."
    echo "Usage: ./gyre-l0.sh .GYRE OUTPUT"
    echo
    echo "Will make a directory called OUTPUT and place OUTPUT.dat in there."
    echo "In absence of OUTPUT, will be created from the GYRE filename."
    exit
fi

## Check that the first input (GYRE file) exists
if [ ! -e "$input" ]; then
    echo "Error: Cannot locate GYRE file $input"
    exit 1
fi

## Pull out the name of the GYRE file
bname="$(basename $input)"
fname="${bname%%.*}-freqs"
pname="${bname::-5}"

## If the second (OUTPUT) argument doesn't exist, create one from the first
if [ -z ${2+x} ]; then
    path=$(dirname "$input")/"$fname"
  else
    path="$output"
fi

## Create a directory for the results and go there
mkdir -p "$path" 
#cp "$input" "$path" 
cd "$path" 

logfile="gyre-l0.log"
#exec > $logfile 2>&1

## Create a gyre.in file to find the large frequency separation
echo "&model
    model_type = 'EVOL'
    file = '../$bname'
    file_format = 'MESA'
/

&constants
    G_GRAVITY = 6.67408d-8
    M_SUN = 1.988475d33
    R_SUN = 6.957d10
    L_SUN = 3.828d33
/

&mode
    l=0
/

&mode
    l=1
/

&mode
    l=2
/

&mode
    l=3
/

&osc
    outer_bound = 'JCD'
    variables_set = 'JCD'
    inertia_norm = 'BOTH'
    !nonadiabatic = .true.
/

&num
    diff_scheme = 'MAGNUS_GL6'
/

&scan
    grid_type = 'LINEAR'
    freq_min_units = 'UHZ'
    freq_max_units = 'ACOUSTIC_CUTOFF'
    freq_min = 0.01
    freq_max = 1
    n_freq = 1000
/

&grid
    alpha_osc = 10
    alpha_exp = 10
/

&ad_output
    summary_file = '$fname.dat'
    summary_file_format = 'TXT'
    summary_item_list = 'l,n_pg,n_p,n_g,freq,E,E_norm'
    freq_units = 'UHZ'
/

&nad_output
/

" >| "l0.in"

## Run GYRE
$GYRE_DIR/bin/gyre l0.in &>l0.out

#exit

### Hooray!
cp "$fname.dat" ..
rm -rf *
currdir=$(pwd)
cd ..
rm -rf "$currdir"
echo "Conversion complete. Results can be found in $fname.dat"
exit 0

