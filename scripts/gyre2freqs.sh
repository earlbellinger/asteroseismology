#!/bin/bash

#### Converter for .GYRE file to oscillation mode frequencies with GYRE
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Updated: March 2018

### Parse command line tokens 

HELP=0
EIGENF=0
SAVE=0
RADIAL=0
FGONG=0
OMP_NUM_THREADS=1

while [ "$#" -gt 0 ]; do
  case "$1" in
    -h) HELP=1; break;;
    -i) INPUT="$2"; shift 2;;
    -o) OUTPUT="$2"; shift 2;;
    -t) OMP_NUM_THREADS="$2"; shift 2;;
    -r) RADIAL=1; shift 1;;
    -e) EIGENF=1;SAVE=1; shift 1;;
    -f) FGONG=1; shift 1;;
    -s) SAVE=1; shift 1;;
    
    *) if [ -z "$INPUT" ]; then 
           INPUT="$1"
           shift 1
         else 
           echo "unknown option: $1" >&2
           exit 1
       fi;
  esac
done

if [ $HELP -gt 0 ] || [ -z "$INPUT" ]; then
    echo "Converter for .GYRE files to oscillation mode frequencies."
    echo "Usage: ./gyre2freqs.sh -i input -o output -t threads -e -r"
    echo "Flags: -s : save calculations directory"
    echo "       -e : calculate eigenfunctions (automatically turns on -s)"
    echo "       -r : only calculate radial modes"
    echo "       -f : FGONG file format"
    exit
fi

## Check that the first input (GYRE file) exists
if [ ! -e "$INPUT" ]; then
    echo "Error: Cannot locate GYRE file $INPUT"
    exit 1
fi

## Pull out the name of the GYRE file
bname="$(basename $INPUT)"
fname="${bname%%.*}-freqs"
pname="${bname::-5}"

## If the OUTPUT argument doesn't exist, create a path from the filename 
if [ -z ${OUTPUT+x} ]; then
    path=$(dirname "$INPUT")/"$fname"
  else
    path="$OUTPUT"
fi

MODES="
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
"
if [ $RADIAL -gt 0 ]; then
    MODES="
&mode
    l=0
/
"
fi

if [ $FGONG -gt 0 ]; then 
    FORMAT="'FGONG'
    data_format = '(1P5E16.9,x)'"
else
    FORMAT="'MESA'"
fi

MODE_ITEM_LIST=''
if [ $EIGENF -gt 0 ]; then
    MODE_ITEM_LIST="mode_file_format = 'TXT'
    mode_template = '%J-%L_%N'
    mode_item_list = 'M_star,R_star,l,n_pg,n_p,n_g,freq,E,E_p,E_g,E_norm,M_r,x,xi_r,xi_h'"
fi

## Create a directory for the results and go there
mkdir -p "$path" 
#cp "$INPUT" "$path" 
cd "$path" 

logfile="gyre-l0.log"
#exec > $logfile 2>&1

## Create a gyre.in file to find the large frequency separation
echo "&model
    model_type = 'EVOL'
    file = '../$bname'
    file_format = $FORMAT
/

&constants
    G_GRAVITY = 6.67408d-8
    M_SUN = 1.988475d33
    R_SUN = 6.957d10
    L_SUN = 3.828d33
!    G_GRAVITY = 6.672320000d-08
!    M_SUN = 1.989000000d+33
!    R_SUN = 6.959906258d+10
!    L_SUN = 3.845999350d+33
/

$MODES

&osc
    outer_bound = 'JCD'
    variables_set = 'JCD'
    inertia_norm = 'BOTH'
/

&num
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
/

&ad_output
    summary_file = '$fname.dat'
    summary_file_format = 'TXT'
    summary_item_list = 'l,n_pg,n_p,n_g,freq,E_norm'
    freq_units = 'UHZ'
    $MODE_ITEM_LIST
/

&nad_output
/

" >| "gyre.in"

## Run GYRE
$GYRE_DIR/bin/gyre gyre.in &>gyre.out

### Hooray!
cp "$fname.dat" ..
echo "Conversion complete. Results can be found in $fname.dat"
if [ $SAVE -gt 0 ]; then exit 0; fi
rm -rf *
currdir=$(pwd)
cd ..
rm -rf "$currdir"
exit 0

