#!/bin/bash

#### Converter for FGONG to oscillation mode frequencies using GYRE 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

### Parse command line tokens 
## -h and --help will display help instructions 
if [ -z ${1+x} ] || [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "Converter for .GYRE to oscillation mode frequencies using GYRE."
    echo "Usage: ./fgong2freqs-gyre.sh FGONG OUTPUT"
    echo
    echo "Will make a directory called OUTPUT and place OUTPUT.dat in there."
    echo "In absence of OUTPUT, will be created from the FGONG filename."
    exit
fi

## Check that the first input (FGONG file) exists
if [ ! -e "$1" ]; then
    echo "Error: Cannot locate FGONG file $1"
    exit 1
fi

## Pull out the name of the FGONG file
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

logfile="fgong2freqs-gyre.log"
#exec > $logfile 2>&1

OMP_NUM_THREADS=1

## Calculate nu_max from the corresponding profile file 
echo "scaling_nu_max <- function(R, M, Teff, Teff_sun=5777, nu_max_sun=3090)
    M * nu_max_sun / ( R**2 * sqrt(Teff/Teff_sun) )

DF <- read.table(file.path('..', '$pname'), header=1, skip=1, nrow=1)
with(DF, cat(scaling_nu_max(photosphere_r, star_mass, Teff)))
" >| "nu_max.R"
numax=$(Rscript nu_max.R)

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

&osc
    outer_bound = 'JCD'
    variables_set = 'JCD'
    inertia_norm = 'BOTH'
    !reduce_order = .FALSE.
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

&shoot_grid
    op_type = 'CREATE_CLONE'
/

&recon_grid
    op_type = 'CREATE_CLONE'
/

&shoot_grid
    op_type = 'RESAMP_DISPERSION'    ! Resample the grid based on the local dispersion relation
    alpha_osc = 5            ! At least alpha points per oscillatory wavelength
    alpha_exp = 3            ! At least alpha points per exponential 'wavelength'
/

&shoot_grid
    op_type = 'RESAMP_CENTER'    ! Resample the grid at the center
    n = 20                ! At least n points in the evanescent region
/

&recon_grid
    op_type = 'RESAMP_DISPERSION'    ! Resample the grid based on the local dispersion relation
    alpha_osc = 5            ! At least alpha points per oscillatory wavelength
    alpha_exp = 3            ! At least alpha point per exponential 'wavelength'
/

&ad_output
    summary_file = 'Dnu.dat'
    summary_file_format = 'TXT'
    summary_item_list = 'freq'
    freq_units = 'UHZ'
/

&nad_output
/

" >| "Dnu.in"

## Run GYRE
$GYRE_DIR/bin/gyre Dnu.in &>Dnu.out

## Calculate the large frequency separation
echo "nus <- read.table('Dnu.dat', header=1, skip=5)[['Re.freq.']]
cat(weighted.mean( diff(nus), dnorm(nus[-1], $numax, 
      (0.66*$numax**0.88) / (2*sqrt(2*log(2)))
    )))
" >| "Dnu.R"
Dnu=$(Rscript Dnu.R)

## Calculate nu_max-5*Dnu and nu_max+5*Dnu
lower=$(echo "$numax - 8 * $Dnu" | bc -l)
upper=$(echo "$numax + 8 * $Dnu" | bc -l)
above=$(echo "$lower > 3" | bc -l)
lower=$([ $above == 1 ] && echo "$lower" || echo "3")

echo "Searching between $lower and $upper"

## Create a gyre.in file to search from nu_max-5*Dnu to nu_max+5*Dnu
## for ONLY l=0 and l=2 
echo "
&model
    model_type = 'EVOL'
    file = '../$bname'
    file_format = 'MESA'
    ! file_format = 'FGONG'
    ! data_format = '(1P5E16.9,x)'
/

&constants
    G_GRAVITY = 6.67408d-8
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
/

&num
    !ivp_solver = 'MAGNUS_GL2'
/

&scan
    grid_type = 'INVERSE'
    freq_units = 'UHZ'
    freq_min = $lower
    freq_max = $upper
    n_freq = 1000
/

&shoot_grid
    op_type = 'CREATE_CLONE'
/

&recon_grid
    op_type = 'CREATE_CLONE'
/

&shoot_grid
    op_type = 'RESAMP_DISPERSION'    ! Resample the grid based on the local dispersion relation
    alpha_osc = 5            ! At least alpha points per oscillatory wavelength
    alpha_exp = 3            ! At least alpha points per exponential 'wavelength'
/

&shoot_grid
    op_type = 'RESAMP_CENTER'    ! Resample the grid at the center
    n = 20                ! At least n points in the evanescent region
/

&recon_grid
    op_type = 'RESAMP_DISPERSION'    ! Resample the grid based on the local dispersion relation
    alpha_osc = 5            ! At least alpha points per oscillatory wavelength
    alpha_exp = 3            ! At least alpha point per exponential 'wavelength'
/

&ad_output
    summary_file = 'Dnu.dat'
    summary_file_format = 'TXT'
    summary_item_list = 'freq'
    freq_units = 'UHZ'
/

&nad_output
/

" >| "l02.in"

## Run GYRE
timeout 3600 $GYRE_DIR/bin/gyre l02.in &>l02.out

RETVAL=$?
if [ $RETVAL -eq 124 ]; then
    echo "l=0,2 search failed"
    exit 1
fi

## Create a gyre.in file to search from nu_max-5*Dnu to nu_max+5*Dnu
## for ONLY l=1
echo "
&model
    model_type = 'EVOL'
    file = '../$bname'
    file_format = 'MESA'
    ! file_format = 'FGONG'
    ! data_format = '(1P5E16.9,x)'
/

&constants
    G_GRAVITY = 6.67408d-8
/

&mode
    l = 1          ! Harmonic degree
/

&osc
    outer_bound = 'JCD'
    variables_set = 'JCD'
    inertia_norm = 'BOTH'
    !reduce_order = .FALSE.
/

&num
    !ivp_solver = 'MAGNUS_GL2'
/

&scan
    grid_type = 'INVERSE'
    freq_units = 'UHZ'
    freq_min = $lower
    freq_max = $upper
    n_freq = 1000
/

&shoot_grid
    op_type = 'CREATE_CLONE'
/

&recon_grid
    op_type = 'CREATE_CLONE'
/

&shoot_grid
    op_type = 'RESAMP_DISPERSION'    ! Resample the grid based on the local dispersion relation
    alpha_osc = 5            ! At least alpha points per oscillatory wavelength
    alpha_exp = 3            ! At least alpha points per exponential 'wavelength'
/

&shoot_grid
    op_type = 'RESAMP_CENTER'    ! Resample the grid at the center
    n = 20                ! At least n points in the evanescent region
/

&recon_grid
    op_type = 'RESAMP_DISPERSION'    ! Resample the grid based on the local dispersion relation
    alpha_osc = 5            ! At least alpha points per oscillatory wavelength
    alpha_exp = 3            ! At least alpha point per exponential 'wavelength'
/

&ad_output
    summary_file = 'Dnu.dat'
    summary_file_format = 'TXT'
    summary_item_list = 'freq'
    freq_units = 'UHZ'
/

&nad_output
/

" >| "l1.in"

## Run GYRE
timeout 3600 $GYRE_DIR/bin/gyre l1.in &>l1.out

RETVAL=$?
if [ $RETVAL -eq 124 ]; then
    echo "l=1 search failed"
else
    tail -n +7 l1.dat >> "$fname".dat
fi

### Hooray!
exit 

cp "$fname.dat" ..
rm -rf *
currdir=$(pwd)
cd ..
rm -rf "$currdir"
echo "Conversion complete. Results can be found in $fname.dat"
exit 0

