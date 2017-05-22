#!/bin/bash

#### Converter for .GYRE file to oscillation mode frequencies with GYRE
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

### Parse command line tokens 
## -h and --help will display help instructions 
if [ -z ${1+x} ] || [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "Obtain large frequency separation from .GYRE file using GYRE."
    echo "Usage: ./gyre-Dnu.sh .GYRE OUTPUT"
    echo
    echo "Will make a directory called OUTPUT and place OUTPUT.dat in there."
    echo "In absence of OUTPUT, will be created from the GYRE filename."
    exit
fi

## Check that the first input (GYRE file) exists
if [ ! -e "$1" ]; then
    echo "Error: Cannot locate GYRE file $1"
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

logfile="gyre-Dnu.log"
#exec > $logfile 2>&1

echo "
scaling_nu_max_Viani <- function(R, M, Teff, mu, 
        Teff_sun=5777.739, nu_max_sun=3090, mu_sun=1.260134) { 
    M * nu_max_sun * sqrt(mu / mu_sun) / ( R**2 * sqrt(Teff/Teff_sun) ) 
} 

scaling_nu_ac_Viani <- function(R, M, Teff, mu, 
        Teff_sun=5777.739, nu_ac_sun=5300, mu_sun=1.260134) { 
    M * nu_ac_sun * sqrt(mu / mu_sun) / ( R**2 * sqrt(Teff/Teff_sun) ) 
} 

DF.head <- read.table(file.path('..', '$pname'), header=1, skip=1, nrow=1) 
DF.body <- read.table(file.path('..', '$pname'), header=1, skip=5, nrow=1) 

#hstry <- read.table(file.path('..', 'history.data'), header=1, skip=5) 
#h.row <- hstry[hstry['model_number'] == as.numeric(DF.head['model_number']),] 
#nu_ac <- as.numeric(h.row['acoustic_cutoff']) 

nu_ac <- scaling_nu_ac_Viani( 
    DF.head[['photosphere_r']], 
    DF.head[['star_mass']], 
    DF.head[['Teff']], 
    DF.body[['mu']] 
) 

nu_max <- scaling_nu_max_Viani( 
    DF.head[['photosphere_r']], 
    DF.head[['star_mass']], 
    DF.head[['Teff']], 
    DF.body[['mu']] 
) 

cat(paste0('freq_min = ', max(0.1, nu_max-(nu_ac-nu_max)), ' 
    freq_max = ', nu_ac)) 
" >| "nu_max.R" 
freqsearch=$(Rscript nu_max.R) 

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
    ivp_solver = 'MAGNUS_GL2'
/

&scan
    grid_type = 'LINEAR'
    ! freq_units = 'ACOUSTIC_CUTOFF'
    ! freq_min = 0.01
    ! freq_max = 1
    freq_units = 'UHZ'
    $freqsearch
    n_freq = 1000
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

&output
    summary_file = 'Dnu.dat'
    summary_file_format = 'TXT'
    summary_item_list = 'l,n_g,freq'
    freq_units = 'UHZ'
/

" >| "Dnu.in"

## Run GYRE
OMP_NUM_THREADS=1
$GYRE_DIR/bin/gyre_ad Dnu.in &>Dnu.out

## Calculate the large frequency separation
echo "
scaling_nu_max <- function(R, M, Teff, Teff_sun=5777, nu_max_sun=3090) 
    M * nu_max_sun / ( R**2 * sqrt(Teff/Teff_sun) ) 

scaling_nu_max_Viani <- function(R, M, Teff, mu, 
        Teff_sun=5777.739, nu_max_sun=3090, mu_sun=1.260134) { 
    M * nu_max_sun * sqrt(mu / mu_sun) / ( R**2 * sqrt(Teff/Teff_sun) ) 
} 

DF.head <- read.table(file.path('..', '$pname'), header=1, skip=1, nrow=1) 
DF.body <- read.table(file.path('..', '$pname'), header=1, skip=5, nrow=1) 

nu_max <- scaling_nu_max_Viani( 
    DF.head[['photosphere_r']], 
    DF.head[['star_mass']], 
    DF.head[['Teff']], 
    DF.body[['mu']] 
) 

nu_max_classic <- scaling_nu_max( 
    DF.head[['photosphere_r']], 
    DF.head[['star_mass']], 
    DF.head[['Teff']] 
) 

nus <- read.table('Dnu.dat', header=1, skip=5)
nus <- nus[nus['l'] == 0 & nus['n_g'] == 0,][['Re.freq.']] 

Dnu <- weighted.mean( diff(nus), dnorm(nus[-1], nu_max, 
      (0.66*nu_max**0.88) / (2*sqrt(2*log(2))) )) 

Dnu_classic <- weighted.mean( diff(nus), dnorm(nus[-1], nu_max_classic, 
      (0.66*nu_max_classic**0.88) / (2*sqrt(2*log(2))) )) 

write.table(
    data.frame( nu_max=nu_max, nu_max_classic=nu_max_classic, 
                Dnu=Dnu, Dnu_classic=Dnu_classic ), 
    '$fname.dat', row.names=F, quote=F) 
" >| "Dnu.R" 
Rscript Dnu.R >> "$fname.dat" 

### Hooray!
cp "$fname.dat" ..
rm -rf *
currdir=$(pwd)
cd ..
rm -rf "$currdir"
echo "Conversion complete. Results can be found in $fname.dat"
exit 0

