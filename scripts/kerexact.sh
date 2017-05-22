#!/bin/bash

#### Obtain kernel functions for a given FGONG file 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

### Parse command line tokens 
## -h and --help will display help instructions 
if [ -z ${1+x} ] || [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "Obtain kernel functions for a given FGONG file."
    echo "Usage: ./kerexact.sh FGONG OUTPUT KERTYPE NSEL FRQWN IORWN"
    echo
    echo "Will make a directory called OUTPUT and place output files in there."
    echo "In absence of OUTPUT, will be created from the FGONG filename."
    echo "Outputs spherical degrees l=NSEL2..(NSEL1-1) (default NSEL=4,0)."
    echo "(or something like that)"
    echo "(that is stupid but it's not my fault)"
    echo "Searches through frequencies FRQWN1,FRQWN2 (default FRQWN=0,9999)"
    echo "and radial orders IORWN1,IORWN2 (default IORWN=0,50)"
    echo
    echo ' Kernel types'
    echo '   1.  Eulerian ln(c^2) and ln(rho)'
    echo '   2.  Eulerian ln(Gamma1) and ln(rho)'
    echo '   3.  Eulerian ln(c^2) and ln(Gamma1)'
    echo '   4.  Eulerian ln(u) and ln(Gamma1)'
    echo '   5.  Eulerian ln(u) and Y'
    echo '   6.  Eulerian ln(rho) and Y'
    echo '   7.  Eulerian ln(c) and ln(Gamma1/c)'
    echo '  11.  Lagrangian ln(c^2) and ln(rho)'
    echo '  12.  Lagrangian ln(Gamma1) and ln(rho)'
    echo '  13.  Lagrangian ln(c^2) and ln(Gamma1)'
    echo '  14.  Lagrangian ln(u) and ln(Gamma1)'
    echo '  17.  Lagrangian ln(c) and ln(Gamma1/c)'
    exit
fi

## Check that the first input (FGONG file) exists
if [ ! -e "$1" ]; then
    echo "Error: Cannot locate FGONG file $1"
    exit 1
fi
FGONG=$1

## Pull out the name of the FGONG file
bname="$(basename $FGONG)"
fname="${bname%%.*}-freqs"

## If the second (OUTPUT) argument doesn't exist, create one from the first
if [ -z ${2+x} ]; then
    path=$(dirname "$FGONG")/"$fname"
  else
    path="$2"
fi

## If the third (KERTYPE) argument doesn't exist, use "5" for u-Y kernels 
if [ -z ${3+x} ]; then
    kpair=5 # default to u-Y kernel
  else
    kpair="$3"
fi

## If the fourth (NSEL) argument doesn't exist, set it to 4
if [ -z ${4+x} ]; then
    nsel0=4
    nsel1=0
  else
    nsel=(${4//,/ })
    nsel0=${nsel[0]}
    nsel1=${nsel[1]}
    #nsel="$4"
fi

lmin=$nsel1
lmax=$(($nsel1 + $nsel0 - 1))
echo $lmin
echo $lmax

precision=2
lstep=1
truncate=-1 
amdl="$fname".model # "$fname".amdl # 
dgam=gamder.out
amde="$fname".amde
ichkr1=ichkr-"$kpair".1
ichkr2=ichkr-"$kpair".2
ichpsi='ichpsi'

## We need AMDL and AMDE files. Check if ADIPLS has already gotten them for us.
if [ ! -f "$path"/$amdl -o ! -f "$path"/$amde ]; then
    echo "Calling ADIPLS on $FGONG"
    fgong2freqs.sh $FGONG $path $4 $5 $6
fi

cd $path

## Convert FGONG to ASCII (for outside convenience)
if [ ! -f "$FGONG".dat ]; then
    python3 $SCRIPTS_DIR/fgong2ascii.py -i $FGONG
fi

#logfile="kerexact.log"
#exec > $logfile 2>&1

## If we are calculating Y kernels then we need EOS variables 
if [ $kpair -eq 5 ] || [ $kpair -eq 6 ]; then
    if [ ! -f $dgam ]; then
        echo "Extracting Gamma derivatives from $FGONG.dat"
        gamder.sh "$FGONG".dat
        #python3 $SCRIPTS_DIR/kern/dgam.py -i $FGONG -o "$fname".dgam
    fi
fi

echo "$precision

$lmin $lmax $lstep

$kpair

$truncate

10 $amdl

14 $amde

11 $ichkr1

12 $ichkr2

13 $ichpsi

9 $dgam

-1 -1" >| kerexact.in 

## Call kerexact with these settings 
echo "Calling kerexact"
cat kerexact.in | $SCRIPTS_DIR/kern/kerexact.x >> kerexact.log

## Now we have the kernels in binary, transform them to ASCII
echo "Transforming binary kernels to ASCII"
python3 $SCRIPTS_DIR/kern/ker_extractor.py -i $ichkr1 -o ichkr-"$kpair"_1/ -p >> kerexact.log
python3 $SCRIPTS_DIR/kern/ker_extractor.py -i $ichkr2 -o ichkr-"$kpair"_2/ -p >> kerexact.log

case $kpair in
    1|3|11|13) kpair1='c2';;
         2|12) kpair1='Gamma1';;
       4|5|14) kpair1='u';;
            6) kpair1='rho';;
         7|17) kpair1='c';;
esac
case $kpair in
    1|2|11|12) kpair2='rho';;
    3|4|13|14) kpair2='Gamma1';;
          5|6) kpair2='Y';;
         7|17) kpair2='Gamma1_over_c';;
esac
case $kpair in
    [1-7]) ktype='E_K';;
    [11-17]) ktype='L_K';;
esac

# Finally collate them into two files and 
echo "Collating kernels into two files"
paste ichkr-"$kpair"_1/* >| "$ktype"_"$kpair1"-"$kpair2".dat
paste ichkr-"$kpair"_2/* >| "$ktype"_"$kpair2"-"$kpair1".dat
#Rscript $SCRIPTS_DIR/kern/ker_combine.R $kpair

#exit

echo "Cleaning up"
rm -rf $ichkr1 $ichkr2 ichkr-"$kpair"_1 ichkr-"$kpair"_2

if [ $kpair -eq 5 ]; then
    python3 $SCRIPTS_DIR/kern/psi_extractor.py -i $ichpsi -o psi -p >> kerexact.log
    paste psi/* >| psi.dat
fi
rm -rf psi 
rm -rf ichpsi 


#cd - 

echo "Calculations complete. Results can be found in ""$path"
exit 0

