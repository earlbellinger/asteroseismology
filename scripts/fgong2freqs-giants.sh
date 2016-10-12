#!/bin/bash

#### Converter for FGONG to oscillation mode frequencies using ADIPLS 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

### Parse command line tokens 
## -h and --help will display help instructions 
if [ -z ${1+x} ] || [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "Converter for FGONG to oscillation mode frequencies using ADIPLS."
    echo "Usage: ./fgong2freqs.sh FGONG OUTPUT NPOINTS"
    echo
    echo "Will make a directory called OUTPUT and place OUTPUT.dat in there."
    echo "In absence of OUTPUT, will be created from the FGONG filename."
    echo "NPOINTS defaults to 40000."
    exit
fi

## Check that the first input (FGONG file) exists
if [ ! -e "$1" ]; then
    echo "Error: Cannot locate FGONG file $1"
    exit 1
fi

## Figure out how many points to use 
if [ ! -e "$3" ]; then
    NPOINTS=40000
else 
    NPOINTS="$3"
fi

## Pull out the name of the FGONG file
bname="$(basename $1)"
fname="${bname%%.*}-freqs"

## If the second (OUTPUT) argument doesn't exist, create one from the first
if [ -z ${2+x} ]; then
    path=$(dirname "$1")/"$fname"
  else
    path="$2"
fi

## Create a directory for the results and go there
mkdir -p "$path"
cp "$1" "$path" 
cd "$path"

logfile="fgong2freqs.log"
exec > $logfile 2>&1

## Convert the FGONG file to AMDL format
(echo "$bname"; echo "$fname.amdl"; echo 6.67408; echo 1) | \
    $aprgdir/adiajobs/fgong-amdl.d.x

## Check that the amdl file was created 
if [ ! -e "$fname.amdl" ]
  then
    echo "Error: Conversion of $1 to .amdl format failed"
    exit 1
fi

## Function to create a redistribute file and pass it to ADIPLS' redistrb 
run_redistrb() {
        echo "
2 '$fname.amdl'    @    
3 '$fname.model'    @
-1 ''        @
nn  ,icnmsh
$NPOINTS,      ,,  @
icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
12   ,      ,     ,0.013,      ,5.    ,      , @
201   ,      ,        64,      ,5.    ,      , 
cg,  cx,ca  ,sig1,sig2,lmax,alphsf,adda4,accrm
1.,0.05,0.05,0.  ,0.  ,   2,      , 0.02, 0.01, @ 
nout,cn,irsu,unew   
60  ,  ,    ,    ,,,,,,, @
nmodel,kmodel,itsaml,ioldex
      ,      ,      ,      ,,,,,,,,  @
" > "redistrb-$fname.in"
    (egrep '# *$|@ *$' "redistrb-$fname.in" | sed -e 's/ *[#,@] *$//') | \
        $aprgdir/adiajobs/redistrb.c.d.x
}

## Function to run ADIPLS 
run_adipls() { 
    (egrep '# *$|@ *$' adipls-"$fname".in | sed -e 's/ *[#,@] *$//') | \
        $aprgdir/adipls/adipls.c.d.x

    ## Check that the .agsm file was created successfully
    if [ ! -e "$fname.agsm" ]; then
        echo "Error: Failed to generate $fname.agsm"
        exit 1
    fi
    ## Convert .agsm to .dat using set-obs
    (echo 16; echo "$fname.agsm"; echo "$fname.dat"; echo "2") | \
        $aprgdir/adiajobs/set-obs.d.x

    ## Check that the frequencies were created, and remove the extra text 
    if [ ! -e "$fname.dat" ]; then
        echo "Error: Failed to generate frequency information $fname.dat"
        exit 1
    fi
    cp "$fname.dat" "$fname.dat.bak"
}

## Create an adipls.in file with some decent (?) settings 
echo "
 2  '$fname.model' @
 9  '$fname.log'  @
 11 '$fname.agsm' @
 4  '$fname.amde' @
 13 '$fname.g1k'  @
 15 '$fname.ssm'  @
 -1 ''   @
  cntrd,
mod.osc.cst.int.out.dgn     @
mod:
  ifind,xmod,imlds,in,irname,nprmod,
       ,    ,     ,  ,      ,      ,  @
  ntrnct,ntrnsf,imdmod,
       ,       ,      , @
osc:
  el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2,
   0,   3,   0,   1,      ,      ,     ,     , @
  itrsig,sig1,istsig,inomde,itrds,
   1,4.000000,      ,      ,     , @
  dfsig,nsig,iscan,sig2,
      0,  10,   10,3000.000000, @
  eltrw1, eltrw2, sgtrw1, sgtrw2,
        ,       ,       ,       ,  @
cst:
  cgrav
  6.67408e-8  @
int:
  iplneq,iturpr,icow,alb,
        ,      ,    ,   , @
  istsbc,fctsbc,ibotbc,fcttbc,
        ,      ,      ,      , @
  mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre,
       5,     0, 0.1,      ,1.0d-9,   ,   15,      , @
  fsig,dsigmx,irsevn,xmnevn,nftmax,itsord,
      ,      ,      ,      ,      ,      , @
out:
  istdpr,nout,nprcen,irsord,iekinr,
       9,    ,      ,    20,     1, @
  iper,ivarf,kvarf,npvarf,nfmode,
     1,    1,    2,     0,     1, @
  irotkr,nprtkr,igm1kr,npgmkr,ispcpr,
       0,     0,     0,     1,     0, @
  icaswn, sigwn1, sigwn2, frqwn1, frqwn2,iorwn1, iorwn2, frlwn1, frlwn2
   10010,       ,       ,       ,       ,      ,       ,       ,       , @
dgn:
  itssol,idgtss,moddet,iprdet,npout,
        ,      ,      ,      ,     , @
  imstsl,imissl,imjssl,idgnrk,
        ,      ,      ,      , @
" > "adipls-$fname.in"

# rerun with fewer points if redistribution or adipls fails 
run_redistrb
if [ ! -e "$fname.model" ]; then 
    echo "Error: Redistribution of $fname failed"
    exit 1
fi
run_adipls
cp "$fname.dat" "$fname.dat.bak"

cat "$fname.dat.bak" | cut -b 1-36 | \
    awk -v FIELDWIDTHS="5 7 10 13" -v OFS=, '{print $1,$2,$3,$4}' | \
    sed "s/,/ /g" | tr -s ' ' >| "$fname.dat"

### Hooray!
cp "$fname.dat" ..
rm -rf *
currdir=$(pwd)
cd ..
rm -rf "$currdir"
echo "Conversion complete. Results can be found in $fname.dat"
exit 0

