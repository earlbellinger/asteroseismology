#!/bin/bash

#### Converter for FGONG to oscillation mode frequencies using ADIPLS 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

### Parse command line tokens 
## -h and --help will display help instructions 
if [ -z ${1+x} ] || [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "Converter for FGONG to oscillation mode frequencies using ADIPLS."
    echo "Usage: ./fgong2freqs.sh FGONG OUTPUT"
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

## Create a redistribute file and pass it to ADIPLS' redistrb 
# rerun with fewer points if the redistribution failed 
for nn in 9602 8602 7602 6602 5602 4602 3602 2602 1602; do
    echo "
2 '$fname.amdl'    @    
3 '$fname.model'    @
-1 ''        @
nn  ,icnmsh
$nn,      ,,  @
icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
211  ,      ,     ,0.013,      ,5.    ,      , 
11   ,      ,     ,   64,      ,5.    ,      , @
cg,  cx,ca  ,sig1,sig2,lmax,alphsf,adda4,accrm
1.,0.05,0.05,0.  ,0.  ,   2,      , 0.02, 0.01, @ 
nout,cn,irsu,unew   
60  ,  ,    ,    ,,,,,,, @
nmodel,kmodel,itsaml,ioldex
      ,      ,      ,     2,,,,,,,,  @
" > "redistrb-$fname.in"
    (egrep '# *$|@ *$' "redistrb-$fname.in" | sed -e 's/ *[#,@] *$//') | \
        $aprgdir/adiajobs/redistrb.c.d.x
    if [ -e "$fname.model" ]; then break; fi
    echo "Failed with nn = $nn, retrying with smaller nn"
    sleep 1
done

## Check that the redistribution was successful
if [ ! -e "$fname.model" ]; then
    echo "Error: Redistribution of $fname failed"
    exit 1
fi

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
    ,   4,   0,   1,     0,      ,    1,    1, @
  itrsig,sig1,istsig,inomde,itrds,
       1,   5,     1,     2,     , @
  dfsig,nsig,iscan,sig2,
       ,  2,  5000,2500, @
  eltrw1, eltrw2, sgtrw1, sgtrw2,
       0,     -1,      0,     -1,  @
cst:
  cgrav
  6.67408e-8  @
int:
  iplneq,iturpr,icow,alb,
        ,      ,    ,   , @
  istsbc,fctsbc,ibotbc,fcttbc,
        ,      ,      ,      , @
  mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre,
        ,     1, 0.9,      ,   ,      ,   16,      , @
  fsig,dsigmx,irsevn,xmnevn,nftmax,itsord,
      ,      ,      ,      ,      ,      , @
out:
  istdpr,nout,nprcen,irsord,iekinr,
       9,    ,      ,    20,      , @
  iper,ivarf,kvarf,npvarf,nfmode,
     1,    1,     ,      ,     0, @
  irotkr,nprtkr,igm1kr,npgmkr,ispcpr,
       0,      ,     0,      ,     0, @
  icaswn, sigwn1, sigwn2, frqwn1, frqwn2,iorwn1, iorwn2, frlwn1, frlwn2
        ,      0,     -1,      0,   9999,     0,     50,      0,     -1, @
dgn:
  itssol,idgtss,moddet,iprdet,npout,
        ,      ,      ,      ,     , @
  imstsl,imissl,imjssl,idgnrk,
        ,      ,      ,      , @
" > "adipls-$fname.in"

### Time to Run ADIPLS!
(egrep '# *$|@ *$' adipls-"$fname".in | sed -e 's/ *[#,@] *$//') | \
    $aprgdir/adipls/adipls.c.d.x

## Check that the .agsm file was created successfully, and convert it to freqs 
if [ ! -e "$fname.agsm" ]; then
    echo "Error: Failed to generate $fname.agsm"
    exit 1
fi
(echo 16; echo "$fname.agsm"; echo "$fname.dat"; echo "2") | \
	$aprgdir/adiajobs/set-obs.d.x

## Check that the frequencies were created, and remove the extra text 
if [ ! -e "$fname.dat" ]; then
    echo "Error: Failed to generate frequency information $fname.dat"
    exit 1
fi
cp "$fname.dat" "$fname.dat.bak"
cat "$fname.dat.bak" | cut -b 1-36 | \
    awk -v FIELDWIDTHS="5 7 10 13" -v OFS=, '{print $1,$2,$3,$4}' | \
    sed "s\,\ \g" | tr -s ' ' >| "$fname.dat"

### Hooray!
cp "$fname.dat" ..
echo "Conversion complete. Results can be found in $path/$fname.dat"
exit 0

