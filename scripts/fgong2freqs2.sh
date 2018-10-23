#!/bin/bash

#### Converter for FGONG to oscillation mode frequencies using ADIPLS 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

cgrav=6.67408 #6.67428 #6.67408 #  #6.67232 #

### Parse command line tokens 
## -h and --help will display help instructions 
if [ -z ${1+x} ] || [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "Converter for FGONG to oscillation mode frequencies using ADIPLS."
    echo "Usage: ./fgong2freqs.sh FGONG OUTPUT NSEL FRQWN IORWN"
    echo
    echo "Will make a directory called OUTPUT and place OUTPUT.dat in there."
    echo "In absence of OUTPUT, will be created from the FGONG filename."
    echo "Outputs spherical degrees l=NSEL2..(NSEL1-1) (default NSEL=4,0)."
    echo "Searches through frequencies FRQWN1,FRQWN2 (default FRQWN=0,9999)"
    echo "and radial orders IORWN1,IORWN2 (default IORWN=0,50)"
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

## If the third (NSEL) argument doesn't exist, set it to 4
if [ -z ${3+x} ]; then
    nsel="4,0"
  else
    nsel="$3"
fi

## If the fourth (FRQWN) argument doesn't exist, set it to 0 and 9999
if [ -z ${4+x} ]; then
    frqwn="0,9999"
    #frqwn1=0
    #frqwn2=9999
  else
    frqwn=$4
    #frqwn=(${4//,/ })
    #frqwn1=${frqwn[0]}
    #frqwn2=${frqwn[1]}
fi

## If the fifth (IORWN) argument doesn't exist, set it to 0 and 50
if [ -z ${5+x} ]; then
    iorwn="0,50"
    #iorwn1=0
    #iorwn2=50
  else
    iorwn=$5
    #iorwn=(${5//,/ })
    #iorwn1=${iorwn[0]}
    #iorwn2=${iorwn[1]}
fi

#if [ -z {$6+x} ]; then
#    cgrav=$6
#fi

## Create a directory for the results and go there
mkdir -p "$path"
cp "$1" "$path" 
cd "$path"

logfile="fgong2freqs.log"
exec > $logfile 2>&1

## Convert the FGONG file to AMDL format
(echo "$bname"; echo "$fname.amdl"; echo $cgrav; echo 1) | \
    $aprgdir/adiajobs/fgong-amdl.d.x

## Check that the amdl file was created 
if [ ! -e "$fname.amdl" ]
  then
    echo "Error: Conversion of $1 to .amdl format failed"
    exit 1
fi

## Redistribute points
#if [ 1 -eq 0 ]; then
echo "
2 '$fname.amdl'    @
3 '$fname.model'   @
-1 ''        @
nn  ,icnmsh
4800,      ,,  
4000,      ,,  @
10000,      ,,  
icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
211  ,      ,     ,0.013,      ,5.    ,      , 
11   ,      ,     ,   64,      ,5.    ,      , old
24   ,      ,     ,     ,      ,      ,      , @  
23   ,      ,     ,     ,      ,      ,      , 
cg,  cx,ca  ,sig1,sig2,lmax,alphsf,adda4,accrm
  ,    ,    ,    ,    ,    ,      ,     ,     , @
1.,0.05,0.05,0.  ,0.  ,   2,      , 0.02, 0.01, old
,,,,,,,,,,,,,,,,, 
nout,cn,irsu,unew   
60  ,  ,    ,    ,,,,,,, old
30  ,  ,    ,    ,,,,,,, @
20,,,,,,,,,  
nmodel,kmodel,itsaml,ioldex
      ,      ,      ,     2,,,,,,,,  old
      ,      ,      ,      ,,,,,,,,  @
" > "redistrb-$fname.in"
(egrep '# *$|@ *$' "redistrb-$fname.in" | sed -e 's/ *[#,@] *$//') | \
    $aprgdir/adiajobs/redistrb.c.d.x
    
## Check that the redistribution was successful
if [ ! -e "$fname.model" ]; then
    echo "Error: Redistribution of $fname failed"
    exit 1
fi
#fi

## Create an adipls.in file with some decent (?) settings 
echo "
input model
 2  '$fname.amdl' 
 2  '$fname.model' @

outputs
 4  '$fname.amde' @
 9  '$fname.log'  @
 11 '$fname.agsm' @
 13 '$fname.g1k'  @
 15 '$fname.ssm'  @
 12 '$fname.rotk' @ 
 -1 ''   @

cntrd:
mod.osc.cst.int.out.dgn     @

mod:
  ifind,xmod,imlds,in,irname,nprmod,
       ,    ,     ,  ,      ,      ,  @
  xtrnct,ntrnsf,imdmod,
        ,      ,      , @

osc:
  el,   nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2,
    ,     "$nsel",    ,      ,      ,     ,     , @
  itrsig,sig1,istsig,inomde,itrds,
       1,   5,     1,     2,     , 
       1,   5,     1,      ,     , @
       1, 4.0,      ,     1,     , 
       1,   1,      ,     1,   10, 
  dfsig,nsig,iscan,sig2,
       ,   2, 5000,2500, @
       ,   2, 1000,2500, 
       ,   2,   90,2600, 
       ,    ,     ,    , 
  eltrw1, eltrw2, sgtrw1, sgtrw2,
       0,     -1,      0,     -1, @ 
        ,       ,       ,       , 

cst:
  cgrav
  "$cgrav"e-8 @
  6.67428e-8 

int:
  iplneq,iturpr,icow,alb,
        ,      ,    ,   , @
        ,     1,   0,   , 
        ,     0,   0,   , 
  istsbc,fctsbc,ibotbc,fcttbc,
        ,      ,      ,      , 
       1,      ,      ,      , 
       1,     1,      ,      , 
       0,     0,      ,      , @
  mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre,
       1,     1, 0.5,      ,   ,      ,  100,      , 
       3,     1,0.99,      ,   ,      ,  100,      , 
       5,      ,0.99,      ,   ,      , 1000,      , 
       3,      ,0.99,      ,   ,      , 1000,      , @
       3,     1,0.99,      ,   ,      , 1000,      , 
       5,     1,0.99,      ,   ,      , 1000,      , 
       3,     1, 0.5,      ,   ,      , 1000,      , 
       5,      , 0.5,      ,   ,      , 1000,      , 
       3,     1,0.99,      ,   ,      ,   15,      , 
       5,     0, 0.5,      ,   ,      ,    8,      , 
  fsig,dsigmx,irsevn,xmnevn,nftmax,itsord,
      ,      ,      ,      ,      ,      , 
      ,      ,      ,      ,      ,    -1, @

out:
  istdpr,nout,nprcen,irsord,iekinr,
       9,    ,      ,    20,     1, @
       9,  10,      ,    20,      , 
       9,    ,      ,      ,      , 
  iper,ivarf,kvarf,npvarf,nfmode,
     1,    1,     ,      ,     3, 
     1,    3,     ,      ,     3, @
     1,    2,     ,      ,     3, 
     1,    1,    1,      ,     3, 
     1,     ,     ,      ,     3, 
  irotkr,nprtkr,igm1kr,npgmkr,ispcpr,
       1,      ,     0,      ,     0, @
  icaswn, sigwn1, sigwn2, frqwn1, frqwn2,iorwn1, iorwn2, frlwn1, frlwn2
        ,      0,     -1,      0,   9999,     0,     50,      0,     -1, 
        ,      0,     -1,       "$frqwn",      "$iorwn",      0,     -1, @
   10010,       ,       ,       ,       , -5000,    100,       ,       , 

dgn:
  itssol,idgtss,moddet,iprdet,npout,
        ,      ,      ,      ,     , @
  imstsl,imissl,imjssl,idgnrk,
        ,      ,      ,      , @
" >| "adipls-$fname.in"

## Run ADIPLS
(egrep '# *$|@ *$' adipls-"$fname".in | sed -e 's/ *[#,@] *$//') | \
    $aprgdir/adipls/adipls.c.d.x

## Check that the .agsm file was created successfully
if [ ! -e "$fname.agsm" ]; then
    echo "Error: Failed to generate $fname.agsm"
    exit 1
fi
## Convert .agsm to .dat using set-obs
(echo 16; echo "$fname.agsm"; echo "$fname.dat"; echo "5") | \
    $aprgdir/adiajobs/set-obs.d.x

## Check that the frequencies were created, and remove the extra text 
if [ ! -e "$fname.dat" ]; then
    echo "Error: Failed to generate frequency information $fname.dat"
    exit 1
fi
cp "$fname.dat" "$fname.dat.bak"

## Remove the text so that the data can be parsed 
cat "$fname.dat.bak" | cut -b 1-38 | \
    awk -v FIELDWIDTHS="5 7 12 14" -v OFS=, '{print $1,$2,$3,$4}' | \
    sed "s/,/ /g" | tr -s ' ' >| "$fname.dat"


## Also calculate the variational frequencies of the model 
(echo 1; echo "$fname.agsm"; echo "$fname.var"; echo "5") | \
    $aprgdir/adiajobs/set-obs.d.x

## Check that the frequencies were created, and remove the extra text 
if [ ! -e "$fname.var" ]; then
    echo "Error: Failed to generate frequency information $fname.dat"
    exit 1
fi
cp "$fname.var" "$fname.var.bak"

## Remove the text so that the data can be parsed 
cat "$fname.var.bak" | cut -b 1-26 | \
    awk -v FIELDWIDTHS="5 7 12" -v OFS=, '{print $1,$2,$3,$4}' | \
    sed "s/,/ /g" | tr -s ' ' >| "$fname.var"


python3 $SCRIPTS_DIR/rotk_extractor.py -i $fname.rotk -o rotk -p
paste rotk/* >| rotk.dat
rm -rf "$fname.rotk" rotk/

#rm -rf "$fname.amde"

### Hooray!
cp "$fname.dat" ..
echo "Conversion complete. Results can be found in $path/$fname.dat"
exit 0

