gfortran -o kerexact.x kerexact.f vintk.d.f zero.f derivk.d.f nrk.d.f store.d.f leq.d.f ofiles.f rdaa5.d.f lir.d.f skpcom.f calci.f vinta.d.f

precision=2
lmin=0
lmax=3
lstep=1
kpair=5
truncate=-1
amdl="modelS.amdl"
amdl2="modelS.amdl2"
amde="modelS.amde"
ichkr1="ichkr1"
ichkr2="ichkr2"

./kerexact.x <<< "$precision

$lmin $lmax $lstep

$kpair

$truncate

10 $amdl

14 $amde

11 $ichkr1

12 $ichkr2

9 $amdl2

-1 -1"






