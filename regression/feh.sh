maybe_sub.sh -p 8 -m 55000M -i 75000000 python3 learn.py ../grid/SG_US_step.dat final 
maybe_sub.sh -p 8 -m 55000M -i 75000000 python3 learn.py ../grid/SG_US_step.dat feh 

for allname in $(ls perturb | grep feh_); do
    maybe_sub.sh -p 4 -m 55000M -i 75000000 python3 learn.py ../grid/SG_US_step.dat $allname 
done 

for allname in $(ls perturb | grep teff_); do
    maybe_sub.sh -p 4 -m 55000M -i 75000000 python3 learn.py ../grid/SG_US_step.dat $allname 
done 

for allname in $(ls perturb | grep both_); do
    maybe_sub.sh -p 4 -m 55000M -i 75000000 python3 learn.py ../grid/SG_US_step.dat $allname 
done 



for allname in $(ls -1t perturb | head -45); do
    maybe_sub.sh -p 4 -m 55000M -i 75000000 python3 learn.py ../grid/SG_US_step.dat $allname 
done 


