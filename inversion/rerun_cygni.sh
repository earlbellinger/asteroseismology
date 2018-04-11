cd models/perturb #models/array
for ii in $(ls); do
    echo $ii
    cd $ii/LOGS_MS
    rm -rf profile1-freqs 
    rm -rf profile1-freqs.dat 
    maybe_sub.sh -p 1 kerexact.sh profile1.data.FGONG 
    cd -
done

