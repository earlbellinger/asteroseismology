#cd $1
ls *.FGONG | xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
    '#echo start {};
    a={};
    bname=$(basename $a); 
    fname="${bname%%.*}-freqs";
    if [ ! -d $fname ] || [ ! -f $fname/E_K_u-Y.dat ]; then
        echo $fname;
        maybe_sub.sh kerexact.sh {} $fname 5;
    fi; 
    if [ ! -d $fname ] || [ ! -f $fname/E_K_c2-rho.dat ]; then
        echo $fname;
        maybe_sub.sh kerexact.sh {} $fname 1;
    fi; 
    if [ ! -d $fname ] || [ ! -f $fname/E_K_Gamma1-rho.dat ]; then
        echo $fname;
        maybe_sub.sh kerexact.sh {} $fname 2;
    fi; 
    #echo end {};
    '



