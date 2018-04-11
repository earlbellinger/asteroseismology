#### Takes a directory name and calculates the u,Y kernels for the models there 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

for dir in $(ls models/$1); do
    #continue 
    cd models/$1/$dir/LOGS_MS
    rm -rf profile1-freqs
    maybe_sub.sh -p 1 kerexact.sh profile1.data.FGONG profile1-freqs 5
    cd -
    #exit 
done 

exit 
change() { 
    # Modifies a parameter in the current inlist. 
    # args: ($1) name of parameter 
    #       ($2) new value 
    #       ($3) filename of inlist where change should occur 
    # Additionally changes the 'inlist_0all' inlist. 
    # example command: change initial_mass 1.3 
    # example command: change log_directory 'LOGS_MS' 
    # example command: change do_element_diffusion .true. 
    param=$1 
    newval=$2 
    filename=$3 
    escapedParam=$(sed 's/[^^]/[&]/g; s/\^/\\^/g' <<< "$param")
    search="^\s*\!*\s*$escapedParam\s*=.+$" 
    replace="      $param = $newval" 
    sed -r -i.bak -e "s/$search/$replace/g" $filename 
    if [ ! "$filename" == 'inlist_0all' ]; then 
        change $param $newval 'inlist_0all'
    fi 
} 


for dir in $(ls models/$1); do
    cd models/$1/$dir
    change 'write_pulsation_plot_data' '.false.' 'inlist_0all'
    change 'add_center_point_to_pulse_info' '.false.' 'inlist_0all'
    change 'keep_surface_point_for_pulse_info' '.false.' 'inlist_0all'
    change 'add_double_points_to_pulse_info' '.false.' 'inlist_0all'
    export MESA_DIR=/scratch/seismo/bellinger/MESA/mesa-r8118
    maybe_sub.sh -p 2 ./rn 
    cd -
    #exit
    #exit 
done 


