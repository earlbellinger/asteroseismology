for directory in $(ls simulations); do
    if [ -d simulations/$directory/LOGS_MS ]; then
        cd simulations/$directory/LOGS_MS
        ls *.FGONG | 
            xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
                "echo start {}; 
                 python3 $SCRIPTS_DIR/fgong2ascii.py -i {};
                 echo end {}"
        find . -d -type d -exec rm -rf '{}' \;
        cd -
    fi
    
    if [ -d simulations/$directory/LOGS_SG ]; then
        cd simulations/$directory/LOGS_SG
        ls *.FGONG | 
            xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
                "echo start {}; 
                 python3 $SCRIPTS_DIR/fgong2ascii.py -i {};
                 echo end {}"
        find . -d -type d -exec rm -rf '{}' \;
        cd -
    fi
done
