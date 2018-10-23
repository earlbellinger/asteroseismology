cd LOGS_3MS
ls | xargs -i --max-procs=$OMP_NUM_THREADS bash -c \
     "echo start {}; gyre-classical.sh -i {} -f; echo end {}"
