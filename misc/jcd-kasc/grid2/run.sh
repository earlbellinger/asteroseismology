./rn
ls LOGS/*.FGONG | xargs -i --max-procs=1 bash -c \
    "echo start {}; fgong2freqs.sh {}; echo end {}"
