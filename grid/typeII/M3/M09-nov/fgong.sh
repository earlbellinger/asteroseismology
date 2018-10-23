cd LOGS_3MS
ls *.FGONG | xargs -L1 -I{} python3 $SCRIPTS_DIR/fgong2ascii.py -i {}
