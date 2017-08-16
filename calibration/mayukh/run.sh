maybe_sub.sh -p 4 Rscript spectral_calibrate.R F3V 4.301 6900 1.434460 187200067

maybe_sub.sh -p 4 Rscript spectral_calibrate.R G2V 4.438 5800 1.0
maybe_sub.sh -p 4 Rscript spectral_calibrate.R K0V 4.609 4900 0.9
maybe_sub.sh -p 4 Rscript spectral_calibrate.R K5V 4.699 4350 0.7

maybe_sub.sh -p 4 Rscript spectral_calibrate.R M0V 4.826 3900 0.5

maybe_sub.sh -p 4 Rscript spectral_calibrate-alpha-ov.R M2V 4.826 3700 0.5 

(cd models/F3V/LOGS_MS; python3 $SCRIPTS_DIR/fgong2ascii.py -i profile1.data.FGONG)
(cd models/G2V/LOGS_MS; python3 $SCRIPTS_DIR/fgong2ascii.py -i profile1.data.FGONG)
(cd models/K0V/LOGS_MS; python3 $SCRIPTS_DIR/fgong2ascii.py -i profile1.data.FGONG)
(cd models/K5V/LOGS_MS; python3 $SCRIPTS_DIR/fgong2ascii.py -i profile1.data.FGONG)
(cd models/M0V/LOGS_MS; python3 $SCRIPTS_DIR/fgong2ascii.py -i profile1.data.FGONG)
(cd models/M2V/LOGS_MS; python3 $SCRIPTS_DIR/fgong2ascii.py -i profile1.data.FGONG)

Rscript atmosphere.R F3V 
Rscript atmosphere.R G2V 
Rscript atmosphere.R K0V 
Rscript atmosphere.R K5V 
Rscript atmosphere.R M0V 
Rscript atmosphere.R M2V 

Rscript kiel.R 

