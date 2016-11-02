#python3 eig_extractor.py -i modelS/5b/modelS.amde -o modelS_eig -f
#python3 ker_extractor.py -f modelS/modelS.freq -g modelS/modelS.fgong -i modelS_eig -o modelS_ker
python3 eig_extractor.py -i modelS/modelS/ModelS-freqs/ModelS-freqs.amde -o modelS_eig -f -p
python3 ker_extractor2.py -f modelS/modelS/ModelS-freqs.dat -g modelS/modelS/ModelS.FGONG -i modelS_eig -o modelS_ker -p


python3 eig_extractor.py -i modelS/centre/ModelS_centre-freqs/ModelS_centre-freqs.amde -o modelS_centre_eig -f -p
python3 ker_extractor2.py -f modelS/centre/ModelS_centre-freqs.dat -g modelS/centre/ModelS_centre.FGONG -i modelS_centre_eig -o modelS_centre_ker
