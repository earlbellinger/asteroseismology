
python3 sobol_dispatcher.py -p 1 -N 128 -d hon    -Y 0.272957 0.272957 -Z 0.018583 0.018583 -a 1.836759 1.836759 -o 0 0 -D 1 1 -g 1 1 -e 0 0 -S -t
python3 sobol_dispatcher.py -p 1 -N 128 -d hon_ce -Y 1.4345 1.4345 -e 0 0 -o 0 0 -S -t -c


python3 sobol_dispatcher.py -p 1 -N 1024 -i 128 -d hon_ce2 -Y 1.4345 1.4345 -e 0 0 -D 1 1 -g 1 1 -o 0 0 -oe 0 0 -u 0 0 -ue 0 0 -S -t -c


python3 sobol_dispatcher.py -p 1 -N 128 -d hannah -M 0.85 1.45 -Y 0.272957 0.272957 -Z 0.018583 0.018583 -a 1.836759 1.836759 -o 0 0 -D 1 1 -g 1 1 -e 0 0 -S -t -rotk



python3 sobol_dispatcher.py -p 1 -N 128 -d honM    -Y 0.272957 0.272957 -Z 0.018583 0.018583 -a 1.836759 1.836759 -o 0 0 -D 1 1 -g 1 1 -e 0 0 -S -t


python3 sobol_dispatcher.py -N 512  -p 1 -d hon-new-M  -Y 0.272957 0.272957 -Z 0.018583 0.018583 -a 1.836759 1.836759  -o 0 0 -oe 0 0 -u 0 0 -ue 0 0 -e 0 0 -D 1 1 -S -t -C -n
python3 sobol_dispatcher.py -N 2048 -p 1 -d hon-new-ce -Y 1.4276221707417 1.4276221707417                              -o 0 0 -oe 0 0 -u 0 0 -ue 0 0 -e 0 0 -D 1 1 -S -t -C -n -c 



python3 sobol_dispatcher.py -N 8 -d grid_SG_US_exp  -M 0.7 1.8 -o  0 0 -u  0 0 -e 0 0 -S -p 4
python3 sobol_dispatcher.py -N 8 -d grid_SG_US_step -M 0.7 1.8 -oe 0 0 -ue 0 0 -e 0 0 -S -p 4


python3 sobol_dispatcher.py -N 7168 -s 21024 -d SG_US_exp  -M 0.7 1.8 -o  0 0 -u  0 0 -e 0 0 -S -p 1 -i 32 -r -C -n 
python3 sobol_dispatcher.py -N 7168 -s 21024 -d SG_US_step -M 0.7 1.8 -oe 0 0 -ue 0 0 -e 0 0 -S -p 1 -i 32 -r -C -n


maybe_sub.sh -p 2 ./dispatch.sh -Y 0.264382896490601 -a 1.68690375735786 -Z 0.0164895937976894 -D 0 -g 0 -c 4572000000 -d calibrate -n no_diffusion_final -f
maybe_sub.sh -p 2 ./dispatch.sh -Y 0.272804507715001 -a 1.83454527402773 -Z 0.0185654918074233 -D 1 -g 1 -c 4572000000 -d calibrate -n diffusion_final -f

 
