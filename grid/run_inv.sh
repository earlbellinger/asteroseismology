python3 sobol_dispatcher.py -d ms-ov   -n -N 8192 -MS -i 64 -r 

-net "'basic.net'"

python3 sobol_dispatcher.py -d inv     -n -N 1024 -i 32 -r -MS    -o 0   0   -oe 0 0 -u 0 0 -ue 0 0 -D 0 0 -g 0 0 -e 0 0 
python3 sobol_dispatcher.py -d inv-ov  -n -N 1024 -i 32 -r -MS    -o 0.2 0.2 -oe 0 0 -u 0 0 -ue 0 0 -D 0 0 -g 0 0 -e 0 0 
python3 sobol_dispatcher.py -d inv-D   -n -N 1024 -i 32 -r -MS -t -o 0   0   -oe 0 0 -u 0 0 -ue 0 0 -D 1 1 -g 1 1 -e 0 0 
python3 sobol_dispatcher.py -d inv-Dov -n -N 1024 -i 32 -r -MS -t -o 0.2 0.2 -oe 0 0 -u 0 0 -ue 0 0 -D 1 1 -g 1 1 -e 0 0 

maybe_sub.sh -p 16 Rscript collate.R inv 
maybe_sub.sh -p 16 Rscript collate.R inv-ov 
maybe_sub.sh -p 16 Rscript collate.R inv-D 
maybe_sub.sh -p 16 Rscript collate.R inv-Dov 


maybe_sub.sh -e -p 8 python3 learn.py ../grid/inv.dat inv
maybe_sub.sh -e -p 8 python3 learn.py ../grid/inv-ov.dat inv-ov
maybe_sub.sh -e -p 8 python3 learn.py ../grid/inv-D.dat inv-D
maybe_sub.sh -e -p 8 python3 learn.py ../grid/inv-Dov.dat inv-Dov



#python3 sobol_dispatcher.py -d inv-ov  -n -s 21024 -N 7168 -i 32 -r -MS    -o 0.2 0.2 -oe 0 0 -u 0 0 -ue 0 0 -D 0 0 -g 0 0 -e 0 0 
#python3 sobol_dispatcher.py -d inv-D   -n -s 21024 -N 7168 -i 32 -r -MS -t -o 0   0   -oe 0 0 -u 0 0 -ue 0 0 -D 1 1 -g 1 1 -e 0 0 
#python3 sobol_dispatcher.py -d inv-Dov -n -s 21024 -N 7168 -i 32 -r -MS -t -o 0.2 0.2 -oe 0 0 -u 0 0 -ue 0 0 -D 1 1 -g 1 1 -e 0 0 


#maybe_sub.sh -p 16 ./dispatch.sh -d testseis -M 1.235342 -Y 0.292697
