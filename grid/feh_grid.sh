#-0.3
#-0.15
#0
#0.2

for M in `seq 0.5 0.1 3`; do 
    maybe_sub.sh -n -e -p 1 ./dispatch.sh -B -s -t -d feh_grid0   -M $M 
    maybe_sub.sh -n -e -p 1 ./dispatch.sh -B -s -t -d feh_grid-15 -M $M -Y 0.26310468 -Z 0.0117711
    maybe_sub.sh -n -e -p 1 ./dispatch.sh -B -s -t -d feh_grid-3  -M $M -Y 0.25833000 -Z 0.0084266
    maybe_sub.sh -n -e -p 1 ./dispatch.sh -B -s -t -d feh_grid02  -M $M -Y 0.28223339 -Z 0.0251701
done


for M in `seq 0.5 0.1 5`; do 
    maybe_sub.sh -n -e -p 1 ./dispatch.sh -MS -t -d feh_grid2_0   -M $M 
    maybe_sub.sh -n -e -p 1 ./dispatch.sh -MS -t -d feh_grid2_-15 -M $M -Y 0.26310468 -Z 0.0117711
    maybe_sub.sh -n -e -p 1 ./dispatch.sh -MS -t -d feh_grid2_-3  -M $M -Y 0.25833000 -Z 0.0084266
    maybe_sub.sh -n -e -p 1 ./dispatch.sh -MS -t -d feh_grid2_02  -M $M -Y 0.28223339 -Z 0.0251701
done

