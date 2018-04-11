KICs='12258514
6225718
10068307
6116048
6106415
3632418
7510397
8938364
8760414
7970740
7940546
10162436
8379927
8228742
8694723
9414417
4914923
12317678
10454113
9139163
7680114
10516096
10963065
12009504
9098294
5773345
8394589
6679371
7103006
5774694
5184732
6106415
8006161
12069449
12069424'

for KIC in $KICs; do
    Rscript make_models.R $KIC learn/covs-simulations/legacy2 legacy
done

./check_cal.sh

for KIC in $KICs; do
    echo $KIC
    for model in $(ls maybe_sub_logs | grep calibrate | grep $KIC); do
        echo $model
        result=$(tail -n 1 maybe_sub_logs/$model/condor.out | cut -c 5-)
        echo $result
    done
done

## rerun models that didn't work 
for KIC in $KICs; do
    echo $KIC
    for model in $(ls maybe_sub_logs | grep calibrate | grep $KIC); do
        result=$(tail -n 1 maybe_sub_logs/$model/condor.out | cut -c 5-)
        #echo $result
        if [[ $result > 0.00001 ]]; then
            cd maybe_sub_logs/$model
            
            job=$(cat condor.job | grep "Executable" | cut -c 15-)
            #echo $job
            cmd=$(cat $job | grep Rscript)
            read -ra cmdarr <<<"$cmd"
            Y0=${cmdarr[7]}
            alpha0=${cmdarr[11]}
            #echo $Y0
            #echo $alpha0
            
            line=$(tail -n 10 condor.out | grep "Nelder-Mead")
            read -ra arr <<<"$line"
            #echo $model
            Y=${arr[1]}
            alpha=${arr[2]}
            if [[ $alpha < 1.5 ]]; then alpha=1.5; fi
            if [[ $Y < 0.15 ]]; then Y=0.15; fi
            if [[ $Y > 0.45 ]]; then Y=0.45; fi
            if [[ $alpha > 3 ]]; then alpha=3; fi
            #echo $Y
            #echo $alpha
            
            #sed -r -i.bak -e "s/$Y0/$Y/g" $job
            #sed -r -i.bak -e "s/$alpha0/$alpha/g" $job
            #condor_submit condor.job
            cd -
        fi
    done
done

re='^[0-9]+([.][0-9]+)?$'
for KIC in $KICs; do
    convs=0
    for model in $(ls maybe_sub_logs | grep $KIC | grep calibrate); do
        conv=$(tail -n 1 maybe_sub_logs/$model/condor.out | cut -c 5-)
        #echo $conv
        if ! [[ $conv =~ $re ]] ; then
            conv=1
        fi
        convs=$(echo "$convs + $conv" | bc -l)
    done
    if [[ $convs < 0.00001 ]]; then
        echo $KIC
    fi
done 

KICs='6225718
6116048
8938364
7970740
8379927
8694723
4914923
10454113
9139163
7680114
10516096
10963065
9098294
8394589
5774694
5184732
8006161
12069449
12069424'

for KIC in $KICs; do
    ./make_kernels.sh $KIC
done

for KIC in $KICs; do
    Rscript run_optimize.R $KIC learn/covs-simulations/legacy2 legacy
done

