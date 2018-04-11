maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_alpha -a 1.4367593860737
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_alpha -a 1.5367593860737
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_alpha -a 1.6367593860737
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_alpha -a 1.7367593860737
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_alpha -a 1.8367593860737
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_alpha -a 1.9367593860737
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_alpha -a 2.0367593860737
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_alpha -a 2.1367593860737
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_alpha -a 2.2367593860737

maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Y -Y 0.242957104887671
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Y -Y 0.250457104887671
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Y -Y 0.257957104887671
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Y -Y 0.265457104887671
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Y -Y 0.272957104887671
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Y -Y 0.280457104887671
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Y -Y 0.287957104887671
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Y -Y 0.295457104887671
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Y -Y 0.302957104887671

maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Z -Z 0.0145827894799524
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Z -Z 0.0155827894799524
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Z -Z 0.0165827894799524
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Z -Z 0.0175827894799524
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Z -Z 0.0185827894799524
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Z -Z 0.0195827894799524
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Z -Z 0.0205827894799524
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Z -Z 0.0215827894799524
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_Z -Z 0.0225827894799524

maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_mass -M 0.94
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_mass -M 0.955
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_mass -M 0.97
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_mass -M 0.985
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_mass -M 1
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_mass -M 1.015
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_mass -M 1.03
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_mass -M 1.045
maybe_sub.sh -e -p 1 ./dispatch.sh -B -s -t -d thesis_mass -M 1.06


for ii in $(ls); do cd $ii; ../../combine_hist.sh; cd -; done
for ii in $(ls); do cd $ii; cp history.data ../$ii.dat; cd -; done

