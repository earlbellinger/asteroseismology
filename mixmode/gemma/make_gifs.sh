
mkdir pres 
mkdir pres/ms 
mkdir pres/hook 
mkdir pres/crossing 
mkdir pres/sg 
mkdir pres/brgb 

anim_dir=animate-u-Y-l.1_n.11

for ii in `seq 640 674`; do 
    cp $anim_dir/"000"$ii".png" pres/ms
done

cd pres/ms
convert -delay 20 -loop 0 ms/*.png ms.gif
convert -delay 20 -loop 0 ms.gif -delay 100 ms/000670.png ms_pause.gif
cd -

for ii in `seq 674 6`; do 
    cp $anim_dir/"000"$ii".png" pres/ms
done


