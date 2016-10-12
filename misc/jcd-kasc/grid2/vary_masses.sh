#!/bin/bash
cp -R $MESA_DIR/star/work .
(cd work; ./mk)
mkdir mass
for i in $(seq 0.8 0.1 10); do 
    echo $i
    cp -R work mass/$i
    cp isochrone_template mass/$i/inlist_project
    cp run.sh mass/$i
    cd mass/$i
    sed -i.bak "s/initial_mass = 1/initial_mass = $i/g" inlist_project
    maybe_sub.sh -p 4 ./run.sh
    cd -
done
rm -rf work

