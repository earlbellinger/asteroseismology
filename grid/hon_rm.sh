cd hon_ce
for ii in $(ls); do
   if [[ ! -d $ii/LOGS_SG ]]; then continue; fi
   echo $ii
   sleep 1
   cd $ii/LOGS_SG
   for jj in $(ls); do 
       test -d $jj && echo $jj && rm -rf $jj
   done
   cd - 
done


for ii in $(ls | grep -v ".dat"); do
    cd $ii
    echo $ii
    rm -rf `ls | grep -v ".dat\|track\|.mod$\|inlist_0all$\|inlist_1relax$\|inlist_2pms$\|inlist_3ms$\|^inlist$\|LOGS_MS\|LOGS_SG"`
    if [ -d "LOGS_MS" ]; then
        (cd LOGS_MS && rm -f `ls | grep -v "history.data\|freqs.dat\|profiles.index\|track\|.mod$\|inlist_0all$\|inlist_1relax$\|inlist_2pms$\|inlist_3ms$\|^inlist$\|LOGS_MS\|LOGS_SG"`)
    fi
    if [ -d "LOGS_SG" ]; then
        (cd LOGS_SG && rm -f `ls | grep -v "history.data\|freqs.dat\|profiles.index\|track\|.mod$\|inlist_0all$\|inlist_1relax$\|inlist_2pms$\|inlist_3ms$\|^inlist$\|LOGS_MS\|LOGS_SG"`)
    fi
    cd -
    #sleep 1
done


