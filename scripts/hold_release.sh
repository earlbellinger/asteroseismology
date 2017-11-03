#!/bin/bash
# export CONDOR_LOCATION=/opt/condor/default
# export CONDOR_CONFIG=/opt/condor/default/etc/condor_config

# ## remove old jobs 
# /opt/condor/default/bin/condor_q bellinger -run | # list all my jobs 
    # grep 4\+[0-9] | # list my jobs older than 4 days 
    # grep -o ^[0-9]*.0 | # grab their id numbers 
    # xargs -i bash -c "/opt/condor/default/bin/condor_rm {}" # remove them 

# sleep 30

# ## hold and release jobs if other people are in the queue 
# #if /opt/condor/default/bin/condor_q | # list all jobs 
# #        grep -v "bellinger" | # find the ones that aren't mine 
# #        grep -q " I "; then # find if any are inactive 
# #    /opt/condor/default/bin/condor_hold -all
# #    sleep 30
# #    /opt/condor/default/bin/condor_release -all
# #fi

