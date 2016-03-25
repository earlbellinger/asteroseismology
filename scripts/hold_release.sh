#!/bin/bash
export CONDOR_LOCATION=/opt/condor/default
export CONDOR_CONFIG=/opt/condor/default/etc/condor_config
if /opt/condor/default/bin/condor_q | grep -v "bellinger" | grep -q " I "; then
    /opt/condor/default/bin/condor_hold -all
    sleep 30
    /opt/condor/default/bin/condor_release -all
fi
