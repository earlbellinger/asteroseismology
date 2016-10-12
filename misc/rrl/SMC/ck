#!/bin/bash

if [ -f checks.md5 ]; then
    if [ -s checks.md5 ]; then
	output=`md5sum -c checks.md5`
	if [ $? -ne 0 ]; then
	    echo "$output"
	    exit 1
	fi
    else
	echo "Empty checks.md5 file"
	exit 1
    fi
else
    echo "Missing checks.md5 file"
    exit 1
fi
