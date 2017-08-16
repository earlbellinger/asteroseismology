#!/bin/bash

#### A script to install GYRE 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung & Yale University 

#################
### Variables ###
#################
# Adjust this if you want a different version 
GYRE_VER=5.0

#########################################
### Download and install GYRE ###
#########################################
mkdir GYRE
cd GYRE

## Download and unpack GYRE
wget https://bitbucket.org/rhdtownsend/gyre/downloads/gyre-"$GYRE_VER".tar.gz 
mkdir gyre-"$GYRE_VER"
tar xvfz gyre-"$GYRE_VER".tar.gz  -C gyre-"$GYRE_VER"
ln -sfn gyre-"$GYRE_VER"/gyre gyre
export GYRE_DIR=$(pwd)/gyre
echo export GYRE_DIR="$GYRE_DIR" >> ~/.bash_profile

## Install GYRE
cd $GYRE_DIR
make

