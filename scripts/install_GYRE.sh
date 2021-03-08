#!/bin/bash

#### A script to install GYRE 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung & Yale University 

#################
### Variables ###
#################
# Adjust this if you want a different version 
GYRE_VER=6.0.1

#########################################
### Download and install GYRE ###
#########################################
mkdir GYRE
cd GYRE

## Download and unpack GYRE
wget https://github.com/rhdtownsend/gyre/archive/v"$GYRE_VER".tar.gz
tar xvfz v"$GYRE_VER".tar.gz
ln -sfn gyre-"$GYRE_VER" gyre
export GYRE_DIR=$(pwd)/gyre
echo export GYRE_DIR="$GYRE_DIR" >> ~/.bash_profile

## Install GYRE
cd $GYRE_DIR
make
