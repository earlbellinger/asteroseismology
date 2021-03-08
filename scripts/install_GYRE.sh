#!/bin/bash

#### A script to install GYRE 
#### Author: Earl Patrick Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre, Aarhus University, Denmark 

#################
### Variables ###
#################
# Adjust this if you want a different version 
GYRE_VER=6.0.1

#################################
### Download and install GYRE ###
#################################
mkdir GYRE
cd GYRE

## Download and unpack GYRE
wget https://github.com/rhdtownsend/gyre/archive/v"$GYRE_VER".tar.gz
tar xvfz v"$GYRE_VER".tar.gz
ln -sfn gyre-"$GYRE_VER" gyre
export GYRE_DIR=$(pwd)/gyre
echo export GYRE_DIR="$GYRE_DIR" >> ~/.bashrc

## Install GYRE
cd $GYRE_DIR
make
