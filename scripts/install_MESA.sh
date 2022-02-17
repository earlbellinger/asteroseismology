#!/bin/bash

#### A script to install the MESA SDK and MESA to a 64-bit Linux system
#### Author: Earl Patrick Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre, Aarhus University, Denmark 

#################
### Variables ###
#################
# Adjust these if you want a different version 
MESA_VER=21.12.1
SDK_VER=20.12.1

# Adjust this if you want to use more than one processor
export OMP_NUM_THREADS=1 

######################
### Pre-requisites ###
######################
# this command is needed for my lubuntu laptop 
#sudo apt install libc-dev binutils make perl libx11-dev tcsh zlib1g-dev libncurses5-dev 

#########################################
### Download and install the MESA SDK ###
#########################################
mkdir MESA
cd MESA

## Download SDK
SDK_REV=mesasdk-x86_64-linux-"$SDK_VER"
curl --remote-name http://www.astro.wisc.edu/~townsend/resource/download/mesasdk/"$SDK_REV".tar.gz
mkdir $SDK_REV
tar xvfz "$SDK_REV".tar.gz -C "$SDK_REV"
ln -sfn "$SDK_REV"/mesasdk mesasdk
export MESASDK_ROOT=$(pwd)/mesasdk
echo export MESASDK_ROOT="$MESASDK_ROOT" >> ~/.bashrc 

## Run SDK
source $MESASDK_ROOT/bin/mesasdk_init.sh
echo source "$MESASDK_ROOT"/bin/mesasdk_init.sh >> ~/.bashrc 

#################################
### Download and install MESA ###
#################################
## Download MESA
MESA_REV=mesa-r"$MESA_VER"
wget https://zenodo.org/record/5798242/files/"$MESA_REV".zip
unzip "$MESA_REV".zip
ln -sfn $(pwd)/"$MESA_REV" mesa
export MESA_DIR=$(pwd)/mesa
echo export MESA_DIR="$MESA_DIR" >> ~/.bashrc
echo export OMP_NUM_THREADS="$OMP_NUM_THREADS" >> ~/.bashrc

## Install MESA
cd $MESA_DIR
./install

