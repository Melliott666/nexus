#!/bin/bash

# Geant4 Path, edit G4Install path to where the main geant4 code folder your downloaded
export G4INSTALL=/Users/$USER/packages/geant4-v11/geant4-v11.3.2/install/;
export PATH=$G4INSTALL/bin:$PATH;
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$G4INSTALL/lib;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4INSTALL/lib;

cd $G4INSTALL/bin;
source geant4.sh;
cd -;

# Path to GSL
export GSL_PATH=$(brew --cellar gsl)/2.8;
export PATH=$GSL_PATH/bin:$PATH;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSL_PATH/lib;
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$GSL_PATH/lib;

# Path to HDF5
export HDF5_PATH=$(brew --cellar hdf5)/1.14.6/;
export HDF5_LIB=$(brew --cellar hdf5)/1.14.6/lib/;
export HDF5_INC=$(brew --cellar hdf5)/1.14.6/include/;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_LIB;
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$HDF5_LIB;

# Source the default nexus setup file to set the nexusdir environmental variable
source scripts/nexus_setup.sh

# Add the nexus exe to the path
export PATH=$PATH:/Users/$USER/packages/nexus/bin/
echo "Setup Nexus is complete!!"
