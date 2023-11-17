#!/bin/bash
# installs software that are multi-command

mkdir /tools
cd /tools

##### dssp
git clone https://github.com/PDB-REDO/libcifpp
cd libcifpp
mkdir build
cd build
cmake .. 
cmake --build . --config Release
cmake --install .


cd /tools
git clone https://github.com/mhekkel/libmcfp.git
cd libmcfp
mkdir build
cd build
cmake ..
cmake --build .
cmake --install .


cd /tools
git clone https://github.com/PDB-REDO/dssp.git
cd dssp
mkdir build
cd build
cmake ..
cmake --build . --config Release
cmake --install .


