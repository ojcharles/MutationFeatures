#!/bin/bash
# installs software that are multi-command

mkdir /tools
cd /tools

##### dssp
git clone https://github.com/ojcharles/libcifpp.git
cd libcifpp
git checkout ojcharles-patch-1
mkdir build
cd build
cmake .. 
cmake --build . --config Release
cmake --install .

cd /tools
git clone https://github.com/PDB-REDO/dssp.git
cd dssp
mkdir build
cd build
cmake ..
cmake --build . --config Release
cmake --install .


##### p2rank
cd /tools
wget https://github.com/rdk/p2rank/releases/download/2.4/p2rank_2.4.tar.gz
gunzip p2rank_2.4.tar.gz
tar -xvf p2rank_2.4.tar

