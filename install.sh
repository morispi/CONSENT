#!/bin/bash

mkdir -p bin
cd BMEAN;
./install.sh;
cd ..;
cd minimap2;
make;
cd ..;
make;
