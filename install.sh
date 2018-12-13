#!/bin/bash

cd BMEAN;
./install.sh;
cd ..;
cd minimap2;
make;
cd ..;
make;
