#!/bin/bash

cd poaV2;
make poa;
cp poa ../bin/poa;
cd ..;
cd minimap;
make;
cd ..;
make;
