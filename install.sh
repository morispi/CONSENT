#!/bin/bash

cd poaV2;
make poa;
cd ..;
cd minimap2;
make;
cd ..;
make;
