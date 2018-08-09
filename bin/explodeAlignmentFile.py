#!/usr/bin/env python3

import sys
import subprocess

f = open(sys.argv[1])
prevHeader = ""
finalString = ""


line = f.readline()
while line != '':
	t = line.split("\t");
	curHeader = t[0]
	if prevHeader == "" or curHeader == prevHeader:
		finalString += line
	else:
		outF = open(sys.argv[2] + prevHeader, 'w')
		outF.write(finalString)
		outF.close()
		finalString = ""
	prevHeader = curHeader
	line = f.readline()
outF = open(sys.argv[2] + prevHeader, 'w')
outF.write(finalString)
outF.close
f.close()