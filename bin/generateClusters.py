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
        if prevHeader == "" :
                finalString = curHeader + "\n" + t[5] + "\n"
        elif curHeader == prevHeader:
                finalString += t[5] + "\n"
        else:
                outF = open(sys.argv[2] + prevHeader, 'w')
                outF.write(finalString)
                outF.close()
                finalString = curHeader + "\n" + t[5] + "\n"
        prevHeader = curHeader
        line = f.readline()
outF = open(sys.argv[2] + prevHeader, 'w')
outF.write(finalString)
outF.close()
f.close()

