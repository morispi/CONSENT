#!/usr/bin/env python3

import sys
import re
import csv

f = open(sys.argv[1])

line = f.readline()
while line != '':
	t = line.split("\t")
	print(t[5] + "\t" + t[6] + "\t" + t[7] + "\t" + t[8] + "\t" + t[4] + "\t" + t[0] + "\t" + t[1] + "\t" + t[2] + "\t" + t[3] + "\t" + t[9] + "\t" + t[10] + "\t" + t[11] + "\t" + t[12] + "\t" + t[13] + "\t" + t[14])
	line = f.readline()