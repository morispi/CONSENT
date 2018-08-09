#!/usr/bin/env python3

import sys
import csv

f = open(sys.argv[1])

line = f.readline()
while line != "":
	t = line.split("\t")
	qName = t[0]
	qLength = t[1]
	qStart = t[2]
	qEnd = t[3]
	strand = t[4]
	tName = t[5]
	tLength = t[6]
	tStart = t[7]
	tEnd = t[8]
	resMatches = t[9]
	alBlockLen = t[10]
	mapQual = t[11]

	print(qName + "\t" + qLength + "\t" + qStart + "\t" + qEnd + "\t" + strand + "\t" + tName + "\t" + tLength + "\t" + tStart + "\t" + tEnd + "\t" + resMatches + "\t" + alBlockLen + "\t" + mapQual)
	print(tName + "\t" + tLength + "\t" + tStart + "\t" + tEnd + "\t" + strand + "\t" + qName + "\t" + qLength + "\t" + qStart + "\t" + qEnd + "\t" + resMatches + "\t" + alBlockLen + "\t" + mapQual)

	line = f.readline()

f.close()
