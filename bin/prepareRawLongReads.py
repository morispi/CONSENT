#!/usr/bin/env python3

import sys
import subprocess

f = open(sys.argv[1])

finalString = ""
id = f.readline()
while id != "":
	seq = f.readline()
	out = open(sys.argv[2] + id[1:-1], "w")
	out.write(id + seq.upper())
	out.close()
	id = f.readline()
f.close()
