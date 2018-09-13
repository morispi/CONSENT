#!/usr/bin/env python3

import sys
import subprocess

f = open(sys.argv[1])
ref = open(sys.argv[3], 'w')
reads = open(sys.argv[4], 'w')

# Process first line, which is the read to correct
line = f.readline()[:-1]
read = open(sys.argv[2] + line)
header = read.readline()
sequence = read.readline()
res = header + sequence
ref.write(res)
ref.close()
res = ""

# Process other lines
line = f.readline()[:-1]
while line != '':
        read = open(sys.argv[2] + line)
        header = read.readline()
        sequence = read.readline()
        res += header + sequence
        line = f.readline()[:-1]
f.close()
reads.write(res)
reads.close()
