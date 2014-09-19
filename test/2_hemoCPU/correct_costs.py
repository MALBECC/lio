#!/usr/bin/python
import fileinput

THRESHOLD = 500000
for line in fileinput.input():
    v = int(line)
    if v < THRESHOLD:
        print THRESHOLD
    else:
        print v
