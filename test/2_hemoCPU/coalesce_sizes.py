#!/usr/bin/env python

import fileinput
import collections

sizes = collections.defaultdict(list)
lines = [line for line in fileinput.input()]
for line in lines[1:]:
    u,_,v = map(int, line.split(" "))
    sizes[u].append(v)

largest = max(map(sum,sizes.values()))

for (k,v) in sizes.items(): 
    print "Thread %d takes %f %% of the largest thread (%d nanounits, %d items)" % (k, 100.0 * sum(v) / largest, sum(v), len(v))
