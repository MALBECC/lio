#!/usr/bin/env python

import fileinput
import collections

sizes = collections.defaultdict(int)
for line in fileinput.input():
    u,_,v = map(int, line.split(" "))
    sizes[u] += v

largest = max(sizes.values())

for (k,v) in sizes.items(): 
    print "Thread %d takes %f %% of the largest thread (%d nanounits)" % (k, 100.0 * v / largest, v)
