#!/usr/bin/python

import fileinput
import re
import collections

SEC = 1000000
times = collections.defaultdict(list)

def avg(l):
    return sum(l) / len(l)

for line in fileinput.input():
    parts = re.match(r'Workload (\d+) took (\d+)s (\d+)ms and it has (\d+) elements \((\d+) nanounits\)', line)
    part, sec, usec, size, cost = parts.groups()
    time = int(sec or 0) * SEC + int(usec)
    times[int(part)].append(time)

avgs = dict([(k,avg(v)) for (k,v) in times.items()])
for (_,v) in avgs.items():
    print "%f" % v

vals = avgs.values()
mt = max(vals)
diff = max(vals) - min(vals)

print "Max difference: %d" % diff
print "Time for worst thread: %ds %dms" % (mt / SEC, mt % SEC)
print "Percent: %f %%" % ((100.0 * diff) / (2.5 * SEC))
