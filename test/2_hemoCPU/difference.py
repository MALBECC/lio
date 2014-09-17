#!/usr/bin/python

import fileinput
import re
import collections

SEC = 1000000
times = collections.defaultdict(list)

def avg(l):
    return sum(l) / len(l)

for line in fileinput.input():
    parts = re.match(r'(Sphere|Cube) piece (\d+) took this much time: (?:(\d+)s\. )?(\d+)us\. and it has (\d+) elements \((\d+) nanounits\)', line)
    _, part, sec, usec, size, cost = parts.groups()
    time = int(sec or 0) * SEC + int(usec)
    times[int(part)].append(time)

avgs = dict([(k,avg(v)) for (k,v) in times.items()])
for (k,v) in avgs.items():
    print "Average for %d: %f" % (k, v)

vals = avgs.values()
diff = max(vals) - min(vals)
print "Max difference: %d" % diff
print "Percent: %f %%" % ((100.0 * diff) / (2.5 * SEC))
