#!/usr/bin/env python

import os
import subprocess
import re

MSECS_IN_SEC = 1000000.0
def time2nanos(spec):
    groups = re.search("(?:(\d+)s\. )?(\d+)us\.", spec)
    secs = groups.group(1) or 1
    msecs = groups.group(2)
    return MSECS_IN_SEC * float(secs) + float(msecs)

progname = "correr-profile.sh"
testdirs = [subdir for subdir, dirs, files in os.walk(".") if progname in set(files)]

flags_and_options = {
    "OMP_NUM_THREADS":      [1,12,12],
    "MKL_NUM_THREADS":      [1,1,12],
    "LIO_INNER_THREADS":    [1,1,12],
    "LIO_OUTER_THREADS":    [1,12,1],
}

always_set = [("OMP_NESTED", "true"), ("MKL_DYNAMIC","false")]
enviros = [[(key, str(vals[i])) for key,vals in flags_and_options.items()] for i in xrange(0,3)] 

for directory in testdirs:
    prog = os.path.join(directory, progname)
    path = os.path.join(os.path.abspath("."), os.path.dirname(prog))

    times = []
    print "Corriendo %s" % prog
    for enviro in enviros:
        env = os.environ.copy()
        for key,val in enviro + always_set:
            env[key] = val
        
        data = subprocess.Popen(["./" + progname], env=env, cwd=path, stdout=subprocess.PIPE).communicate()[0]
        measures = []

        for line in data.splitlines():
            groups = re.search("^iteracion total: (.*)$", line)
            if groups: 
                measures.append(groups.group(0))
            groups = re.search("^-->", line)
            if groups:
                print line

        print "{0} ({1}) => {2}".format(directory, ','.join(["{0} = {1}".format(k,v) for k,v in enviro]), measures[-2])
        times.append(time2nanos(measures[-2]))

    base = times[0]
    print "speedups: %s" % " - ".join([str(base / t) for t in times])
