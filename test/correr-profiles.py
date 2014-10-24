#!/usr/bin/env python

import os
import subprocess
import re

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

        print measures
        print "{0} ({1}) => {2}".format(directory, ','.join(["{0} = {1}".format(k,v) for k,v in enviro]), measures[-2])
