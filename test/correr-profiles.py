#!/usr/bin/env python2.7

import os
import subprocess
import re
import multiprocessing
import optparse

import numpy as np
import matplotlib.pyplot as plt

import itertools as it

MSECS_IN_SEC = 1000000.0
def time2nanos(spec):
    groups = re.search("(?:(\d+)s\. )?(\d+)us\.", spec)
    secs = groups.group(1) or 0
    msecs = groups.group(2)
    return MSECS_IN_SEC * float(secs) + float(msecs)

progname = "correr-profile.sh"
env_flags = [
    "OMP_NUM_THREADS",
    "LIO_INNER_THREADS",
    "LIO_OUTER_THREADS",
]

def subdirs_with_benchmark(basedir):
    return [subdir for subdir,dirs,files in os.walk(basedir) if progname in set(files)]

def tuples2str(tuples):
    return ','.join(["{0} = {1}".format(k,v) for k,v in tuples])

def print_file(filename):
    with open(filename) as f:
        for line in f.readlines():
            print "--> " + line,

def kmp_affinity_value():
    return ("KMP_AFFINITY", "granularity=fine,scatter")

def get_enviroments(threadlist):
    res = []
    for thread in threadlist:
        res.append([(key,thread) for key in env_flags] + [kmp_affinity_value()])
    return res

def process_lio_output(output):
    measures, info = [],[]
    for line in output.splitlines():
        groups = re.search("^iteracion total: (.*)$", line)

        if groups: 
            measures.append(groups.group(0))

        if re.search("^-->", line):
            info.append(line)
    return measures, '\n'.join(info)

def plot_scalability(speedups):
    plt.xlabel("Cantidad de threads")
    plt.ylabel("Speedup")
    threads = xrange(1,len(speedups)+1)
    plt.plot(threads, speedups, label="Experimental")
    plt.plot(threads, threads, label="Teorico")
    plt.legend()
    plt.show()

def benchmark(regex, gpu_opts, threads, threadscale):
    """ 
    Correr los correr-profile.sh de todas las carpetas que lo posean, 
    y generar un reporte de cada uno.
    """
    
    testrx = re.compile(regex)
    testdirs = filter(testrx.search, subdirs_with_benchmark("."))
    threadlist = xrange(1,threads+1) if threadscale else [1,threads]

    for directory in testdirs:
        prog = os.path.join(directory, progname)
        path = os.path.join(os.path.abspath("."), os.path.dirname(prog))

        times = []
        print "Corriendo %s..." % prog
        for enviro in get_enviroments(threadlist):
            env = os.environ.copy()
            for key,val in enviro:
                env[key] = str(val)

            env["LIO_OPTIONS_FILE"] = gpu_opts 

            print_file(os.path.join(directory,gpu_opts))

            data = subprocess.Popen([os.path.join(".",progname)],
                    env=env, cwd=path, stdout=subprocess.PIPE).communicate()[0]
            measures,info = process_lio_output(data)

            print info

            if len(measures) < 3:
                print "No hay resultados para %s, revise el test" % prog
                break

            measure = measures[-3]
            print "{0} ({1}) => {2}".format(directory, tuples2str(enviro), measure)
            times.append(time2nanos(measure))
        
        if len(times) == 0:
            print "Ha habido un error en los resultados"
            break

        speedups = [max(times) / t for t in times]
        print "speedups: %s" % " - ".join(map(str, speedups))
        if threadscale:
            plot_scalability(speedups)

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-r", "--regex", dest="regex", help="Filtrar los test con expresion regular")
    parser.add_option("-g", "--gpu_opts", dest="gpu_opts", help="Archivo gpu_options a usar (local a la carpeta)")
    parser.add_option("-t", "--threads", dest="threads",default=multiprocessing.cpu_count(), help="Maxima cantidad de threads a usar")
    parser.add_option("-a", "--threadscale", action="store_true", default=False, help="Graficar escalabilidad intermedia")
    (options, args) = parser.parse_args()
    benchmark(options.regex,options.gpu_opts,int(options.threads), options.threadscale)
