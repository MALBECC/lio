#!/usr/bin/env python2.7

import os
import subprocess
import re
import multiprocessing
import optparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import itertools as it

from datetime import datetime

MSECS_IN_SEC = 1000000.0
def time2nanos(spec):
    groups = re.search("(?:(\d+)s\. )?(\d+)us\.", spec)
    secs = groups.group(1) or 0
    msecs = groups.group(2)
    return MSECS_IN_SEC * float(secs) + float(msecs)

progname = "correr-profile.sh"

def subdirs_with_benchmark(basedir):
    return [subdir for subdir,dirs,files in os.walk(basedir) if progname in set(files)]

def tuples2str(tuples):
    return ','.join(["{0} = {1}".format(k,v) for k,v in tuples])

def print_file(filename):
    with open(filename) as f:
        for line in f.readlines():
            print "--> " + line,

def process_lio_output(output):
    measures, info = [],[]
    for line in output.splitlines():
        groups = re.search("^iteracion total: (.*)$", line)

        if groups: 
            measures.append(groups.group(0))

        if re.search("^-->", line):
            info.append(line)
    return measures, '\n'.join(info)

def timestamp():
    return '-'.join(str(datetime.today()).split(" "))

def replace_special(s):
    return re.sub("[.:-]","-", s)

def plot_scalability(speedups, expname):
    plt.clf()
    plt.xlabel("Cantidad de threads")
    plt.ylabel("Speedup")
    plt.title("Prueba de escalabilidad en cores para %s" % expname)
    threads = xrange(1,len(speedups)+1)
    plt.plot(threads, speedups, "ro-", label="Experimental")
    plt.plot(threads, threads, "bo-", label="Teorico")
    plt.legend(loc=2)

    name = replace_special("escalabilidad-%s-%s.png" % (expname,timestamp()))
    plt.savefig(os.path.join("escalabilidad",name))

def get_enviroments(dic, keylist=None):
    if keylist is None:
        keylist = dic.keys()
    
    if keylist == []: 
        return [[]]

    key, rest, ret = keylist[0], keylist[1:], []
    sublists = get_enviroments(dic, rest)
    for value in dic[key]:
        for sublist in sublists:
            ret.append([(key,value)] + sublist)
    return ret 

def average(l):
    return sum(l) / len(l)

def benchmark(regex, gpu_opts, threadlist, thresholdlist, offsetlist, plot_scale):
    """ 
    Correr los correr-profile.sh de todas las carpetas que lo posean, 
    y generar un reporte de cada uno.
    """
    
    testrx = re.compile(regex)
    testdirs = filter(testrx.search, subdirs_with_benchmark("."))

    env_flags = {
        "OMP_NUM_THREADS": threadlist,
        "KMP_AFFINITY": ["granularity=fine,scatter"],
        "LIO_OPTIONS_FILE": [gpu_opts],
        "LIO_SPLIT_THRESHOLD": thresholdlist,
        "LIO_MINCOST_OFFSET": offsetlist,
    }

    for directory in testdirs:
        prog = os.path.join(directory, progname)
        path = os.path.join(os.path.abspath("."), os.path.dirname(prog))

        times = []
        print "Corriendo %s..." % prog
        for enviro in get_enviroments(env_flags):
            env = os.environ.copy()
            for key,val in enviro:
                env[key] = str(val)

            print_file(os.path.join(directory,gpu_opts))

            data, _ = subprocess.Popen([os.path.join(".",progname)],
                    env=env, cwd=path, stdout=subprocess.PIPE).communicate()
            measures,info = process_lio_output(data)

            print info

            if len(measures) < 2:
                print "No hay resultados para %s, revise el test" % prog
                break

            measure = average(measures[1:-2])
            print "{0} ({1}) => {2}".format(directory, tuples2str(enviro), measure)
            times.append(time2nanos(measure))
        
        if len(times) == 0:
            print "Ha habido un error en los resultados"
            break

        speedups = [max(times) / t for t in times]
        print "speedups: %s" % " - ".join(map(str, speedups))
        if plot_scale:
            plot_scalability(speedups, directory[2:])

def parse_range(rangestr):
    parts = rangestr.split(":")
    if len(parts) < 3:
        parts.append("1")
    start,end,skip = map(int,parts)
    return range(start,end+1,skip)

def parse_range_string(rangestr):
    res = []
    if rangestr[0] == "[":
        parts = rangestr[1:-1].split(",")
        for part in parts:
            if part.isdigit():
                res.append(int(part))
            else:
                res += parse_range(part)
    elif rangestr.isdigit():
        res = [int(rangestr)]
    else:
        res += parse_range(rangestr)

    return sorted(res)

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-r", "--regex", dest="regex", default=".*", help="Filtrar los test con expresion regular")
    parser.add_option("-g", "--gpu_opts", dest="gpu_opts", help="Archivo gpu_options a usar (local a la carpeta)")
    parser.add_option("-t", "--threads", dest="threads",default="[1,%s]" % multiprocessing.cpu_count(), help="Lista de threads a usar, como rango (N:M:J)")
    parser.add_option("-a", "--plot_threadscale", action="store_true", default=False, help="Graficar escalabilidad intermedia")
    parser.add_option("-s", "--thresholds", dest="thresholds", default="100:100", help="Threshold para considerar un grupo como chico, como rango (N:M:J)")
    parser.add_option("-o", "--offsets", dest="offsets", default="50000:50000", help="Compensacion de costo para cubos chicos, como rango (N:M:J)")
    (options, args) = parser.parse_args()

    threadlist = parse_range_string(options.threads)
    thresholdlist = parse_range_string(options.thresholds)
    offsetlist = parse_range_string(options.offsets)

    benchmark(options.regex,options.gpu_opts, threadlist, thresholdlist, offsetlist, options.plot_threadscale)
