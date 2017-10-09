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
    return "\n" + '\n'.join(["\t{0} = {1}".format(k,v) for k,v in tuples]) + "\n"

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

def replace_special(s):
    return re.sub("[.:-]","-", s)

def timestamp():
    return replace_special('-'.join(str(datetime.today()).split(" ")))

def plot_scalability(speedups, expname):
    plt.clf()
    plt.xlabel("Cantidad de threads")
    plt.ylabel("Speedup")
    plt.title("Prueba de escalabilidad en cores para %s" % expname)
    threads = xrange(1,len(speedups)+1)
    plt.plot(threads, speedups, "ro-", label="Experimental")
    plt.plot(threads, threads, "bo-", label="Teorico")
    plt.legend(loc=2)

    name = replace_special("escalabilidad-%s-%s" % (expname,timestamp()))
    plt.savefig(os.path.join("escalabilidad","%s.png" % name))

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

def computeresult(l):
    return sum(l) / len(l)

def save_lio_output(name, output):
    filename = os.path.join("outputs","experiment-%s-%s.txt" % (name, timestamp()))
    with open(filename, "w") as f:
        print >> f, output

def benchmark(regex, gpu_opts, threadlist, thresholdlist, offsetlist, plot_scale, g2gpath, save_outputs):
    """ 
    Correr los correr-profile.sh de todas las carpetas que lo posean, 
    y generar un reporte de cada uno.
    """
    
    testrx = re.compile(regex)
    testdirs = filter(testrx.search, subdirs_with_benchmark("."))

    ldpath = ":".join([g2gpath,os.environ.get("LD_LIBRARY_PATH")])

    env_flags = {
        "OMP_NUM_THREADS": threadlist,
        "KMP_AFFINITY": ["granularity=fine,scatter"],
        "LIO_OPTIONS_FILE": [gpu_opts],
        "LIO_SPLIT_THRESHOLD": thresholdlist,
        "LIO_MINCOST_OFFSET": offsetlist,
        "LD_LIBRARY_PATH": [ldpath],
        "OMP_NESTED": ["false"],
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

            if save_outputs:
                save_lio_output(directory[2:], data)

            print info

            if len(measures) < 2:
                print "No hay resultados para %s, revise el test" % prog
                break

            measure = computeresult(map(time2nanos,[measures[-2]]))
            print "{0} ({1}) => {2} mus.".format(directory, tuples2str(enviro), measure)
            times.append(measure)
        
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
    parser.add_option("-s", "--thresholds", dest="thresholds", default="150:150", help="Threshold para considerar un grupo como chico, como rango (N:M:J)")
    parser.add_option("-o", "--offsets", dest="offsets", default="200000:200000", help="Compensacion de costo para cubos chicos, como rango (N:M:J)")
    parser.add_option("-d", "--directory", dest="g2g_path", default=("./lio-g2gs/latest"), help="Direccion donde encontrar los .so de G2G y lioamber")
    parser.add_option("-q", "--save_outputs", action="store_true", default=False, help="Guardar los outputs de corridas")
    (options, args) = parser.parse_args()

    threadlist = parse_range_string(options.threads)
    thresholdlist = parse_range_string(options.thresholds)
    offsetlist = parse_range_string(options.offsets)

    g2g_path = os.path.abspath(options.g2g_path)
    benchmark(options.regex,options.gpu_opts, threadlist, thresholdlist,\
        offsetlist, options.plot_threadscale, g2g_path, options.save_outputs)
