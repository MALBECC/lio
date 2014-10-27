#!/usr/bin/env python2.7

import os
import subprocess
import re
import click

import itertools as it

MSECS_IN_SEC = 1000000.0
def time2nanos(spec):
    groups = re.search("(?:(\d+)s\. )?(\d+)us\.", spec)
    secs = groups.group(1) or 0
    msecs = groups.group(2)
    return MSECS_IN_SEC * float(secs) + float(msecs)

progname = "correr-profile.sh"
flags_and_options = {
    "OMP_NUM_THREADS":      [1,12],
    "LIO_INNER_THREADS":    [1,12],
    "LIO_OUTER_THREADS":    [1,12],
}

def subdirs_with_benchmark(basedir):
    return [subdir for subdir,dirs,files in os.walk(basedir) if progname in set(files)]

def tuples2str(tuples):
    return ','.join(["{0} = {1}".format(k,v) for k,v in tuples])

def print_file(filename):
    with open(filename) as f:
        for line in f.readlines():
            print "--> " + line,

def get_enviroments():
    return zip(*map(lambda (k,v): zip(it.repeat(k,len(v)),v), flags_and_options.items()))

@click.command()
@click.option("--regex", default=".*", help="Filtro para los tests")
@click.option("--gpu_opts", default="gpu_options", help="Archivo gpu_options a usar (local a la carpeta)")
def benchmark(regex, gpu_opts):
    """ 
    Correr los correr-profile.sh de todas las carpetas que lo posean, 
    y generar un reporte de cada uno.
    """

    testrx = re.compile(regex)
    testdirs = filter(testrx.search, subdirs_with_benchmark("."))

    print "Corriendo %d tests" % len(testdirs)

    for directory in testdirs:
        prog = os.path.join(directory, progname)
        path = os.path.join(os.path.abspath("."), os.path.dirname(prog))

        times = []
        print "Corriendo %s..." % prog
        for enviro in get_enviroments():
            env = os.environ.copy()
            for key,val in enviro:
                env[key] = str(val)
            env["LIO_OPTIONS_FILE"] = gpu_opts            

            print_file(os.path.join(directory,gpu_opts))

            data = subprocess.Popen([os.path.join(".",progname)], env=env, cwd=path, stdout=subprocess.PIPE).communicate()[0]
            measures = []

            for line in data.splitlines():
                groups = re.search("^iteracion total: (.*)$", line)

                if groups: 
                    measures.append(groups.group(0))

                if re.search("^-->", line):
                    print line

            if len(measures) < 2:
                print "No hay resultados para %s, revise el test" % prog
                break

            print "{0} ({1}) => {2}".format(directory, tuples2str(enviro), measures[-2])
            times.append(time2nanos(measures[-2]))
        
        if len(times) == 0:
            print "Ha habido un error en los resultados"
            break

        base = times[0]
        print "speedups: %s" % " - ".join([str(base / t) for t in times])

if __name__ == "__main__":
    benchmark()
