#!/usr/bin/env python

import fileinput
import os
import re
import subprocess
import sys

from collections import namedtuple

Summary = namedtuple('Summary', [
    'iterations','converged','total_time','avg_time','xc_energy'
])

SEC_TO_USEC = 1000*1000

def avg(l):
    return (1. * sum(l)) / len(l)

def get_statistics(out_file):
    "Get statistics for the LIO run out file"

    iterations = []
    xc_energy = []
    iteration_time = []
    convergence_at = []

    for line in out_file.readlines():
        # Iteration output line
        m = re.match("\s+iter\s+(\d+)", line)
        if m:
            iterations.append(float(m.group(1)))

        # Correlation Energy output line
        m = re.match("XC energy: ([0-9.-]+)", line)
        if m:
            xc_energy.append(float(m.group(1)))

        # Iteration time output line
        m = re.match("TIMER \[Total iter\]: (?:(\d+)s. )?(\d+)us.",line)
        if m:
            sec,usec = m.groups()
            sec = sec or 0
            iteration_time.append(float(sec)*SEC_TO_USEC+float(usec))

        # Iteration convergence output line
        m = re.match(" CONVERGED AT \s+(\d+) ITERATIONS", line)
        if m:
            convergence_at.append(float(m.group(1)))

    return Summary(
            iterations=max(iterations),\
            converged=len(convergence_at) > 0,\
            total_time=sum(iteration_time),\
            avg_time=avg(iteration_time),\
            xc_energy=xc_energy[-1])

# Tolerance parameters
EPS = 1e-5
OVERTIME = 1000*1000

def print_test_summary(run_summary, ok_summary):
    print "\n\tResult = %r" % (run_summary,)
    print "\tExpected = %r\n" % (ok_summary,)
    print "\tXC Diff: %f" % abs(run_summary.xc_energy - ok_summary.xc_energy)

    per = abs(run_summary.total_time - ok_summary.total_time) / ok_summary.total_time
    print "\tTime increase: %f %%" % (100.0 * per)

def acceptable(run_summary, ok_summary):
    "Returns whether the result is within the test bounds"
    if run_summary.avg_time > ok_summary.avg_time + OVERTIME:
        return False

    if abs(run_summary.xc_energy - ok_summary.xc_energy) > EPS:
        return False

    return True

def lio_env():
    lioenv = os.environ.copy()
    lioenv["LIOBIN"] = os.path.abspath("../liosolo/liosolo")
    return lioenv

def lio_run(dir, lioenv):
    execpath = ["./correr.sh"]
    process = subprocess.Popen(execpath, env=lioenv, cwd=os.path.abspath(dir))
    process.wait()
    return process.returncode

def recompile_lio(dir, lioenv):
    g2gpath = os.path.abspath("../g2g")
    try:
        devnull = open(os.devnull, "w")
        opts = open(os.path.join(dir,"compile_options"),"r").read()
        print "\tOptions: %s" % opts.rstrip()
        cmd = ["./compilar.sh", g2gpath, opts]
        process = subprocess.Popen(cmd, stdout=devnull, stderr=devnull)
        process.wait()
        return process.returncode
    except IOError:
        print "\tNo compile-options!"
        return -1

def run_tests(dirs_with_tests):
    "Run all tests in the folder"

    res = []
    lioenv = lio_env()
    for dir in dirs_with_tests:
        print("Running %s..." % dir)

        print "\tRecompiling..."
        errcode = recompile_lio(dir, lioenv)
        if errcode != 0:
            print "\tFailed to recompile"
            break
        print "\tRecompiling Done\n\tRunning..."

        errcode = lio_run(dir, lioenv)
        if errcode != 0:
            print "\tFailed to run with errcode %d" % errcode
            break

        ok_summary = None
        with open(os.path.join(dir, "salida.ok"),"r") as f:
            ok_summary = get_statistics(f)

        out_summary = None
        with open(os.path.join(dir, "salida"),"r") as f:
            out_summary = get_statistics(f)

        print_test_summary(out_summary, ok_summary)
        if not acceptable(out_summary, ok_summary):
            print "\tFailed because not acceptable result"

        print "\tPassed\n"

if __name__ == '__main__':
    filterrx = ".*"
    if len(sys.argv) > 1:
        filterrx = sys.argv[1]

    subdirs = list(os.walk('.'))[0][1]
    dirs_with_tests = sorted([d for d in subdirs if re.match(filterrx,d)])

    run_tests(dirs_with_tests)
