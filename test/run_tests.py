#!/usr/bin/env python2.7

import argparse
import fileinput
import os
import re
import subprocess
import sys

from collections import namedtuple

Summary = namedtuple('Summary', [
    'iterations','converged','total_time','avg_time','scf_energy'
])

SEC_TO_USEC = 1000*1000

def avg(l):
    return (1. * sum(l)) / len(l)

def get_statistics(out_file):
    "Get statistics for the LIO run out file"

    iterations = []
    scf_energy = []
    iteration_time = []
    convergence_at = []

    for line in out_file.readlines():
        # Iteration output line
        m = re.match("\s+iter\s+(\d+)", line)
        if m:
            iterations.append(float(m.group(1)))

        # Correlation Energy output line
        m = re.match(r"\s+SCF ENRGY=\s+([0-9.-]+)", line)
        if m:
            scf_energy.append(float(m.group(1)))

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

    if len(scf_energy) < 1:
        return None

    return Summary(
            iterations=max(iterations),\
            converged=len(convergence_at) > 0,\
            total_time=sum(iteration_time),\
            avg_time=avg(iteration_time),\
            scf_energy=scf_energy[-1])

# Tolerance parameters
OVERTIME = 1000*1000

def print_test_summary(run_summary, ok_summary):
    print "\n\tResult = %r" % (run_summary,)
    print "\tExpected = %r\n" % (ok_summary,)
    print "\tSCF Diff: %f" % abs(run_summary.scf_energy - ok_summary.scf_energy)

    per = (run_summary.total_time - ok_summary.total_time) / ok_summary.total_time
    print "\tTime increase: %f %%" % (100.0 * per)

def acceptable(run_summary, ok_summary):
    "Returns whether the result is within the test bounds"

    if abs(run_summary.scf_energy - ok_summary.scf_energy)*627 > 0.2:
        return "invalid numerical result"

    if not run_summary.converged:
        return "no convergence"

    if run_summary.iterations > 2 * ok_summary.iterations:
        return "took too long to converge"

    return None

def lio_env():
    """"
    Set lio enviroment variables, including adding g2g and
    lioamber to LD_LIBRARY_PATH.
    """
    lioenv = os.environ.copy()
    lioenv["LIOBIN"] = os.path.abspath("../liosolo/liosolo")
    prev = lioenv["LD_LIBRARY_PATH"]
    dirs = ["../g2g", "../lioamber"]
    lioenv["LD_LIBRARY_PATH"] = ":".join([prev] + [os.path.abspath(p) for p in dirs])
    return lioenv

def lio_run(dir, lioenv):
    "Run Lio script"
    execpath = ["./correr.sh"]
    process = subprocess.Popen(execpath, env=lioenv, cwd=os.path.abspath(dir))
    process.wait()
    return process.returncode

def recompile_lio(dir, lioenv):
    "Rebuild lio with the enviroment variables given"
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
    failed = 0

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

        if not os.path.isfile(os.path.join(dir, "salida.ok")):
            print "\tFailed because no ideal output in directory"
            return

        ok_summary = None
        with open(os.path.join(dir, "salida.ok"),"r") as f:
            ok_summary = get_statistics(f)
            if not ok_summary:
                print "\tFailed because ideal output doesn't have XC energy info ('XC energy: ' lines)"
                return

        out_summary = None
        with open(os.path.join(dir, "salida"),"r") as f:
            out_summary = get_statistics(f)
            if not out_summary:
                print "\tFailed because output of run doesn't have XC energy info ('XC energy: ' lines)"
                return

        print_test_summary(out_summary, ok_summary)
        veredict = acceptable(out_summary, ok_summary)
        if veredict:
            print "\tFailed because not acceptable result: %s" % veredict
            failed += 1
        else:
            print "\tPassed\n"

    return failed

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filter_rx", help="Expresion regular para filtrar que tests se corren", default=".*")
    args = parser.parse_args()
    filterrx = args.filter_rx

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    subdirs = list(os.walk('.'))[0][1]
    dirs_with_tests = sorted([d for d in subdirs if re.search(filterrx,d)])

    failed = run_tests(dirs_with_tests)
    if failed > 0:
        print "%d tests fallaron..." % failed
        sys.exit(1)
