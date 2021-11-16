#!/usr/bin/env python3

import os
import subprocess
import itertools


def set_options(flag):
    opt = ["%s=%s" % (k, v) for (k, v) in flag_set.items()]
    opt_str = "%s" % (" ".join(opt),)
    options = opt_str.rstrip()
    return options


def compile_lio(options):
    liodir = os.path.abspath("../")
    devnull = open(os.devnull, "w")
    cmd = ["tests_engine/build.sh", liodir, options]
    process = subprocess.Popen(cmd, stdout=devnull, stderr=devnull)
    process.wait()
    return process.returncode


def run_lio():
    cmd = ["./new_tests.py"]
    process = subprocess.Popen(cmd, cwd=os.path.abspath("."))
    process.wait()
    return process.returncode


# Verifies if CUDA is installed.


def cuda_is_installed():
    devnull = open(os.devnull, "wb")
    process = subprocess.Popen(
        ["nvcc --version"], shell=True, stdout=devnull, stderr=devnull
    )
    try:
        stdout, stderr = process.communicate()
    except BaseException:
        process.kill()
        process.wait()
        raise
    retcode = process.poll()

    is_installed = False
    if retcode == 0:
        is_installed = True

    return is_installed


# Checks if Intel compilers are present.


def intel_is_installed():
    devnull = open(os.devnull, "wb")
    process = subprocess.Popen(
        ["icc --version"], shell=True, stdout=devnull, stderr=devnull
    )
    try:
        stdout, stderr = process.communicate()
    except BaseException:
        process.kill()
        process.wait()
        raise
    retcode = process.poll()

    is_installed = False
    if retcode == 0:
        is_installed = True

    return is_installed


if __name__ == "__main__":
    comp = ["precision"]

    if cuda_is_installed():
        comp.append("cuda")
    if intel_is_installed():
        comp.append("intel")
    seq = list(itertools.product(["0", "1"], repeat=len(comp)))
    all_sets = []

    for cases in seq:
        compile_opts = dict([(comp[i], cases[i]) for i in range(0, len(comp))])
        if cuda_is_installed():
            if compile_opts["cuda"] == "1":
                compile_opts["cuda"] = "2"

        if intel_is_installed():
            if compile_opts["intel"] == "1":
                compile_opts["intel"] = "2"

        all_sets.append(compile_opts)

    for flag_set in all_sets:
        opts = set_options(flag_set)
        if not cuda_is_installed():
            opts = opts + " cuda=0 "
        print("Compiling LIO with Options: %s" % opts.rstrip())
        error = compile_lio(opts)
        if not error:
            print("\tSuccessfully compiled.")
        else:
            print("\tError!")

        error = run_lio()
        if not error:
            print("run_lio successfully finished.")
        else:
            print("run_lio Error!")
            exit(-1)
