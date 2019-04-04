#!/usr/bin/env python2.7

import re
import os
import argparse
import subprocess

def lio_env():
    """"
    Sets lio enviroment variables, adding g2g and
    lioamber to LD_LIBRARY_PATH.
    """

    lioenv = os.environ.copy()
    lioenv["LIOBIN"] = os.path.abspath("../liosolo/liosolo")
    prev = lioenv["LD_LIBRARY_PATH"]
    dirs = ["../g2g", "../lioamber"]
    lioenv["LD_LIBRARY_PATH"] = ":".join([prev] + [os.path.abspath(p) for p in dirs])
    lioenv["LIOHOME"] = os.path.abspath("../")
    return lioenv


def run_lio(dirs_with_tests):
   "Runs all Tests"
   lioenv = lio_env()

   for dir in dirs_with_tests:
      is_file = os.path.isfile(os.path.abspath(dir) + "/run.sh")
      print "Running LIO in",dir

      if is_file:
        execpath = ["./run.sh"]
        process = subprocess.Popen(execpath, env=lioenv, cwd=os.path.abspath(dir))
        process.wait()
        if process.returncode != 0:
           print "Error in this folder."
           continue
        else:
           execpath = ["./check_test.py"]
           check = subprocess.Popen(execpath, env=lioenv, cwd=os.path.abspath(dir))
           check.wait()
      else:
	print("Nothing to do.")


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument("--filter_rx", help="RegExp used to filter which tests are run.", default=".*")
   args = parser.parse_args()
   filterrx = args.filter_rx

   # This obtain tests folder
   subdirs = list(os.walk('LIO_test/'))[0][1]
   dirs_with_tests = sorted([d for d in subdirs if re.search(filterrx,d)])
   total = len(dirs_with_tests)
   for i in range(0,total):
      dirs_with_tests[i] = "LIO_test/" + dirs_with_tests[i]

   # Run lio
   filed = run_lio(dirs_with_tests)

