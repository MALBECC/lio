#!/usr/bin/env python

import os
import subprocess
import itertools

def set_options(flag):
   opt = ["%s=%s" % (k,v) for (k,v) in flag_set.items()]
   opt_str = "%s" % (" ".join(opt),)
   options = opt_str.rstrip()
   return options

def compilar(options):
   liodir = os.path.abspath("../../")
   devnull = open(os.devnull, "w")
   cmd = ["./build.sh", liodir, options]
   process = subprocess.Popen(cmd, stdout=devnull, stderr=devnull)
   process.wait()
   return process.returncode

def run_lio():
   cmd = ["./new_tests.py"]
   process = subprocess.Popen(cmd, cwd=os.path.abspath("."))
   process.wait()
   return process.returncode

if __name__ == "__main__":
   comp = ["cuda","intel","precision"]
   seq = list(itertools.product(["0","1"],repeat=len(comp)))
   all_sets = []

   for cases in seq:
      compile_opts = dict([(comp[i],cases[i]) for i in xrange(0,len(comp))])
      if compile_opts["cuda"] == "1":
         compile_opts["cuda"] = "2"

      all_sets.append(compile_opts)

   for flag_set in all_sets:
      opts = set_options(flag_set)
      print "Compiling LIO with Options: %s" % opts.rstrip()
      error = compilar(opts)
      if not error:
         print "\tSuccessfully compiled"
      else:
         print "\tError!"

      error = run_lio()
      if not error:
         print "run_lio successfully finished"
      else:
         print "run_lio Error!"
         exit(-1)
