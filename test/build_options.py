#!/usr/bin/env python
import itertools
import random
import sys
if __name__ == '__main__':
    opts = ["cpu","cuda","magma","cublas","libxc"]

    seq = list(itertools.product(["0","1"],repeat=len(opts)))
    all_sets = []
    for cases in seq:
        compile_opts = dict([(opts[i],cases[i]) for i in xrange(0,len(opts))]) 
        if compile_opts["cpu"] == compile_opts["cuda"]:
            continue
        all_sets.append(compile_opts)
 
    if len(sys.argv) > 1:
        all_sets = random.sample(all_sets,int(sys.argv[1]))
    
    for flag_set in all_sets:
        opts_str = ["%s=%s" % (k,v) for (k,v) in flag_set.items()]
        print "%s" % (" ".join(opts_str),)
