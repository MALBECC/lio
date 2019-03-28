#!/usr/bin/env python2.7
import sys

sys.path.insert(0,"../tests_engine")
import energy 
import forces

print "dentro de check tests"

energy.Check()
forces.Check()

