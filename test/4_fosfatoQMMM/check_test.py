#!/usr/bin/env python2.7
import sys

sys.path.insert(0,"../tests_engine")
import energy 
import forces
import mulliken

energy.Check()
forces.Check()
mulliken.Check()

