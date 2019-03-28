#!/usr/bin/env python2.7
import sys

sys.path.insert(0,"../tests_engine")
import energy 
import mulliken
import forces
import restart

energy.Check()
restart.Check()
forces.Check()
mulliken.Check()
