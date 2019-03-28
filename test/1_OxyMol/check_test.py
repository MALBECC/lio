#!/usr/bin/env python2.7
import sys

sys.path.insert(0,"../tests_engine")
import energy 
import mulliken
import forces

energy.Check()
mulliken.Check()
forces.Check()

