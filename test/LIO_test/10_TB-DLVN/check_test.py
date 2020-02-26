#!/usr/bin/env python2.7
import sys

sys.path.insert(0,"../../tests_engine")
import energy 
import dipole

energy.Check()
dipole.Check("td")
