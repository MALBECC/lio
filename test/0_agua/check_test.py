#!/usr/bin/env python2.7
import sys

sys.path.insert(0,"../tests_engine")
import energy 
import fukui
import mulliken
import forces
import dipole

energy.Check()
fukui.Check()
mulliken.Check()
forces.Check()
dipole.Check()

