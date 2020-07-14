#!/usr/bin/env python3
import sys

sys.path.insert(0,"../../tests_engine")
import energy 
import forces
import mulliken
import dipole

energy.Check()
forces.Check()
mulliken.Check()
dipole.Check()
