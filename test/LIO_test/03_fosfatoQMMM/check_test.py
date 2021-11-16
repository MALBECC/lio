#!/usr/bin/env python3
import energy
import dipole
import mulliken
import forces
import sys

sys.path.insert(0, "../../tests_engine")

energy.Check()
forces.Check()
mulliken.Check()
dipole.Check()
