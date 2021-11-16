#!/usr/bin/env python3
import energy
import dipole
import forces
import mulliken
import fukui
import sys

sys.path.insert(0, "../../tests_engine")

energy.Check()
fukui.Check()
mulliken.Check()
forces.Check()
dipole.Check()
