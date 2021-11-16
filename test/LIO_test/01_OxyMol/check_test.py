#!/usr/bin/env python3
import energy
import forces
import mulliken
import sys

sys.path.insert(0, "../../tests_engine")

energy.Check()
mulliken.Check()
forces.Check()
