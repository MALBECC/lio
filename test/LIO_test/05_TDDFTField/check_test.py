#!/usr/bin/env python3
import restart
import dipole
import sys

sys.path.insert(0, "../../tests_engine")

dipole.Check("td")
restart.Check()
