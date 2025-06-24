#!/usr/bin/env python3
import sys

sys.path.insert(0,"../../tests_engine")
import dipole
import restart

dipole.Check("td")
restart.Check()

