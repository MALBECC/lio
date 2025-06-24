#!/usr/bin/env python3
import sys

sys.path.insert(0,"../../tests_engine")
TB_thresholds = {
    'default': 0.05,  # Loosen the default for all non-specified energies
    'total': 0.05    # Also loosen the total energy threshold
}
import energy 

energy.Check(custom_thresholds=TB_thresholds)

