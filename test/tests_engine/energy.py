#!/usr/bin/env python3
import re
from typing import Dict, List, TextIO, Optional # Import Optional

# Configuration is at the top, making it easy to see the defaults.
ENERGY_PATTERNS = [
    ('total',         r"\s+Total energy =\s+([0-9.-]+)"),
    ('one_electron',  r"\s+One electron =\s+([0-9.-]+)"),
    ('coulomb',       r"\s+Coulomb\s+ =\s+([0-9.-]+)"),
    ('nuclear',       r"\s+Nuclear\s+ =\s+([0-9.-]+)"),
    ('exch_corr',     r"\s+Exch. Corr.\s+ =\s+([0-9.-]+)"),
    ('exact_exch',    r"\s+Exact. Exc.\s+ =\s+([0-9.-]+)"),
    ('qmmm_nuclear',  r"\s+QM-MM nuc.\s+ =\s+([0-9.-]+)"),
    ('qmmm_electron', r"\s+QM-MM elec.\s+ =\s+([0-9.-]+)"),
    ('dftd3',         r"\s+DFTD3 Energy =\s+([0-9.-]+)"),
]

# Default thresholds are still defined here.
THRESHOLDS = {
    'total': 1.6e-4,
    'default': 1.5e-2
}

def obtain_energies(file_in: TextIO) -> Dict[str, float]:
    """Parses an output file and extracts energy values."""
    energies = {key: 0.0 for key, _ in ENERGY_PATTERNS}
    compiled_patterns = [(key, re.compile(pattern)) for key, pattern in ENERGY_PATTERNS]

    for line in file_in:
        for key, pattern in compiled_patterns:
            if match := pattern.match(line):
                energies[key] = float(match.group(1))
                break 
    return energies

# CHANGE 1: `compare_energies` now takes a `thresholds` dictionary as an argument.
def compare_energies(
    test_energies: Dict[str, float], 
    ref_energies: Dict[str, float], 
    thresholds: Dict[str, float]
) -> List[str]:
    """
    Compares two dictionaries of energies using a given set of thresholds.

    Returns a list of error messages for any values that are outside the
    defined thresholds.
    """
    discrepancies = []
    
    if test_energies.keys() != ref_energies.keys():
        discrepancies.append("The set of energy types in the two files is different.")
        return discrepancies

    for key, test_value in test_energies.items():
        ref_value = ref_energies.get(key, 0.0)
        # It now uses the 'thresholds' argument, not the global constant.
        threshold = thresholds.get(key, thresholds['default'])
        
        if abs(test_value - ref_value) > threshold:
            msg = (
                f"Error in '{key}':\n"
                f"  ├─ Test value:     {test_value}\n"
                f"  └─ Reference value:  {ref_value}"
            )
            discrepancies.append(msg)
            
    return discrepancies

# CHANGE 2: `Check` now accepts an optional `custom_thresholds` argument.
def Check(
    test_filepath: str = "output", 
    ref_filepath: str = "output.ok", 
    custom_thresholds: Optional[Dict[str, float]] = None
):
    """
    Checks energy values in an output file against a reference file.

    Args:
        test_filepath (str): Path to the output file to be checked.
        ref_filepath (str): Path to the reference '.ok' file.
        custom_thresholds (Optional[Dict]): A dictionary of thresholds to override
                                             the defaults. E.g., {'total': 1e-5}.
    """
    try:
        with open(test_filepath, 'r') as f_test:
            test_energies = obtain_energies(f_test)
        
        with open(ref_filepath, 'r') as f_ref:
            ref_energies = obtain_energies(f_ref)

    except FileNotFoundError as e:
        print(f"File not found: {e.filename}")
        print("Test Energy:     ERROR")
        return

    # CHANGE 3: Logic to merge default and custom thresholds.
    # Start with a copy of the defaults.
    effective_thresholds = THRESHOLDS.copy()
    if custom_thresholds:
        # Update the copy with any user-provided values.
        # .update() will overwrite existing keys and add new ones.
        effective_thresholds.update(custom_thresholds)

    # Pass the final, effective thresholds to the comparison function.
    errors = compare_energies(test_energies, ref_energies, effective_thresholds)

    if not errors:
        print("Test Energy:     OK")
    else:
        print("Test Energy:     ERROR")
        for error_msg in errors:
            print(error_msg)

if __name__ == '__main__':
    print("Running energy check as a standalone script...")
    Check()
