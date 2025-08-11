#!/usr/bin/env python3
"""
Simple PyLIO - Simplified Python wrapper for LIO
Direct access to basic LIO functionality using ctypes
"""

import ctypes
import numpy as np
import os
from pathlib import Path
from typing import Optional, Dict, List, Tuple

class SimplePyLIO:
    """Simplified Python interface to LIO libraries"""
    
    def __init__(self):
        """Inicializar el wrapper de LIO"""
        print("Inicializando SimplePyLIO...")
        self.lib_path = self._find_libraries()
        self.libg2g = None
        self.liblio = None
        self._load_libraries()
    
    def __enter__(self):
        """Context manager entry"""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        try:
            if self.libg2g and hasattr(self.libg2g, 'g2g_deinit_'):
                self.libg2g.g2g_deinit_()
        except:
            pass
    
    def _find_libraries(self):
        """Find LIO installation directory"""
        # Find lio directory from the current path or environment
        current_dir = Path.cwd()
        
        # Try to find lio in path hierarchy
        for parent in [current_dir] + list(current_dir.parents):
            lio_dir = parent / "lio"
            if lio_dir.exists():
                print(f"Found LIO directory at: {lio_dir}")
                self.lio_path = lio_dir
                return lio_dir
            if parent.name == "lio":
                print(f"Found LIO directory at: {parent}")
                self.lio_path = parent
                return parent
        
        # If not found, assume we're in a subdirectory of lio
        lio_dir = Path("/home/nano/lio")
        if lio_dir.exists():
            print(f"Using default LIO directory at: {lio_dir}")
            self.lio_path = lio_dir
            return lio_dir
        
        raise RuntimeError("Cannot find LIO installation directory")
    
    def _load_libraries(self):
        """Load shared libraries"""
        try:
            # Try different possible locations for libraries
            possible_paths = [
                # New build system location
                (self.lio_path / "build" / "lib" / "libg2g.so", 
                 self.lio_path / "build" / "lib" / "liblio-g2g.so"),
                # Original locations
                (self.lio_path / "g2g" / "libg2g.so", 
                 self.lio_path / "lioamber" / "liblio-g2g.so")
            ]
            
            # Load g2g library
            for g2g_path, lio_path in possible_paths:
                if g2g_path.exists():
                    self.libg2g = ctypes.CDLL(str(g2g_path))
                    print(f"✓ Loaded libg2g from {g2g_path}")
                    break
            else:
                print(f"⚠ libg2g.so not found in any location")
            
            # Load lio library  
            for g2g_path, lio_path in possible_paths:
                if lio_path.exists():
                    self.liblio = ctypes.CDLL(str(lio_path))
                    print(f"✓ Loaded liblio from {lio_path}")
                    break
            else:
                print(f"⚠ liblio-g2g.so not found in any location")
                
        except Exception as e:
            print(f"❌ Error loading libraries: {e}")
    
    def available_functions(self):
        """List available functions in loaded libraries"""
        functions = {}
        
        if self.libg2g:
            # Try to access some known g2g functions
            g2g_funcs = []
            for func_name in [
                'g2g_init_', 'g2g_deinit_', 'g2g_parameter_init_',
                'g2g_solve_groups_', 'g2g_reload_atom_positions_',
                'g2g_new_grid_', 'g2g_extern_functional_'
            ]:
                try:
                    func = getattr(self.libg2g, func_name)
                    g2g_funcs.append(func_name)
                except AttributeError:
                    pass
            functions['g2g'] = g2g_funcs
        
        return functions
    
    def call_g2g_init(self):
        """Initialize G2G library"""
        if not self.libg2g:
            raise RuntimeError("libg2g not loaded")
        
        try:
            # Setup function signature
            self.libg2g.g2g_init_.argtypes = []
            self.libg2g.g2g_init_.restype = None
            
            # Call function
            self.libg2g.g2g_init_()
            print("✓ g2g_init_() called successfully")
            return True
            
        except Exception as e:
            print(f"❌ Error calling g2g_init_: {e}")
            return False
    
    def call_g2g_deinit(self):
        """Deinitialize G2G library"""
        if not self.libg2g:
            raise RuntimeError("libg2g not loaded")
        
        try:
            self.libg2g.g2g_deinit_.argtypes = []
            self.libg2g.g2g_deinit_.restype = None
            
            self.libg2g.g2g_deinit_()
            print("✓ g2g_deinit_() called successfully")
            return True
            
        except Exception as e:
            print(f"❌ Error calling g2g_deinit_: {e}")
            return False
    
    def test_basic_workflow(self):
        """Test basic G2G workflow"""
        print("Testing basic G2G workflow...")
        
        success = True
        
        # Test initialization
        if not self.call_g2g_init():
            success = False
        
        # Test deinitialization 
        if not self.call_g2g_deinit():
            success = False
        
        if success:
            print("✓ Basic workflow test passed")
        else:
            print("❌ Basic workflow test failed")
        
        return success


class LIOInputGenerator:
    """Generate LIO input files from Python"""
    
    @staticmethod
    def create_input_file(filename: str, 
                         atoms: List[str],
                         coordinates: np.ndarray,
                         basis: str = "sto-3g",
                         functional: str = "pbe",
                         charge: int = 0,
                         multiplicity: int = 1,
                         **kwargs):
        """
        Create LIO input file
        
        Args:
            filename: Output filename
            atoms: List of atomic symbols
            coordinates: Coordinates in Angstrom (N,3)
            basis: Basis set name
            functional: Exchange-correlation functional
            charge: Molecular charge
            multiplicity: Spin multiplicity
        """
        
        # Map functionals to LIO IDs
        functional_map = {
            'lda': 1,
            'pbe': 2, 
            'b3lyp': 3,
            'pbe0': 4
        }
        
        iexch = functional_map.get(functional.lower(), 2)
        
        # Create namelist
        namelist = f"""&lio
 natom={len(atoms)}
 charge={charge}
 multiplicity={multiplicity}
 basis_set='{basis}'
 iexch={iexch}
 nopt=0
 verbose=1
&end
"""
        
        # Create coordinates section
        coords_section = "\n"
        for i, (atom, coord) in enumerate(zip(atoms, coordinates)):
            coords_section += f"{atom:2s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n"
        
        # Write file
        with open(filename, 'w') as f:
            f.write(namelist)
            f.write(coords_section)
        
        print(f"✓ Created LIO input file: {filename}")
        return filename
    
    @staticmethod
    def create_basis_file(filename: str, basis: str = "sto-3g"):
        """Create basis set file"""
        # This would contain actual basis set data
        # For now, just create a placeholder
        basis_content = f"! Basis set: {basis}\n! (Placeholder - would contain real basis data)\n"
        
        with open(filename, 'w') as f:
            f.write(basis_content)
        
        print(f"✓ Created basis file: {filename}")
        return filename


class LIORunner:
    """Run LIO calculations using subprocess"""
    
    def __init__(self, lio_path: Optional[str] = None):
        if lio_path is None:
            lio_path = Path(__file__).parent.parent
        self.lio_path = Path(lio_path)
        self.liosolo_exe = self.lio_path / "liosolo" / "liosolo"
    
    def run_calculation(self, 
                       input_file: str,
                       output_file: Optional[str] = None,
                       basis_file: Optional[str] = None) -> Dict:
        """
        Run LIO calculation
        
        Args:
            input_file: LIO input file
            output_file: Output file (default: input + .out)
            basis_file: Basis set file
            
        Returns:
            Dictionary with results
        """
        
        if not self.liosolo_exe.exists():
            raise FileNotFoundError(f"liosolo executable not found: {self.liosolo_exe}")
        
        if output_file is None:
            output_file = str(Path(input_file).with_suffix('.out'))
        
        # Build command
        cmd = [str(self.liosolo_exe), "-i", input_file]
        if basis_file:
            cmd.extend(["-b", basis_file])
        
        print(f"Running: {' '.join(cmd)}")
        
        try:
            import subprocess
            
            # Run calculation
            with open(output_file, 'w') as outf:
                result = subprocess.run(
                    cmd,
                    stdout=outf,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=300  # 5 minute timeout
                )
            
            # Parse output
            if result.returncode == 0:
                print(f"✓ Calculation completed successfully")
                results = self._parse_output(output_file)
                results['success'] = True
            else:
                print(f"❌ Calculation failed with return code {result.returncode}")
                print(f"Error: {result.stderr}")
                results = {'success': False, 'error': result.stderr}
            
            return results
            
        except subprocess.TimeoutExpired:
            print("❌ Calculation timed out")
            return {'success': False, 'error': 'Timeout'}
        except Exception as e:
            print(f"❌ Error running calculation: {e}")
            return {'success': False, 'error': str(e)}
    
    def _parse_output(self, output_file: str) -> Dict:
        """Parse LIO output file"""
        results = {}
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Look for energy
            for line in content.split('\n'):
                if 'Total Energy' in line or 'SCF Energy' in line:
                    try:
                        energy = float(line.split()[-1])
                        results['energy'] = energy
                    except:
                        pass
                
                if 'Dipole' in line:
                    # Parse dipole moment
                    pass
            
            results['output'] = content
            
        except Exception as e:
            results['parse_error'] = str(e)
        
        return results


def example_usage():
    """Example of using the simplified wrapper"""
    
    print("SimplePyLIO Example")
    print("===================")
    
    # 1. Test library loading
    lio = SimplePyLIO()
    
    # 2. Check available functions
    funcs = lio.available_functions()
    print(f"Available functions: {funcs}")
    
    # 3. Test basic workflow
    lio.test_basic_workflow()
    
    # 4. Create input files
    atoms = ['O', 'H', 'H']
    coordinates = np.array([
        [0.0000,  0.0000,  0.1173],
        [0.0000,  0.7572, -0.4692],
        [0.0000, -0.7572, -0.4692]
    ])
    
    # Generate input
    input_gen = LIOInputGenerator()
    input_file = input_gen.create_input_file(
        "water.inp",
        atoms=atoms,
        coordinates=coordinates,
        basis="sto-3g",
        functional="pbe"
    )
    
    # Create basis file
    basis_file = input_gen.create_basis_file("water.bas", "sto-3g")
    
    # 5. Run calculation (if executable exists)
    runner = LIORunner()
    if runner.liosolo_exe.exists():
        print("\nRunning LIO calculation...")
        results = runner.run_calculation(input_file, basis_file=basis_file)
        
        if results.get('success'):
            print(f"✓ Energy: {results.get('energy', 'N/A')} hartree")
        else:
            print(f"❌ Calculation failed: {results.get('error', 'Unknown error')}")
    else:
        print(f"⚠ liosolo executable not found at {runner.liosolo_exe}")
        print("  Compile LIO first with 'make liosolo'")


if __name__ == "__main__":
    example_usage()
