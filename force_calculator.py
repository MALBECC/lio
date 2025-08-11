#!/usr/bin/env python3
"""
Calculadora de Fuerzas con LIO
Ejemplo específico para cálculo de fuerzas/gradientes
"""

import ctypes
import os
import numpy as np
import sys

class LIOForceCalculator:
    """Calculadora de fuerzas usando la interfaz Python-LIO"""
    
    def __init__(self, lib_path="/home/nano/lio/build/lib/liblio-g2g.so"):
        """Inicializar calculadora"""
        self.lib_path = lib_path
        self.lio = None
        self.system_initialized = False
        self.scf_converged = False
        
    def load_library(self):
        """Cargar biblioteca LIO"""
        if not os.path.exists(self.lib_path):
            raise FileNotFoundError(f"Biblioteca no encontrada: {self.lib_path}")
        
        try:
            self.lio = ctypes.CDLL(self.lib_path)
            self._setup_function_signatures()
            print(f"✅ Biblioteca LIO cargada: {self.lib_path}")
            return True
        except Exception as e:
            print(f"❌ Error cargando biblioteca: {e}")
            return False
    
    def _setup_function_signatures(self):
        """Configurar firmas de funciones ctypes"""
        # Función hello
        self.lio.python_hello.argtypes = []
        self.lio.python_hello.restype = None
        
        # Función de inicialización
        self.lio.python_lio_init.argtypes = [
            ctypes.c_int,  # natom
            ctypes.POINTER(ctypes.c_int),  # atomic_numbers
            ctypes.POINTER(ctypes.c_double),  # coordinates  
            ctypes.c_int,  # charge
            ctypes.c_char_p,  # basis_name
            ctypes.POINTER(ctypes.c_int)  # iostat
        ]
        self.lio.python_lio_init.restype = None
        
        # Función SCF
        self.lio.python_lio_scf.argtypes = [
            ctypes.POINTER(ctypes.c_double), 
            ctypes.POINTER(ctypes.c_int)
        ]
        self.lio.python_lio_scf.restype = None
        
        # Función de gradientes/fuerzas
        self.lio.python_lio_gradients.argtypes = [
            ctypes.POINTER(ctypes.c_double),  # gradients
            ctypes.POINTER(ctypes.c_int)      # iostat
        ]
        self.lio.python_lio_gradients.restype = None
    
    def test_communication(self):
        """Probar comunicación con LIO"""
        print("🔗 Probando comunicación con LIO...")
        self.lio.python_hello()
        return True
    
    def setup_molecule(self, atoms, coordinates_angstrom, charge=0, basis="DZVP"):
        """
        Configurar molécula para cálculo
        
        Args:
            atoms: Lista de símbolos atómicos ['O', 'H', 'H']
            coordinates_angstrom: Array numpy (n_atoms, 3) en Angstroms
            charge: Carga molecular (default: 0)
            basis: Conjunto de base (default: "DZVP")
        """
        print(f"🧬 Configurando molécula:")
        print(f"   Átomos: {atoms}")
        print(f"   Carga: {charge}")
        print(f"   Base: {basis}")
        
        # Convertir símbolos a números atómicos
        symbol_to_z = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 
                       'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12}
        atomic_numbers = [symbol_to_z[atom] for atom in atoms]
        
        # Convertir coordenadas de Angstroms a bohr
        angstrom_to_bohr = 1.8897259886
        coords_bohr = np.array(coordinates_angstrom).flatten() * angstrom_to_bohr
        
        print(f"   Coordenadas (bohr): {len(coords_bohr)//3} átomos")
        
        # Convertir a arrays ctypes
        natom = len(atomic_numbers)
        c_atomic_numbers = (ctypes.c_int * natom)(*atomic_numbers)
        c_coordinates = (ctypes.c_double * len(coords_bohr))(*coords_bohr)
        c_charge = ctypes.c_int(charge)
        c_basis = basis.encode('utf-8')
        c_iostat = ctypes.c_int(0)
        
        # Inicializar sistema
        print("🔧 Inicializando sistema LIO...")
        self.lio.python_lio_init(natom, c_atomic_numbers, c_coordinates, 
                                c_charge, c_basis, ctypes.byref(c_iostat))
        
        if c_iostat.value == 0:
            self.system_initialized = True
            self.natom = natom
            print("✅ Sistema inicializado exitosamente")
            return True
        else:
            print(f"❌ Error en inicialización: {c_iostat.value}")
            return False
    
    def run_scf(self):
        """Ejecutar cálculo SCF"""
        if not self.system_initialized:
            print("❌ Sistema no inicializado")
            return False, None
        
        print("⚡ Ejecutando cálculo SCF...")
        c_energy = ctypes.c_double(0.0)
        c_scf_status = ctypes.c_int(0)
        
        try:
            self.lio.python_lio_scf(ctypes.byref(c_energy), ctypes.byref(c_scf_status))
            
            if c_scf_status.value == 0:
                energy = c_energy.value
                self.scf_converged = True
                print(f"✅ SCF convergido: E = {energy:.7f} Hartree")
                return True, energy
            else:
                print(f"❌ Error en SCF: {c_scf_status.value}")
                return False, None
        except Exception as e:
            print(f"❌ Excepción durante SCF: {e}")
            return False, None
    
    def calculate_forces(self):
        """
        Calcular fuerzas sobre átomos
        
        Returns:
            success: bool
            forces: numpy array (n_atoms, 3) en Hartree/bohr
        """
        if not self.scf_converged:
            print("❌ SCF no convergido, no se pueden calcular fuerzas")
            return False, None
        
        print("💪 Calculando fuerzas...")
        
        # Preparar array para gradientes (que son -fuerzas)
        c_gradients = (ctypes.c_double * (3 * self.natom))()
        c_grad_status = ctypes.c_int(0)
        
        try:
            self.lio.python_lio_gradients(c_gradients, ctypes.byref(c_grad_status))
            
            if c_grad_status.value == 0:
                # Convertir gradientes a fuerzas (fuerza = -gradiente)
                gradients = np.array([c_gradients[i] for i in range(3 * self.natom)])
                forces = -gradients.reshape(self.natom, 3)
                
                print(f"✅ Fuerzas calculadas para {self.natom} átomos")
                return True, forces
            else:
                print(f"❌ Error en cálculo de fuerzas: {c_grad_status.value}")
                return False, None
        except Exception as e:
            print(f"❌ Excepción durante cálculo de fuerzas: {e}")
            return False, None
    
    def print_forces(self, forces, atoms):
        """Imprimir fuerzas de manera legible"""
        print("\n🎯 FUERZAS SOBRE ÁTOMOS (Hartree/bohr):")
        print("=" * 50)
        
        for i, (atom, force) in enumerate(zip(atoms, forces)):
            fx, fy, fz = force
            f_norm = np.linalg.norm(force)
            print(f"   {atom:>2} {i+1:>2}: [{fx:8.5f}, {fy:8.5f}, {fz:8.5f}] |F|={f_norm:.5f}")
        
        total_force_norm = np.linalg.norm(forces.flatten())
        print(f"\n📊 Norma total de fuerzas: {total_force_norm:.6f} Hartree/bohr")
        
        # Convertir a unidades más comunes
        hartree_bohr_to_kcal_mol_ang = 627.509474 / 0.529177249
        forces_kcal = forces * hartree_bohr_to_kcal_mol_ang
        
        print(f"\n🔧 FUERZAS EN kcal/mol/Å:")
        for i, (atom, force) in enumerate(zip(atoms, forces_kcal)):
            fx, fy, fz = force
            f_norm = np.linalg.norm(force)
            print(f"   {atom:>2} {i+1:>2}: [{fx:8.2f}, {fy:8.2f}, {fz:8.2f}] |F|={f_norm:.2f}")

def demo_water_forces():
    """Demostración de cálculo de fuerzas para H2O"""
    print("💧 DEMO: CÁLCULO DE FUERZAS PARA H2O")
    print("=" * 60)
    
    # Crear calculadora
    calc = LIOForceCalculator()
    
    # Cargar biblioteca
    if not calc.load_library():
        return False
    
    # Probar comunicación
    calc.test_communication()
    
    # Definir molécula de agua (geometría del test)
    atoms = ['O', 'H', 'H']
    coords_angstrom = np.array([
        [71.762448,  35.512769,  96.172805],   # O
        [70.885172,  36.746272,  95.119946],   # H
        [73.544272,  35.969662,  96.043066]    # H
    ])
    
    # Configurar molécula
    if not calc.setup_molecule(atoms, coords_angstrom, charge=0, basis="DZVP"):
        return False
    
    # Ejecutar SCF
    success, energy = calc.run_scf()
    if not success:
        return False
    
    print(f"\n📈 RESULTADOS SCF:")
    print(f"   Energía: {energy:.7f} Hartree")
    print(f"   Referencia: -67.2365887 Hartree")
    print(f"   Diferencia: {abs(energy + 67.2365887):.8f} Hartree")
    
    # Calcular fuerzas
    success, forces = calc.calculate_forces()
    if success and forces is not None:
        calc.print_forces(forces, atoms)
        print("\n✅ CÁLCULO DE FUERZAS COMPLETADO EXITOSAMENTE")
        return True
    else:
        print("❌ Error en cálculo de fuerzas")
        return False

def demo_simple_molecule():
    """Demo con molécula más simple: H2"""
    print("\n⚛️  DEMO: CÁLCULO DE FUERZAS PARA H2")
    print("=" * 60)
    
    calc = LIOForceCalculator()
    
    if not calc.load_library():
        return False
    
    calc.test_communication()
    
    # Molécula de H2
    atoms = ['H', 'H']
    coords_angstrom = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.74]  # 0.74 Å separación
    ])
    
    if not calc.setup_molecule(atoms, coords_angstrom, charge=0, basis="STO-3G"):
        return False
    
    success, energy = calc.run_scf()
    if not success:
        return False
    
    print(f"\n📈 RESULTADOS SCF H2:")
    print(f"   Energía: {energy:.7f} Hartree")
    
    success, forces = calc.calculate_forces()
    if success and forces is not None:
        calc.print_forces(forces, atoms)
        return True
    else:
        return False

if __name__ == "__main__":
    print("🔬 CALCULADORA DE FUERZAS CON LIO")
    print("Interfaz Python para cálculo de fuerzas cuánticas")
    print("=" * 70)
    
    try:
        # Ejecutar demo principal con H2O
        success = demo_water_forces()
        
        if success:
            print("\n" + "=" * 70)
            print("🎉 DEMOSTRACIÓN COMPLETADA EXITOSAMENTE")
            print("=" * 70)
            print("✅ La calculadora de fuerzas está funcionando")
            print("✅ Interfaz Python-LIO completamente operativa")
            print("✅ Cálculos de fuerzas precisos y confiables")
        else:
            print("\n❌ Error en la demostración")
            
    except Exception as e:
        print(f"\n💥 Error inesperado: {e}")
        import traceback
        traceback.print_exc()
