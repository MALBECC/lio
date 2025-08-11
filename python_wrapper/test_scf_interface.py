#!/usr/bin/env python3
"""
Test script para la interfaz SCF de LIO
"""

import ctypes
import numpy as np

class LIOSCFInterface:
    def __init__(self, lib_path):
        self.lib = ctypes.CDLL(lib_path)
        
        # Configurar funciones
        self.lib.python_lio_init.argtypes = [
            ctypes.c_int,                    # natom
            ctypes.POINTER(ctypes.c_int),    # atomic_numbers
            ctypes.POINTER(ctypes.c_double), # coordinates
            ctypes.c_int,                    # charge
            ctypes.c_char_p,                 # basis_name
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_init.restype = None
        
        self.lib.python_lio_scf.argtypes = [
            ctypes.POINTER(ctypes.c_double), # total_energy
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_scf.restype = None
        
        self.lib.python_lio_gradients.argtypes = [
            ctypes.POINTER(ctypes.c_double), # gradients (3*natom)
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_gradients.restype = None
        
        self.lib.python_hello.argtypes = []
        self.lib.python_hello.restype = None
    
    def hello(self):
        """Función de prueba"""
        self.lib.python_hello()
    
    def init_system(self, atomic_numbers, coordinates, charge=0, basis_set="DZVP"):
        """Inicializar sistema LIO"""
        natom = len(atomic_numbers)
        
        # Convertir a arrays de ctypes
        atomic_nums = (ctypes.c_int * natom)(*atomic_numbers)
        coords = (ctypes.c_double * (3*natom))(*coordinates)
        iostat = ctypes.c_int()
        
        self.lib.python_lio_init(
            natom,
            atomic_nums,
            coords,
            charge,
            basis_set.encode('utf-8'),
            ctypes.byref(iostat)
        )
        
        return iostat.value
    
    def scf_calculation(self):
        """Ejecutar cálculo SCF"""
        energy = ctypes.c_double()
        iostat = ctypes.c_int()
        
        self.lib.python_lio_scf(
            ctypes.byref(energy),
            ctypes.byref(iostat)
        )
        
        return energy.value, iostat.value
    
    def calculate_gradients(self, natom):
        """Calcular gradientes después de SCF"""
        gradients_array = (ctypes.c_double * (3*natom))()
        iostat = ctypes.c_int()
        
        self.lib.python_lio_gradients(
            gradients_array,
            ctypes.byref(iostat)
        )
        
        # Convertir a lista Python
        gradients = [gradients_array[i] for i in range(3*natom)]
        return gradients, iostat.value

def main():
    print("🧬 TEST INTERFAZ SCF + GRADIENTES PARA LIO")
    print("=" * 60)
    
    # Cargar biblioteca
    lib_path = "/home/nano/lio/build/lib/liblio-g2g.so"
    lio = LIOSCFInterface(lib_path)
    print(f"✅ Biblioteca cargada: {lib_path}")
    
    # Prueba básica
    print("\n1️⃣ Prueba de comunicación:")
    lio.hello()
    
    # Configurar molécula H2O con geometría específica del test
    print("\n2️⃣ Configurando sistema H2O:")
    atomic_numbers = [8, 1, 1]  # O, H, H
    
    # Coordenadas en Angstroms del test
    coords_angstrom = [
        71.762448,  35.512769,  96.172805,   # O
        70.885172,  36.746272,  95.119946,   # H
        73.544272,  35.969662,  96.043066    # H
    ]
    
    # Convertir de Angstroms a bohr (LIO usa unidades atómicas)
    angstrom_to_bohr = 1.8897259886
    coordinates = [coord * angstrom_to_bohr for coord in coords_angstrom]
    
    print(f"   Átomos: {atomic_numbers}")
    print(f"   Coordenadas (Å): {coords_angstrom}")
    print(f"   Coordenadas (bohr): {[f'{coord:.6f}' for coord in coordinates]}")
    
    # Inicializar
    init_status = lio.init_system(atomic_numbers, coordinates, charge=0, basis_set="DZVP")
    print(f"📊 Estado de inicialización: {init_status}")
    
    if init_status == 0:
        print("\n3️⃣ Ejecutando cálculo SCF:")
        try:
            energy, scf_status = lio.scf_calculation()
            print(f"📊 Estado SCF: {scf_status}")
            
            if scf_status == 0:
                print(f"⚡ Energía total SCF: {energy:.8f} Hartree")
                print(f"⚡ Energía en eV: {energy * 27.211385:.4f} eV")
                print("\n🎉 ¡CÁLCULO SCF EXITOSO!")
                
                print("\n4️⃣ Calculando gradientes:")
                try:
                    gradients, grad_status = lio.calculate_gradients(len(atomic_numbers))
                    print(f"📊 Estado gradientes: {grad_status}")
                    
                    if grad_status == 0:
                        # Mostrar gradientes por átomo
                        print("🧮 Gradientes por átomo (Hartree/bohr):")
                        for i in range(len(atomic_numbers)):
                            idx = 3*i
                            gx, gy, gz = gradients[idx], gradients[idx+1], gradients[idx+2]
                            gnorm = (gx*gx + gy*gy + gz*gz)**0.5
                            print(f"   Átomo {i+1}: [{gx:8.5f}, {gy:8.5f}, {gz:8.5f}] |g|={gnorm:.5f}")
                        
                        # Norma total del gradiente
                        total_norm = sum(g*g for g in gradients)**0.5
                        print(f"🎯 Norma total del gradiente: {total_norm:.6f} Hartree/bohr")
                        
                        print("\n🎉 ¡CÁLCULO DE GRADIENTES EXITOSO!")
                    else:
                        print(f"❌ Error en cálculo de gradientes: {grad_status}")
                except Exception as e:
                    print(f"❌ Excepción durante gradientes: {e}")
            else:
                print(f"❌ Error en cálculo SCF: {scf_status}")
        except Exception as e:
            print(f"❌ Excepción durante SCF: {e}")
    else:
        print(f"❌ Error en inicialización: {init_status}")

if __name__ == "__main__":
    main()
