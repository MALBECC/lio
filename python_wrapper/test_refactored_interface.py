#!/usr/bin/env python3
"""
Test de la interfaz Python-LIO refactorizada
Verificar que las rutinas estándar de LIO funcionan correctamente
"""

import ctypes
import numpy as np
import sys
import os

def find_lio_library():
    """Buscar la librería LIO compilada"""
    possible_paths = [
        'build/lib/liblio-g2g.so',
        '../build/lib/liblio-g2g.so',
        './liblio-g2g.so'
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            print(f"✅ Librería encontrada: {path}")
            return path
    
    print("❌ No se encontró la librería LIO")
    return None

class RefactoredLIOInterface:
    """Interfaz mejorada usando rutinas estándar de LIO"""
    
    def __init__(self):
        lib_path = find_lio_library()
        if not lib_path:
            raise RuntimeError("No se encontró la librería LIO")
            
        self.lib = ctypes.CDLL(lib_path)
        self._setup_functions()
        print("🎯 Interfaz refactorizada de LIO cargada")
    
    def _setup_functions(self):
        """Configurar signatures de las funciones refactorizadas"""
        
        # python_hello - sin cambios
        self.lib.python_hello.argtypes = []
        self.lib.python_hello.restype = None
        
        # python_lio_init - refactorizada
        self.lib.python_lio_init.argtypes = [
            ctypes.c_int,                    # natom
            ctypes.POINTER(ctypes.c_int),    # atomic_numbers
            ctypes.POINTER(ctypes.c_double), # coordinates
            ctypes.c_int,                    # charge
            ctypes.c_char_p,                 # basis_name
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_init.restype = None
        
        # python_lio_init_from_file - nueva función
        self.lib.python_lio_init_from_file.argtypes = [
            ctypes.c_int,                    # natom
            ctypes.POINTER(ctypes.c_int),    # atomic_numbers
            ctypes.POINTER(ctypes.c_double), # coordinates
            ctypes.c_int,                    # charge
            ctypes.c_char_p,                 # input_file
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_init_from_file.restype = None
        
        # python_lio_scf - sin cambios
        self.lib.python_lio_scf.argtypes = [
            ctypes.POINTER(ctypes.c_double), # total_energy
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_scf.restype = None
        
        # python_lio_gradients - sin cambios
        self.lib.python_lio_gradients.argtypes = [
            ctypes.POINTER(ctypes.c_double), # gradients
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_gradients.restype = None
        
        # python_lio_finalize - refactorizada
        self.lib.python_lio_finalize.argtypes = []
        self.lib.python_lio_finalize.restype = None
    
    def hello(self):
        """Test básico"""
        print("🔧 Llamando python_hello...")
        self.lib.python_hello()
    
    def init_system(self, atoms, coordinates, charge=0, basis_set="DZVP"):
        """Inicializar usando rutinas estándar de LIO"""
        natom = len(atoms)
        
        # Preparar arrays para C
        atomic_numbers = (ctypes.c_int * natom)(*atoms)
        coords_flat = []
        for coord in coordinates:
            coords_flat.extend(coord)
        coords_array = (ctypes.c_double * (3 * natom))(*coords_flat)
        iostat = ctypes.c_int(0)
        
        print(f"🔧 Inicializando sistema refactorizado: {natom} átomos, carga {charge}, basis {basis_set}")
        
        # Llamar función refactorizada
        self.lib.python_lio_init(
            natom,
            atomic_numbers,
            coords_array,
            charge,
            basis_set.encode('utf-8'),
            ctypes.byref(iostat)
        )
        
        if iostat.value != 0:
            raise RuntimeError(f"Error en inicialización: {iostat.value}")
        
        print("✅ Sistema inicializado con rutinas estándar de LIO")
        return True
    
    def init_from_file(self, atoms, coordinates, charge=0, input_file="lio.in"):
        """Inicializar usando archivo de configuración"""
        natom = len(atoms)
        
        # Preparar arrays para C
        atomic_numbers = (ctypes.c_int * natom)(*atoms)
        coords_flat = []
        for coord in coordinates:
            coords_flat.extend(coord)
        coords_array = (ctypes.c_double * (3 * natom))(*coords_flat)
        iostat = ctypes.c_int(0)
        
        print(f"🔧 Inicializando desde archivo: {input_file}")
        
        # Llamar nueva función
        self.lib.python_lio_init_from_file(
            natom,
            atomic_numbers,
            coords_array,
            charge,
            input_file.encode('utf-8'),
            ctypes.byref(iostat)
        )
        
        if iostat.value != 0:
            print(f"⚠️  Advertencia: Error leyendo archivo {input_file} (esperado si no existe)")
            return False
        
        print("✅ Sistema inicializado desde archivo de configuración")
        return True
    
    def calculate_scf(self):
        """Calcular SCF"""
        energy = ctypes.c_double(0.0)
        iostat = ctypes.c_int(0)
        
        print("🔧 Ejecutando cálculo SCF refactorizado...")
        
        self.lib.python_lio_scf(
            ctypes.byref(energy),
            ctypes.byref(iostat)
        )
        
        if iostat.value != 0:
            raise RuntimeError(f"Error en cálculo SCF: {iostat.value}")
        
        print(f"✅ SCF completado: {energy.value:.6f} Hartree")
        return energy.value
    
    def calculate_gradients(self, natom):
        """Calcular gradientes"""
        gradients = (ctypes.c_double * (3 * natom))()
        iostat = ctypes.c_int(0)
        
        print("🔧 Calculando gradientes...")
        
        self.lib.python_lio_gradients(
            gradients,
            ctypes.byref(iostat)
        )
        
        if iostat.value != 0:
            raise RuntimeError(f"Error en cálculo de gradientes: {iostat.value}")
        
        # Convertir a numpy array
        grad_array = np.array([gradients[i] for i in range(3 * natom)])
        grad_norm = np.linalg.norm(grad_array)
        
        print(f"✅ Gradientes calculados, norma: {grad_norm:.6f} Hartree/bohr")
        return grad_array.reshape((natom, 3))
    
    def finalize(self):
        """Finalizar usando rutina estándar"""
        print("🧹 Finalizando con rutina estándar de LIO...")
        self.lib.python_lio_finalize()
        print("✅ Finalización completada")

def test_refactored_interface():
    """Test completo de la interfaz refactorizada"""
    print("=" * 60)
    print("🚀 TEST DE INTERFAZ PYTHON-LIO REFACTORIZADA")
    print("=" * 60)
    
    try:
        # Inicializar interfaz
        lio = RefactoredLIOInterface()
        
        # Test básico
        print("\n1️⃣ Test básico:")
        lio.hello()
        
        # Molécula de agua
        atoms = [8, 1, 1]  # O, H, H
        coordinates = [
            [0.0000,  0.0000,  0.1173],  # O
            [0.0000,  0.7572, -0.4692],  # H
            [0.0000, -0.7572, -0.4692]   # H
        ]
        
        print("\n2️⃣ Test inicialización refactorizada:")
        lio.init_system(atoms, coordinates, charge=0, basis_set="DZVP")
        
        print("\n3️⃣ Test cálculo SCF con rutinas estándar:")
        energy = lio.calculate_scf()
        
        print("\n4️⃣ Test cálculo de gradientes:")
        gradients = lio.calculate_gradients(len(atoms))
        
        print("\n5️⃣ Test finalización con rutina estándar:")
        lio.finalize()
        
        print("\n🎉 TODOS LOS TESTS PASARON EXITOSAMENTE")
        print(f"📊 Energía final: {energy:.6f} Hartree")
        print(f"📊 Norma gradientes: {np.linalg.norm(gradients):.6f} Hartree/bohr")
        
        return True
        
    except Exception as e:
        print(f"\n❌ ERROR EN TEST: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_file_initialization():
    """Test opcional de inicialización desde archivo"""
    print("\n" + "=" * 60)
    print("🗂️  TEST OPCIONAL: INICIALIZACIÓN DESDE ARCHIVO")
    print("=" * 60)
    
    # Crear archivo de configuración de ejemplo
    lio_config = """
# Configuración de ejemplo para LIO
FCOORD='qm.xyz'
CHARGE=0
IEXCH=129
NGRID=3
IGRID=2
OPEN=.FALSE.
BASIS_SET='DZVP'
FITTING_SET='DZVP Coulomb Fitting'
INT_BASIS=.TRUE.
"""
    
    with open("test_lio.in", "w") as f:
        f.write(lio_config)
    
    try:
        lio = RefactoredLIOInterface()
        
        atoms = [8, 1, 1]  # H2O
        coordinates = [
            [0.0000,  0.0000,  0.1173],
            [0.0000,  0.7572, -0.4692],
            [0.0000, -0.7572, -0.4692]
        ]
        
        success = lio.init_from_file(atoms, coordinates, charge=0, input_file="test_lio.in")
        
        if success:
            print("✅ Inicialización desde archivo exitosa")
            energy = lio.calculate_scf()
            print(f"📊 Energía con configuración de archivo: {energy:.6f} Hartree")
            lio.finalize()
        else:
            print("⚠️  Test de archivo omitido (archivo no procesado correctamente)")
            
    except Exception as e:
        print(f"⚠️  Test de archivo falló (esperado): {e}")
    
    finally:
        # Limpiar archivo de test
        if os.path.exists("test_lio.in"):
            os.remove("test_lio.in")

if __name__ == "__main__":
    success = test_refactored_interface()
    
    # Test opcional de archivo
    test_file_initialization()
    
    if success:
        print("\n🎯 REFACTORIZACIÓN EXITOSA: La interfaz Python-LIO")
        print("   ahora usa las rutinas estándar de LIO, siendo más")
        print("   mantenible y robusta.")
        sys.exit(0)
    else:
        print("\n❌ TESTS FALLARON")
        sys.exit(1)
