#!/usr/bin/env python3
"""
Test simple para debuggear el problema del basis_set path
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

def debug_basis_path():
    """Debug del problema de path del basis set"""
    
    lib_path = find_lio_library()
    if not lib_path:
        return False
        
    lib = ctypes.CDLL(lib_path)
    
    # Setup function
    lib.python_lio_init.argtypes = [
        ctypes.c_int,                    # natom
        ctypes.POINTER(ctypes.c_int),    # atomic_numbers
        ctypes.POINTER(ctypes.c_double), # coordinates
        ctypes.c_int,                    # charge
        ctypes.c_char_p,                 # basis_name
        ctypes.POINTER(ctypes.c_int)     # iostat
    ]
    lib.python_lio_init.restype = None
    
    # Molécula de agua como en el test original
    atoms = [8, 1, 1]  # O, H, H
    coordinates = [
        [0.0000,  0.0000,  0.1173],  # O
        [0.0000,  0.7572, -0.4692],  # H
        [0.0000, -0.7572, -0.4692]   # H
    ]
    
    natom = len(atoms)
    atomic_numbers = (ctypes.c_int * natom)(*atoms)
    coords_flat = []
    for coord in coordinates:
        coords_flat.extend(coord)
    coords_array = (ctypes.c_double * (3 * natom))(*coords_flat)
    iostat = ctypes.c_int(0)
    
    print("🔧 Probando con H2O y basis set 'DZVP'")
    
    try:
        lib.python_lio_init(
            natom,
            atomic_numbers,
            coords_array,
            0,  # charge
            b"DZVP",  # basis_name
            ctypes.byref(iostat)
        )
        
        if iostat.value == 0:
            print("✅ Inicialización exitosa con H2O y DZVP")
            return True
        else:
            print(f"❌ Error en inicialización: {iostat.value}")
            return False
            
    except Exception as e:
        print(f"❌ Excepción: {e}")
        return False

if __name__ == "__main__":
    success = debug_basis_path()
    if success:
        print("✅ Debug completado exitosamente")
    else:
        print("❌ Debug falló")
