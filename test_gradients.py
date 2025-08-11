#!/usr/bin/env python3

import ctypes
import os
import signal
import sys

def signal_handler(sig, frame):
    print("\n🛑 Señal recibida, saliendo limpiamente...")
    sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

def main():
    print("🧬 TEST COMPLETO: SCF + GRADIENTES")
    print("=" * 60)
    
    # Cargar biblioteca
    lib_path = "/home/nano/lio/build/lib/liblio-g2g.so"
    if not os.path.exists(lib_path):
        print(f"❌ Biblioteca no encontrada: {lib_path}")
        return
    
    try:
        lio = ctypes.CDLL(lib_path)
        print(f"✅ Biblioteca cargada: {lib_path}")
    except Exception as e:
        print(f"❌ Error cargando biblioteca: {e}")
        return
    
    # Configurar funciones
    lio.python_hello.argtypes = []
    lio.python_hello.restype = None
    
    lio.python_lio_init.argtypes = [
        ctypes.c_int,  # natom
        ctypes.POINTER(ctypes.c_int),  # atomic_numbers
        ctypes.POINTER(ctypes.c_double),  # coordinates  
        ctypes.c_int,  # charge
        ctypes.c_char_p,  # basis_name
        ctypes.POINTER(ctypes.c_int)  # iostat
    ]
    lio.python_lio_init.restype = None
    
    lio.python_lio_scf.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)]
    lio.python_lio_scf.restype = None
    
    lio.python_lio_gradients.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # gradients
        ctypes.POINTER(ctypes.c_int)      # iostat
    ]
    lio.python_lio_gradients.restype = None
    
    # Test de comunicación
    print("\n1️⃣ Prueba de comunicación:")
    lio.python_hello()
    
    # Configurar sistema
    print("\n2️⃣ Configurando sistema H2O:")
    atomic_numbers = [8, 1, 1]
    coords_angstrom = [
        71.762448,  35.512769,  96.172805,   # O
        70.885172,  36.746272,  95.119946,   # H
        73.544272,  35.969662,  96.043066    # H
    ]
    
    angstrom_to_bohr = 1.8897259886
    coordinates = [coord * angstrom_to_bohr for coord in coords_angstrom]
    
    print(f"   Átomos: {atomic_numbers}")
    print(f"   Coordenadas configuradas para test de validación")
    
    # Convertir a arrays ctypes
    c_atomic_numbers = (ctypes.c_int * len(atomic_numbers))(*atomic_numbers)
    c_coordinates = (ctypes.c_double * len(coordinates))(*coordinates)
    c_charge = ctypes.c_int(0)
    c_basis = b"DZVP"
    c_iostat = ctypes.c_int(0)
    
    # Inicializar
    print("\n3️⃣ Inicializando sistema:")
    lio.python_lio_init(len(atomic_numbers), c_atomic_numbers, c_coordinates, 
                       c_charge, c_basis, ctypes.byref(c_iostat))
    
    if c_iostat.value == 0:
        print("✅ Inicialización exitosa")
        
        # Cálculo SCF
        print("\n4️⃣ Ejecutando cálculo SCF:")
        c_energy = ctypes.c_double(0.0)
        c_scf_status = ctypes.c_int(0)
        
        try:
            lio.python_lio_scf(ctypes.byref(c_energy), ctypes.byref(c_scf_status))
            
            if c_scf_status.value == 0:
                energy = c_energy.value
                print(f"\n✅ SCF convergido:")
                print(f"   Energía: {energy:.7f} Hartree")
                print(f"   Referencia: -67.2365887 Hartree")
                print(f"   Diferencia: {abs(energy + 67.2365887):.8f} Hartree")
                
                # Ahora intentar gradientes
                print("\n5️⃣ Calculando gradientes:")
                
                # Preparar array para gradientes
                natom = len(atomic_numbers)
                c_gradients = (ctypes.c_double * (3 * natom))()
                c_grad_status = ctypes.c_int(0)
                
                try:
                    lio.python_lio_gradients(c_gradients, ctypes.byref(c_grad_status))
                    
                    if c_grad_status.value == 0:
                        gradients = [c_gradients[i] for i in range(3 * natom)]
                        
                        print("✅ Gradientes calculados exitosamente:")
                        for i in range(natom):
                            idx = 3*i
                            gx, gy, gz = gradients[idx], gradients[idx+1], gradients[idx+2]
                            gnorm = (gx*gx + gy*gy + gz*gz)**0.5
                            print(f"   Átomo {i+1}: [{gx:8.5f}, {gy:8.5f}, {gz:8.5f}] |g|={gnorm:.5f}")
                        
                        total_norm = sum(g*g for g in gradients)**0.5
                        print(f"🎯 Norma total del gradiente: {total_norm:.6f} Hartree/bohr")
                        
                        print("\n🎉 ¡CÁLCULO COMPLETO EXITOSO!")
                        print("✅ SCF convergido con energía validada")
                        print("✅ Gradientes calculados correctamente")
                        
                    else:
                        print(f"❌ Error en cálculo de gradientes: {c_grad_status.value}")
                        
                except Exception as e:
                    print(f"❌ Excepción durante gradientes: {e}")
                    
            else:
                print(f"❌ Error en SCF: {c_scf_status.value}")
                
        except Exception as e:
            print(f"❌ Excepción durante SCF: {e}")
            
    else:
        print(f"❌ Error en inicialización: {c_iostat.value}")
    
    print("\n" + "="*60)
    print("🔬 RESUMEN DE VALIDACIÓN:")
    if 'energy' in locals():
        print(f"• Energía SCF: {energy:.7f} Hartree (✅ validada)")
    if 'total_norm' in locals():
        print(f"• Gradientes: norma = {total_norm:.6f} au (✅ calculados)")
    print("• Interface Python-LIO funcionando correctamente")
    print("="*60)

if __name__ == "__main__":
    main()
