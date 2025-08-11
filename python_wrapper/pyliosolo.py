#!/usr/bin/env python3
"""
pyliosolo - Versión Python de liosolo usando rutinas refactorizadas de LIO
Replica la funcionalidad completa de liosolo usando la interfaz Python-LIO optimizada.
"""

import ctypes
import numpy as np
import sys
import os
import argparse
import time
from pathlib import Path

class PyliosoloInterface:
    """Interfaz completa Python-LIO equivalente a liosolo"""
    
    def __init__(self):
        """Inicializar interfaz"""
        self.lib_path = self._find_lio_library()
        if not self.lib_path:
            raise RuntimeError("❌ No se encontró la librería LIO")
            
        self.lib = ctypes.CDLL(self.lib_path)
        self._setup_functions()
        
        # Variables del sistema
        self.natom = 0
        self.atomic_numbers = []
        self.coordinates = []
        self.charge = 0
        self.basis_set = "DZVP"
        self.verbose = False
        
        print("🎯 pyliosolo: Interfaz Python-LIO cargada")
    
    def _find_lio_library(self):
        """Buscar la librería LIO compilada"""
        possible_paths = [
            'build/lib/liblio-g2g.so',
            '../build/lib/liblio-g2g.so',
            './liblio-g2g.so',
            '/home/nano/lio/build/lib/liblio-g2g.so'
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                print(f"✅ Librería LIO encontrada: {path}")
                return path
        
        return None
    
    def _setup_functions(self):
        """Configurar signatures de las funciones LIO"""
        
        # python_lio_init_from_file - para leer archivos .in
        self.lib.python_lio_init_from_file.argtypes = [
            ctypes.c_int,                    # natom
            ctypes.POINTER(ctypes.c_int),    # atomic_numbers
            ctypes.POINTER(ctypes.c_double), # coordinates
            ctypes.c_int,                    # charge
            ctypes.c_char_p,                 # input_file
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_init_from_file.restype = None
        
        # python_lio_scf
        self.lib.python_lio_scf.argtypes = [
            ctypes.POINTER(ctypes.c_double), # total_energy
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_scf.restype = None
        
        # python_lio_gradients
        self.lib.python_lio_gradients.argtypes = [
            ctypes.POINTER(ctypes.c_double), # gradients
            ctypes.POINTER(ctypes.c_int)     # iostat
        ]
        self.lib.python_lio_gradients.restype = None
        
        # python_lio_finalize
        self.lib.python_lio_finalize.argtypes = []
        self.lib.python_lio_finalize.restype = None
    
    def parse_arguments(self):
        """Parsear argumentos de línea de comandos (equivalente a liosolo)"""
        parser = argparse.ArgumentParser(
            description='pyliosolo - Versión Python de liosolo',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='''
Ejemplos de uso:
  pyliosolo.py -i agua.in -c agua.xyz
  pyliosolo.py -i benzene.in -c benzene.xyz -v
  pyliosolo.py -i complex.in -c complex.xyz -b "6-31G*"
            '''
        )
        
        parser.add_argument('-i', '--input', 
                          help='Archivo de input (.in)', required=True)
        parser.add_argument('-c', '--coords', 
                          help='Archivo de coordenadas (.xyz)', required=True)
        parser.add_argument('-b', '--basis', 
                          help='Basis set (por defecto: DZVP)', default='DZVP')
        parser.add_argument('-v', '--verbose', 
                          action='store_true', help='Modo verbose')
        
        args = parser.parse_args()
        
        # Verificar que los archivos existen
        if not os.path.exists(args.input):
            print(f"❌ Error: Archivo de input no encontrado: {args.input}")
            sys.exit(1)
        
        if not os.path.exists(args.coords):
            print(f"❌ Error: Archivo de coordenadas no encontrado: {args.coords}")
            sys.exit(1)
        
        return args
    
    def read_xyz_file(self, xyz_file):
        """Leer archivo XYZ (formato estándar o formato liosolo)"""
        print(f"📂 Leyendo coordenadas desde: {xyz_file}")
        
        atomic_symbols_to_numbers = {
            'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
            'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
            'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
            'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
            'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36
        }
        
        atomic_numbers = []
        coordinates = []
        
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        
        # Detectar formato automáticamente
        first_line = lines[0].strip()
        
        # Formato liosolo: sin número de átomos en primera línea
        # Formato estándar: número de átomos en primera línea
        if first_line.isdigit():
            # Formato estándar XYZ
            natom = int(first_line)
            start_line = 2  # Saltar línea de comentario
            file_format = "estándar"
            print(f"   🔍 Formato detectado: XYZ estándar ({natom} átomos)")
        else:
            # Formato liosolo: cada línea es un átomo
            natom = len([line for line in lines if line.strip()])
            start_line = 0  # Empezar desde la primera línea
            file_format = "liosolo"
            print(f"   🔍 Formato detectado: liosolo ({natom} átomos)")
        
        # Leer átomos según el formato
        for i in range(start_line, start_line + natom):
            if i >= len(lines):
                break
                
            parts = lines[i].strip().split()
            if len(parts) < 4:
                continue
            
            # Puede ser símbolo atómico o número atómico directo
            atom_identifier = parts[0]
            if atom_identifier.isdigit():
                # Es un número atómico directo (formato liosolo)
                atomic_number = int(atom_identifier)
                if self.verbose:
                    atomic_symbol = {v: k for k, v in atomic_symbols_to_numbers.items()}.get(atomic_number, f"Z{atomic_number}")
                    print(f"   📍 Átomo {len(atomic_numbers)+1}: Z={atomic_number} ({atomic_symbol})")
            else:
                # Es un símbolo atómico (formato estándar)
                if atom_identifier not in atomic_symbols_to_numbers:
                    raise ValueError(f"Símbolo atómico desconocido: {atom_identifier}")
                atomic_number = atomic_symbols_to_numbers[atom_identifier]
                if self.verbose:
                    print(f"   📍 Átomo {len(atomic_numbers)+1}: {atom_identifier} (Z={atomic_number})")
            
            atomic_numbers.append(atomic_number)
            
            # Coordenadas (en Angstroms, se convertirán a bohr después)
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            coordinates.append([x, y, z])
        
        print(f"   ✅ {len(atomic_numbers)} átomos leídos (formato {file_format})")
        print(f"   📋 Números atómicos: {atomic_numbers}")
        
        return len(atomic_numbers), atomic_numbers, coordinates
    
    def convert_angstrom_to_bohr(self, coordinates_angstrom):
        """Convertir coordenadas de Angstroms a bohrs"""
        angstrom_to_bohr = 1.8897259886
        coordinates_bohr = []
        
        for coord in coordinates_angstrom:
            coord_bohr = [c * angstrom_to_bohr for c in coord]
            coordinates_bohr.append(coord_bohr)
        
        if self.verbose:
            print(f"   🔄 Conversión Å → bohr (factor: {angstrom_to_bohr})")
        
        return coordinates_bohr
    
    def print_logo(self):
        """Imprimir logo de pyliosolo"""
        print("=" * 70)
        print("  ____        _ _                  _       ")
        print(" |  _ \\ _   _| (_) ___  ___  ___ | | ___  ")
        print(" | |_) | | | | | |/ _ \\/ __|/ _ \\| |/ _ \\ ")
        print(" |  __/| |_| | | | (_) \\__ \\ (_) | | (_) |")
        print(" |_|    \\__, |_|_|\\___/|___/\\___/|_|\\___/ ")
        print("        |___/                            ")
        print("")
        print("  🐍 Python version of liosolo using refactored LIO routines")
        print("  🧬 Quantum chemistry calculations with Gaussian basis sets")
        print("=" * 70)
    
    def print_system_info(self):
        """Imprimir información del sistema"""
        print(f"\n📊 INFORMACIÓN DEL SISTEMA:")
        print(f"   🧮 Átomos: {self.natom}")
        print(f"   ⚡ Carga: {self.charge}")
        print(f"   🔬 Basis set: {self.basis_set}")
        print(f"   📝 Verbose: {'Activado' if self.verbose else 'Desactivado'}")
    
    def run_calculation(self, input_file, coords_file, basis_set="DZVP", verbose=False):
        """Ejecutar cálculo completo (equivalente a liosolo)"""
        
        self.verbose = verbose
        self.basis_set = basis_set
        
        # Logo
        self.print_logo()
        
        # 1. Leer coordenadas
        self.natom, self.atomic_numbers, coords_angstrom = self.read_xyz_file(coords_file)
        
        # 2. Convertir unidades
        self.coordinates = self.convert_angstrom_to_bohr(coords_angstrom)
        
        # 3. Mostrar información del sistema
        self.print_system_info()
        
        # 4. Inicializar LIO con archivo de input
        print(f"\n🚀 INICIALIZANDO SISTEMA LIO")
        print(f"   📂 Input file: {input_file}")
        
        # Preparar arrays para C
        atomic_numbers_array = (ctypes.c_int * self.natom)(*self.atomic_numbers)
        coords_flat = []
        for coord in self.coordinates:
            coords_flat.extend(coord)
        coords_array = (ctypes.c_double * (3 * self.natom))(*coords_flat)
        iostat = ctypes.c_int(0)
        
        # Llamar inicialización desde archivo
        self.lib.python_lio_init_from_file(
            self.natom,
            atomic_numbers_array,
            coords_array,
            self.charge,
            input_file.encode('utf-8'),
            ctypes.byref(iostat)
        )
        
        if iostat.value != 0:
            raise RuntimeError(f"❌ Error en inicialización LIO: {iostat.value}")
        
        print("   ✅ Sistema LIO inicializado exitosamente")
        
        # 5. Ejecutar cálculo SCF
        print(f"\n⚡ EJECUTANDO CÁLCULO SCF")
        start_time = time.time()
        
        energy = ctypes.c_double(0.0)
        iostat = ctypes.c_int(0)
        
        self.lib.python_lio_scf(
            ctypes.byref(energy),
            ctypes.byref(iostat)
        )
        
        if iostat.value != 0:
            raise RuntimeError(f"❌ Error en cálculo SCF: {iostat.value}")
        
        scf_time = time.time() - start_time
        final_energy = energy.value
        
        print(f"   ✅ SCF convergido en {scf_time:.2f} segundos")
        print(f"   🎯 Energía total: {final_energy:.8f} Hartree")
        print(f"   🎯 Energía total: {final_energy * 27.211386:.4f} eV")
        
        # 6. Calcular gradientes/fuerzas
        print(f"\n🧮 CALCULANDO GRADIENTES")
        
        gradients_array = (ctypes.c_double * (3 * self.natom))()
        iostat = ctypes.c_int(0)
        
        self.lib.python_lio_gradients(
            gradients_array,
            ctypes.byref(iostat)
        )
        
        if iostat.value != 0:
            print(f"   ⚠️  Warning: Error en cálculo de gradientes: {iostat.value}")
            gradients = None
        else:
            # Convertir a numpy y mostrar resultados
            grad_array = np.array([gradients_array[i] for i in range(3 * self.natom)])
            gradients = grad_array.reshape((self.natom, 3))
            grad_norm = np.linalg.norm(grad_array)
            
            print(f"   ✅ Gradientes calculados exitosamente")
            print(f"   📊 Norma del gradiente: {grad_norm:.6f} Hartree/bohr")
            
            if verbose:
                print(f"\n   📋 Gradientes por átomo (Hartree/bohr):")
                for i, grad in enumerate(gradients):
                    atom_grad_norm = np.linalg.norm(grad)
                    print(f"      Átomo {i+1} (Z={self.atomic_numbers[i]}): "
                          f"[{grad[0]:8.5f}, {grad[1]:8.5f}, {grad[2]:8.5f}] "
                          f"|g|={atom_grad_norm:.5f}")
        
        # 7. Finalizar
        print(f"\n🧹 FINALIZANDO CÁLCULO")
        self.lib.python_lio_finalize()
        print("   ✅ Sistema LIO finalizado")
        
        # 8. Resumen final
        total_time = time.time() - start_time
        print(f"\n" + "=" * 60)
        print(f"🎉 CÁLCULO COMPLETADO EXITOSAMENTE")
        print(f"=" * 60)
        print(f"⏱️  Tiempo total: {total_time:.2f} segundos")
        print(f"🎯 Energía final: {final_energy:.8f} Hartree")
        if gradients is not None:
            print(f"📊 Norma gradientes: {np.linalg.norm(gradients):.6f} Hartree/bohr")
        print(f"=" * 60)
        
        return {
            'energy': final_energy,
            'gradients': gradients,
            'time': total_time,
            'natom': self.natom,
            'atomic_numbers': self.atomic_numbers
        }

def main():
    """Función principal - equivalente al programa liosolo"""
    
    try:
        # Crear interfaz
        pyliosolo = PyliosoloInterface()
        
        # Parsear argumentos
        args = pyliosolo.parse_arguments()
        
        # Ejecutar cálculo
        results = pyliosolo.run_calculation(
            input_file=args.input,
            coords_file=args.coords,
            basis_set=args.basis,
            verbose=args.verbose
        )
        
        print("\n✅ pyliosolo ejecutado exitosamente")
        return 0
        
    except Exception as e:
        print(f"\n❌ ERROR EN pyliosolo: {e}")
        if 'args' in locals() and hasattr(args, 'verbose') and args.verbose:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
