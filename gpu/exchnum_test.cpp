#include <cmath>
#include <iostream>
#include <cuda_runtime.h>
#include "timer.h"
#include "matrix.h"
#include "cuda/exchnum.h"
using namespace std;
using namespace G2G;

int main(void) {
	unsigned int atoms = 32;
	unsigned int grid_type = 1;
	
	HostMatrixFloat3 atom_positions(atoms);
	HostMatrixUInt types(atoms);
	HostMatrixFloat energy(1);
	
	unsigned int points = (grid_type == 1 ? EXCHNUM_GRIDA_SIZE : EXCHNUM_GRIDB_SIZE);
	HostMatrixFloat3 point_positions(points);
	HostMatrixFloat wang(points);
	
	for (unsigned int i = 0; i < atoms; i++) { atom_positions.data[i] = make_float3(i, i + 1, i + 2); types.data[i] = 1; }
	for (unsigned int i = 0; i < points; i++) { point_positions.data[i] = make_float3(1, 1, 1); wang.data[i] = 0.5; }	

	Timer total, actual;
	for (unsigned int it = 0; it < 100; it++) {
		actual.start();
		calc_energy(atom_positions, types, grid_type, point_positions, energy, wang);
		actual.stop();
		cout << "Resultado: " << energy.data[0] << " (" << actual << ")" << endl;
	}
}
