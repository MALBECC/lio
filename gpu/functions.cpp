#include <iostream>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include "common.h"
#include "init.h"
#include "cuda/double.h"
#include "matrix.h"
#include "cuda/exchnum.h"
#include "cubes.h"
using namespace std;
using namespace G2G;

void assign_significative_functions(LittleCube& cube, const double3& cube_coord, const vector<double>& min_exps)
{
	uint func = 0;

	HostMatrix<double> atom_cube_dists(fortran_vars.atoms);
	for (uint i = 0; i < fortran_vars.atoms; i++) {
		const double3& atom_dist = fortran_vars.atom_positions.get(i);
		double3 dist_vec;

    for (uint j = 0; j < 3; j++) {
      if (atom_dist[j] < cube_coord[j])
        dist_vec[j] = cube_coord[j] - atom_dist[j];
      else if (atom_dist[j] > (cube_coord[j] + little_cube_size))
        dist_vec[j] = atom_dist[j] - (cube_coord[j] + little_cube_size);
      else
        dist_vec[j] = 0;
    }
		
		double len = dist_vec.length();
		atom_cube_dists.get(i) = len * len;
	}
	
	/** S **/
	while (func < fortran_vars.s_funcs) {
    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all_functions || (min_exps[func] * atom_cube_dists.get(atom_nuc)) < max_function_exponent) {
			cube.functions.insert(func); cube.s_functions++;
			cube.nucleii.insert(atom_nuc);
		}
		func++;
	}
	
	/** P **/
	while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) {
    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all_functions || (min_exps[func] * atom_cube_dists.get(atom_nuc)) < max_function_exponent) {
			cube.functions.insert(func); cube.p_functions++;
			cube.nucleii.insert(atom_nuc);
		}

		func += 3;
	}
	
	/** D **/
	while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3 + fortran_vars.d_funcs * 6) {
    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all_functions || (min_exps[func] * atom_cube_dists.get(atom_nuc)) < max_function_exponent) {
			cube.functions.insert(func); cube.d_functions++; 
			cube.nucleii.insert(atom_nuc);
		}
		func += 6;
	}
}
