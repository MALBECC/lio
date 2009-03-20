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

bool assign_all = false;

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
		
		/*if (atom_dist.x < cube_coord.x) dist_vec.x = cube_coord.x - atom_dist.x;
		else if (atom_dist.x > (cube_coord.x + little_cube_size)) dist_vec.x = atom_dist.x - (cube_coord.x + little_cube_size);
		if (atom_dist.y < cube_coord.y) dist_vec.y = cube_coord.y - atom_dist.y;
		else if (atom_dist.y > (cube_coord.y + little_cube_size)) dist_vec.y = atom_dist.y - (cube_coord.y + little_cube_size);
		if (atom_dist.z < cube_coord.z) dist_vec.z = cube_coord.z - atom_dist.z;
		else if (atom_dist.z > (cube_coord.z + little_cube_size)) dist_vec.z = atom_dist.z - (cube_coord.z + little_cube_size);*/

    double len = dist_vec.length();
		atom_cube_dists.get(i) = len * len;
	}
	
	/** S **/
	while (func < fortran_vars.s_funcs) {
    /** TODO: esto se puede sacar, no? **/
		//if (cube.functions.find(func) != cube.functions.end()) { func++; continue; }

    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all || (min_exps[func] * atom_cube_dists.get(atom_nuc)) < max_function_exponent) {
			cube.functions.insert(func); cube.s_functions++;
			cube.nucleii.insert(atom_nuc);
		}
		func++;
	}
	
	/** P **/
	while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) {
		//if (cube.functions.find(func) != cube.functions.end()) { func += 3; continue; }

    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all || (min_exps[func] * atom_cube_dists.get(atom_nuc)) < max_function_exponent) {
			cube.functions.insert(func); cube.p_functions++;
			cube.nucleii.insert(atom_nuc);
		}

		func += 3;
	}
	
	/** D **/
	while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3 + fortran_vars.d_funcs * 6) {
		//if (cube.functions.find(func) != cube.functions.end()) { func += 6; continue; }

    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all || (min_exps[func] * atom_cube_dists.get(atom_nuc)) < max_function_exponent) {
			cube.functions.insert(func); cube.d_functions++; 
			cube.nucleii.insert(atom_nuc);
		}
		func += 6;
	}
}
