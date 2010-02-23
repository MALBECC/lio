#include <iostream>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include "common.h"
#include "init.h"
#include "cuda/double3.h"
#include "matrix.h"
#include "partition.h"
using namespace std;
using namespace G2G;


/*******************************
 * Cube
 *******************************/

void Cube::assign_significative_functions(const double3& cube_coord, const vector<double>& min_exps)
{
	uint func = 0;

  // TODO: esto solo es necesario para los atomos en nucleii, pero no se nucleii todavia, se podria ver si esto cuesta mas que calcularla on-the-fly
	HostMatrix<double> atom_cube_dists(fortran_vars.atoms);
	for (uint i = 0; i < fortran_vars.atoms; i++) {
		const double3& atom_pos = fortran_vars.atom_positions.get(i);
		double3 dist_vec;

    for (uint j = 0; j < 3; j++) {
      if (atom_pos[j] < cube_coord[j])
        dist_vec[j] = cube_coord[j] - atom_pos[j];
      else if (atom_pos[j] > (cube_coord[j] + little_cube_size))
        dist_vec[j] = atom_pos[j] - (cube_coord[j] + little_cube_size);
      else
        dist_vec[j] = 0;
    }
		
		double len = dist_vec.length();
		atom_cube_dists.get(i) = len * len;
	}

  set<uint> functions_set;
	
	/** S **/
	while (func < fortran_vars.s_funcs) {
    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all_functions || (min_exps[func] * atom_cube_dists.get(atom_nuc)) < max_function_exponent) {
			functions_set.insert(func); s_functions++;
			nucleii.insert(atom_nuc);
      nuc_map.push_back(atom_nuc);
		}
		func++;
	}
	
	/** P **/
	while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) {
    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all_functions || (min_exps[func] * atom_cube_dists.get(atom_nuc)) < max_function_exponent) {
			functions_set.insert(func); p_functions++;
			nucleii.insert(atom_nuc);
      nuc_map.push_back(atom_nuc);
		}
		func += 3;
	}
	
	/** D **/
	while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3 + fortran_vars.d_funcs * 6) {
    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all_functions || (min_exps[func] * atom_cube_dists.get(atom_nuc)) < max_function_exponent) {
			functions_set.insert(func); d_functions++;
			nucleii.insert(atom_nuc);
      nuc_map.push_back(atom_nuc);
		}
		func += 6;
	}

  functions.resize(functions_set.size());
  copy(functions_set.begin(), functions_set.end(), functions.begin());
}

/*****************************
 * Sphere
 *****************************/
void Sphere::assign_significative_functions(const std::vector<double>& min_exps) {
 	uint func = 0;

  // TODO: esto solo es necesario para los atomos en nucleii, idem arriba
	HostMatrix<double> atom_sphere_dists(fortran_vars.atoms);
  const double3& own_atom_pos = fortran_vars.atom_positions.get(atom);

	for (uint i = 0; i < fortran_vars.atoms; i++) {
    if (i == atom) atom_sphere_dists.get(i) = 0;
    else {
      const double3& atom_pos = fortran_vars.atom_positions.get(i);
      double3 dist_vec = (atom_pos - own_atom_pos);
			double dist_to_atom = dist_vec.length();
			assert(radius <= dist_to_atom);
			double dist = dist_to_atom - radius;
      atom_sphere_dists.get(i) = dist * dist;
    }
	}

  set<uint> functions_set;

	/** S **/
	while (func < fortran_vars.s_funcs) {
    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all_functions || (min_exps[func] * atom_sphere_dists.get(atom_nuc)) < max_function_exponent) {
			functions_set.insert(func); s_functions++;
			nucleii.insert(atom_nuc);
      nuc_map.push_back(atom_nuc);
		}
		func++;
	}

	/** P **/
	while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) {
    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all_functions || (min_exps[func] * atom_sphere_dists.get(atom_nuc)) < max_function_exponent) {
			functions_set.insert(func); p_functions++;
			nucleii.insert(atom_nuc);
      nuc_map.push_back(atom_nuc);
		}

		func += 3;
	}

	/** D **/
	while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3 + fortran_vars.d_funcs * 6) {
    uint atom_nuc = fortran_vars.nucleii.get(func) - 1;
		if (assign_all_functions || (min_exps[func] * atom_sphere_dists.get(atom_nuc)) < max_function_exponent) {
			functions_set.insert(func); d_functions++;
			nucleii.insert(atom_nuc);
      nuc_map.push_back(atom_nuc);
		}
		func += 6;
	}

  functions.resize(functions_set.size());
  copy(functions_set.begin(), functions_set.end(), functions.begin());
}
