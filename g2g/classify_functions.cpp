#include <iostream>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include "common.h"
#include "init.h"
#include "cuda/cuda_extra.h"
#include "matrix.h"
#include "partition.h"
using namespace std;
using namespace G2G;


/*******************************
 * Cube
 *******************************/

void Cube::assign_significative_functions(const double3& cube_coord, const vector<double>& min_exps, const vector<double>& min_coeff)
{
  uint func = 0;

  HostMatrix<double> atom_cube_dists(fortran_vars.atoms);
  for (uint i = 0; i < fortran_vars.atoms; i++) {
    const double3& atom_pos = fortran_vars.atom_positions(i);
    double3 dist_vec;

    for (uint j = 0; j < 3; j++) {
      if (elem(atom_pos,j) < elem(cube_coord,j))
        elem(dist_vec,j) = elem(cube_coord,j) - elem(atom_pos,j);
      else if (elem(atom_pos,j) > (elem(cube_coord,j) + little_cube_size))
        elem(dist_vec,j) = elem(atom_pos,j) - (elem(cube_coord,j) + little_cube_size);
      else
        elem(dist_vec,j) = 0;
    }
    
    atom_cube_dists(i) = length2(dist_vec);
  }

  set<uint> functions_set;
  set<uint> nucleii_set;
  
  /** S **/
  while (func < fortran_vars.s_funcs) {
    uint atom_nuc = fortran_vars.nucleii(func) - 1;
    if (assign_all_functions || is_significative(FUNCTION_S, min_exps[func], min_coeff[func], atom_cube_dists(atom_nuc))) {
      functions_set.insert(func); s_functions++;
      nucleii_set.insert(atom_nuc);
    }
    func++;
  }
  
  /** P **/
  while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) {
    uint atom_nuc = fortran_vars.nucleii(func) - 1;
    if (assign_all_functions || is_significative(FUNCTION_P, min_exps[func], min_coeff[func], atom_cube_dists(atom_nuc))) {
      functions_set.insert(func); p_functions++;
      nucleii_set.insert(atom_nuc);
    }
    func += 3;
  }
  
  /** D **/
  while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3 + fortran_vars.d_funcs * 6) {
    uint atom_nuc = fortran_vars.nucleii(func) - 1;
    if (assign_all_functions || is_significative(FUNCTION_D, min_exps[func], min_coeff[func], atom_cube_dists(atom_nuc))) {
      functions_set.insert(func); d_functions++;
      nucleii_set.insert(atom_nuc);
    }
    func += 6;
  }

  local2global_func.resize(functions_set.size());
  copy(functions_set.begin(), functions_set.end(), local2global_func.begin());

  local2global_nuc.resize(nucleii_set.size());
  copy(nucleii_set.begin(), nucleii_set.end(), local2global_nuc.begin());

  PointGroup::compute_nucleii_maps();
}

/*****************************
 * Sphere
 *****************************/
void Sphere::assign_significative_functions(const std::vector<double>& min_exps, const std::vector<double>& min_coeff) {
   uint func = 0;

  // TODO: esto solo es necesario para los atomos en nucleii, idem arriba
  HostMatrix<double> atom_sphere_dists(fortran_vars.atoms);
  const double3& own_atom_pos = fortran_vars.atom_positions(atom);

  for (uint i = 0; i < fortran_vars.atoms; i++) {
    if (i == atom) atom_sphere_dists(i) = 0;
    else {
      const double3& atom_pos = fortran_vars.atom_positions(i);
      double3 dist_vec = (atom_pos - own_atom_pos);
      double dist_to_atom = length(dist_vec);
      double dist = (radius > dist_to_atom ? 0 : dist_to_atom - radius);
      atom_sphere_dists(i) = dist * dist;
    }
  }

  set<uint> functions_set;
  set<uint> nucleii_set;

  /** S **/
  while (func < fortran_vars.s_funcs) {
    uint atom_nuc = fortran_vars.nucleii(func) - 1;
    if (assign_all_functions || is_significative(FUNCTION_S, min_exps[func], min_coeff[func], atom_sphere_dists(atom_nuc))) {
      functions_set.insert(func); s_functions++;
      nucleii_set.insert(atom_nuc);
    }
    func++;
  }

  /** P **/
  while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) {
    uint atom_nuc = fortran_vars.nucleii(func) - 1;
    if (assign_all_functions || is_significative(FUNCTION_P, min_exps[func], min_coeff[func], atom_sphere_dists(atom_nuc))) {
      functions_set.insert(func); p_functions++;
      nucleii_set.insert(atom_nuc);
    }

    func += 3;
  }

  /** D **/
  while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3 + fortran_vars.d_funcs * 6) {
    uint atom_nuc = fortran_vars.nucleii(func) - 1;
    if (assign_all_functions || is_significative(FUNCTION_D, min_exps[func], min_coeff[func], atom_sphere_dists(atom_nuc))) {
      functions_set.insert(func); d_functions++;
      nucleii_set.insert(atom_nuc);
    }
    func += 6;
  }

  local2global_func.resize(functions_set.size());
  copy(functions_set.begin(), functions_set.end(), local2global_func.begin());

  local2global_nuc.resize(nucleii_set.size());
  copy(nucleii_set.begin(), nucleii_set.end(), local2global_nuc.begin());

  PointGroup::compute_nucleii_maps();
}
