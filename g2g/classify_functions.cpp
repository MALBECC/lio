#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "cuda/cuda_extra.h"
using namespace std;
using namespace G2G;

template <class scalar_type>
void PointGroup<scalar_type>::assign_functions(
    HostMatrix<double> dists, const std::vector<double>& min_exps,
    const std::vector<double>& min_coeff) {
  uint func = 0;
  set<uint> functions_set;
  set<uint> nucleii_set;

  /** S **/
  while (func < fortran_vars.s_funcs) {
    uint atom_nuc = fortran_vars.nucleii(func) - 1;
    if (assign_all_functions ||
        is_significative(FUNCTION_S, min_exps[func], min_coeff[func],
                         dists(atom_nuc))) {
      functions_set.insert(func);
      s_functions++;
      nucleii_set.insert(atom_nuc);
    }
    func++;
  }

  /** P **/
  while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) {
    uint atom_nuc = fortran_vars.nucleii(func) - 1;
    if (assign_all_functions ||
        is_significative(FUNCTION_P, min_exps[func], min_coeff[func],
                         dists(atom_nuc))) {
      functions_set.insert(func);
      p_functions++;
      nucleii_set.insert(atom_nuc);
    }
    func += 3;
  }

  /** D **/
  while (func < fortran_vars.s_funcs + fortran_vars.p_funcs * 3 +
                    fortran_vars.d_funcs * 6) {
    uint atom_nuc = fortran_vars.nucleii(func) - 1;
    if (assign_all_functions ||
        is_significative(FUNCTION_D, min_exps[func], min_coeff[func],
                         dists(atom_nuc))) {
      functions_set.insert(func);
      d_functions++;
      nucleii_set.insert(atom_nuc);
    }
    func += 6;
  }

  local2global_func.resize(functions_set.size());
  copy(functions_set.begin(), functions_set.end(), local2global_func.begin());

  local2global_nuc.resize(nucleii_set.size());
  copy(nucleii_set.begin(), nucleii_set.end(), local2global_nuc.begin());
}

template <class scalar_type>
void PointGroup<scalar_type>::compute_nucleii_maps(void) {
  if (total_functions_simple() != 0) {
    func2global_nuc.resize(total_functions_simple());
    for (uint i = 0; i < total_functions_simple(); i++) {
      func2global_nuc(i) = fortran_vars.nucleii(local2global_func[i]) - 1;
    }

    func2local_nuc.resize(total_functions());
    uint ii = 0;
    for (uint i = 0; i < total_functions_simple(); i++) {
      uint global_atom = func2global_nuc(i);
      uint local_atom =
          std::distance(local2global_nuc.begin(),
                        std::find(local2global_nuc.begin(),
                                  local2global_nuc.end(), global_atom));
      uint inc = small_function_type(i);
      for (uint k = 0; k < inc; k++, ii++) func2local_nuc(ii) = local_atom;
    }
  }
}

/*******************************
 * Cube
 *******************************/
template <class scalar_type>
void PointGroup<scalar_type>::assign_functions_as_cube(
    const double3& cube_coord, const std::vector<double>& min_exps,
    const std::vector<double>& min_coeff) {
  HostMatrix<double> atom_cube_dists(fortran_vars.atoms);
  for (uint i = 0; i < fortran_vars.atoms; i++) {
    const double3& atom_pos = fortran_vars.atom_positions(i);
    double3 dist_vec;

    for (uint j = 0; j < 3; j++) {
      if (elem(atom_pos, j) < elem(cube_coord, j))
        elem(dist_vec, j) = elem(cube_coord, j) - elem(atom_pos, j);
      else if (elem(atom_pos, j) > (elem(cube_coord, j) + little_cube_size))
        elem(dist_vec, j) =
            elem(atom_pos, j) - (elem(cube_coord, j) + little_cube_size);
      else
        elem(dist_vec, j) = 0;
    }

    atom_cube_dists(i) = length2(dist_vec);
  }
  assign_functions(atom_cube_dists, min_exps, min_coeff);
  compute_nucleii_maps();
}

/*****************************
 * Sphere
 *****************************/
template <class scalar_type>
void PointGroup<scalar_type>::assign_functions_as_sphere(
    uint atom, double radius, const std::vector<double>& min_exps,
    const std::vector<double>& min_coeff) {
  // TODO: esto solo es necesario para los atomos en nucleii, idem arriba
  HostMatrix<double> atom_sphere_dists(fortran_vars.atoms);
  const double3& own_atom_pos = fortran_vars.atom_positions(atom);

  for (uint i = 0; i < fortran_vars.atoms; i++) {
    if (i == atom)
      atom_sphere_dists(i) = 0;
    else {
      const double3& atom_pos = fortran_vars.atom_positions(i);
      double3 dist_vec = (atom_pos - own_atom_pos);
      double dist_to_atom = length(dist_vec);
      double dist = (radius > dist_to_atom ? 0 : dist_to_atom - radius);
      atom_sphere_dists(i) = dist * dist;
    }
  }
  assign_functions(atom_sphere_dists, min_exps, min_coeff);
  compute_nucleii_maps();
}

#if FULL_DOUBLE
template class PointGroup<double>;
#else
template class PointGroup<float>;
#endif
