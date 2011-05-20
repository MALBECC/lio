#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include <algorithm>
#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "timer.h"
using namespace std;
using namespace G2G;

Partition G2G::partition;

/********************
 * Partition
 ********************/

template <bool forces, bool gga> void Partition::compute_functions(void)
{
	Timer t1;
	t1.start_and_sync();
	for (list<PointGroup*>::iterator it = group_list.begin(); it != group_list.end(); ++it)
    (*it)->compute_functions<forces,gga>();
	t1.stop_and_sync();
	cout << "TIMER: funcs: " << t1 << endl;
}

template void Partition::compute_functions<true, true>(void);
template void Partition::compute_functions<true, false>(void);
template void Partition::compute_functions<false, true>(void);
template void Partition::compute_functions<false, false>(void);

/********************
 * PointGroup
 ********************/
void PointGroup::add_point(const Point& p) {
  points.push_back(p);
  number_of_points++;
}

void PointGroup::get_rmm_input(HostMatrix<real>& rmm_input) const {
  rmm_input.zero();
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          rmm_input(ii, jj) = fortran_vars.rmm_input_ndens1.data[big_index];
          rmm_input(jj, ii) = rmm_input(ii, jj);
        }
      }
    }
  }
}

void PointGroup::add_rmm_output(const HostMatrix<real>& rmm_output) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          fortran_vars.rmm_output(big_index) += rmm_output(ii, jj);
        }
      }
    }
  }
}

void PointGroup::compute_nucleii_maps(void)
{
  if (total_functions_simple() != 0) {
    func2global_nuc.resize(total_functions_simple());
    for (uint i = 0; i < total_functions_simple(); i++) {
      func2global_nuc(i) = fortran_vars.nucleii(local2global_func[i]) - 1;
    }

    func2local_nuc.resize(total_functions());
    uint ii = 0;
    for (uint i = 0; i < total_functions_simple(); i++) {
      uint global_atom = func2global_nuc(i);
      uint local_atom = std::distance(local2global_nuc.begin(), std::find(local2global_nuc.begin(), local2global_nuc.end(), global_atom));
      uint inc = small_function_type(i);
      for (uint k = 0; k < inc; k++, ii++) func2local_nuc(ii) = local_atom;
    }
  }
}

#define EXP_PREFACTOR 1.01057089636005 // (2 * pow(4, 1/3.0)) / M_PI

bool PointGroup::is_significative(FunctionType type, double exponent, double coeff, double d2) {
  switch(type) {
    case FUNCTION_S:
      return (exponent * d2 < max_function_exponent-log(pow((2.*exponent/M_PI),3))/4);
    break;
    default:
    {
      double x = 1;
      double delta;
      double e = 0.1;
      double factor = pow((2.0*exponent/M_PI),3);
      factor = sqrt(factor*4.0*exponent) ;
      double norm = (type == FUNCTION_P ? sqrt(factor) : abs(factor)) ;
      do {
	double div = (type == FUNCTION_P ? log(x) : 2 * log(x));
	double x1 = sqrt((max_function_exponent - log(norm) + div) / exponent);
	delta = abs(x-x1);
	x = x1;
      } while (delta > e);
      return (sqrt(d2) < x);
    }
    break;
  }
}

/**********************
 * Sphere
 **********************/
Sphere::Sphere(void) : PointGroup(), atom(0), radius(0) { }
Sphere::Sphere(uint _atom, double _radius) : PointGroup(), atom(_atom), radius(_radius) { }

/**********************
 * Cube
 **********************/
Cube::Cube(void) : PointGroup() {  }


