#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "timer.h"
using namespace std;
using namespace G2G;

void g2g_compute_group_weights(PointGroup& group);

/********************
 * PointGroup
 ********************/
void PointGroup::add_point(const Point& p) {
  points.push_back(p);
  number_of_points++;
}

void PointGroup::compute_weights(void) {
  g2g_compute_group_weights(*this);
}

void PointGroup::get_rmm_input(HostMatrixFloat& rmm_input) const {
  rmm_input.zero();
  for (uint i = 0, ii = 0; i < functions.size(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = functions[i] + k;
      for (uint j = 0, jj = 0; j < functions.size(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = functions[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          rmm_input(ii, jj) = fortran_vars.rmm_input_ndens1.data[big_index];
          rmm_input(jj, ii) = rmm_input(ii, jj);
        }
      }
    }
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
Cube::Cube(void) : PointGroup() { }


