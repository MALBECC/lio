#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../common.h"
//#include "../cuda_includes.h"
#include "../init.h"
#include "../matrix.h"
#include "../partition.h"
using namespace std;

namespace G2G {
template<class scalar_type>
void PointGroup<scalar_type>::compute_functions(bool forces, bool gga)
{
  /* Load group functions */
  uint group_m = total_functions();

  function_values.resize(group_m, number_of_points);
  if (forces || gga) gradient_values.resize(group_m, number_of_points);
  if (gga) hessian_values.resize(group_m * 2, number_of_points);

  std::vector<Point> _points(points.begin(),points.end());

#pragma omp parallel for
  for(int point = 0; point<_points.size(); point++) {
    vec_type3 point_position = vec_type3(_points[point].position.x, _points[point].position.y, _points[point].position.z);

    for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
      // compute exponential
      uint nuc = func2global_nuc(i);
      vec_type3 v(point_position - vec_type3(fortran_vars.atom_positions(nuc)));
      scalar_type dist = v.length2();

      scalar_type t = 0, tg = 0, th = 0;
      uint global_func = local2global_func[i];
      uint contractions = fortran_vars.contractions(global_func);
      for (uint contraction = 0; contraction < contractions; contraction++) {
        scalar_type a = fortran_vars.a_values(global_func, contraction);
        scalar_type c = fortran_vars.c_values(global_func, contraction);
        scalar_type t0 = exp(-(a * dist)) * c;
        t += t0;
        if (forces || gga) tg += t0 * a;
        if (gga) th += t0 * (a * a);
      }

      vec_type3 vxxy, vyzz;
      if (gga) { vxxy = vec_type3(v.x(), v.x(), v.y()); vyzz = vec_type3(v.y(), v.z(), v.z()); }

      // compute s, p, d
      if (i < s_functions) {
        function_values(ii, point) = t;
        if (forces || gga) {
          gradient_values(ii, point) = vec_type3(v * (-2 * tg));
        }
        if (gga) {
          hessian_values(2 * ii + 0, point) = vec_type3((v * v) * 4 * th - 2 * tg); // Fxx, Fxy, Fxz
          hessian_values(2 * ii + 1, point) = vec_type3(vxxy * vyzz * 4 * th); // Fxy, Fxz, Fyz
        }

        ii++;
      }
      else if (i < s_functions + p_functions) {
        function_values(ii + 0, point) = v.x() * t;
        function_values(ii + 1, point) = v.y() * t;
        function_values(ii + 2, point) = v.z() * t;

        if (forces || gga) {
          gradient_values(ii + 0, point) = vec_type3(vec_type3(t, 0, 0) - v * 2 * tg * v.x());
          gradient_values(ii + 1, point) = vec_type3(vec_type3(0, t, 0) - v * 2 * tg * v.y());
          gradient_values(ii + 2, point) = vec_type3(vec_type3(0, 0, t) - v * 2 * tg * v.z());
        }
        if (gga) {
          hessian_values(2 * (ii + 0) + 0, point) = vec_type3((v * v) *       4 * th * v.x() - vec_type3(6, 2, 2)     * tg * v.x());
          hessian_values(2 * (ii + 0) + 1, point) = vec_type3((vxxy * vyzz) * 4 * th * v.x() - vec_type3(v.y(), v.z(), 0) * 2 * tg);
          hessian_values(2 * (ii + 1) + 0, point) = vec_type3((v * v) *       4 * th * v.y() - vec_type3(2, 6, 2)     * tg * v.y());
          hessian_values(2 * (ii + 1) + 1, point) = vec_type3((vxxy * vyzz) * 4 * th * v.y() - vec_type3(v.x(), 0, v.z()) * 2 * tg);
          hessian_values(2 * (ii + 2) + 0, point) = vec_type3((v * v)       * 4 * th * v.z() - vec_type3(2, 2, 6)     * tg * v.z());
          hessian_values(2 * (ii + 2) + 1, point) = vec_type3((vxxy * vyzz) * 4 * th * v.z() - vec_type3(0, v.x(), v.y()) * 2 * tg);
        }

        ii += 3;
      }
      else {
        function_values(ii + 0, point) = t * v.x() * v.x() * fortran_vars.normalization_factor;
        function_values(ii + 1, point) = t * v.y() * v.x();
        function_values(ii + 2, point) = t * v.y() * v.y() * fortran_vars.normalization_factor;
        function_values(ii + 3, point) = t * v.z() * v.x();
        function_values(ii + 4, point) = t * v.z() * v.y();
        function_values(ii + 5, point) = t * v.z() * v.z() * fortran_vars.normalization_factor;

        if (forces || gga) {
          gradient_values(ii + 0, point) = vec_type3((vec_type3(2 * v.x(), 0      , 0      ) * t - v * 2 * tg * v.x() * v.x()) * fortran_vars.normalization_factor);
          gradient_values(ii + 1, point) = vec_type3(vec_type3(v.y()     , v.x()    , 0      ) * t - v * 2 * tg * v.y() * v.x());
          gradient_values(ii + 2, point) = vec_type3((vec_type3(0      , 2 * v.y(), 0      ) * t - v * 2 * tg * v.y() * v.y()) * fortran_vars.normalization_factor);
          gradient_values(ii + 3, point) = vec_type3( vec_type3(v.z()    , 0      , v.x()    ) * t - v * 2 * tg * v.z() * v.x());
          gradient_values(ii + 4, point) = vec_type3(vec_type3(0       , v.z()    , v.y()    ) * t - v * 2 * tg * v.z() * v.y());
          gradient_values(ii + 5, point) = vec_type3((vec_type3(0      , 0      , 2 * v.z()) * t - v * 2 * tg * v.z() * v.z()) * fortran_vars.normalization_factor);
        }

        if (gga) {
          hessian_values(2 * (ii + 0) + 0, point) = vec_type3(((v * v)       * 4 * th * (v.x() * v.x()) - vec_type3(10, 2, 2) * tg * (v.x() * v.x())    + vec_type3(2 * t, 0    , 0)) * fortran_vars.normalization_factor);
          hessian_values(2 * (ii + 0) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (v.x() * v.x()) - vec_type3(4,  4, 0) * tg * (vxxy * vyzz)                                 ) * fortran_vars.normalization_factor);
          hessian_values(2 * (ii + 1) + 0, point) = vec_type3(((v * v)       * 4 * th * (v.x() * v.y()) - vec_type3(6,  6, 2) * tg * (v.x() * v.y())                                   ));
          hessian_values(2 * (ii + 1) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (v.x() * v.y()) - vec_type3(2 * (v.x() * v.x() + v.y() * v.y()), 2 * v.y() * v.z(), 2 * v.x() * v.z()) * tg + vec_type3(t     , 0    , 0)));
          hessian_values(2 * (ii + 2) + 0, point) = vec_type3(((v * v)       * 4 * th * (v.y() * v.y()) - vec_type3(2, 10, 2) * tg * (v.y() * v.y()) + vec_type3(0    , 2 * t, 0)) * fortran_vars.normalization_factor);
          hessian_values(2 * (ii + 2) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (v.y() * v.y()) - vec_type3(4,  0, 4) * tg * (vxxy * vyzz)                                ) * fortran_vars.normalization_factor);
          hessian_values(2 * (ii + 3) + 0, point) = vec_type3(((v * v)       * 4 * th * (v.x() * v.z()) - vec_type3(6,  2, 6) * tg * (v.x() * v.z())                                  ));
          hessian_values(2 * (ii + 3) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (v.x() * v.z()) - vec_type3(2 * v.y() * v.z(), 2 * (v.x() * v.x() + v.z() * v.z()), 2 * v.x() * v.y()) * tg + vec_type3(0,      t,     0)));
          hessian_values(2 * (ii + 4) + 0, point) = vec_type3(((v * v)       * 4 * th * (v.y() * v.z()) - vec_type3(2,  6, 6) * tg * (v.y() * v.z())                                ));
          hessian_values(2 * (ii + 4) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (v.y() * v.z()) - vec_type3(2 * v.x() * v.z(), 2 * v.x() * v.y(), 2 * (v.y() * v.y() + v.z() * v.z())) * tg + vec_type3(0,      0,     t)));
          hessian_values(2 * (ii + 5) + 0, point) = vec_type3(((v * v)       * 4 * th * (v.z() * v.z()) - vec_type3(2,  2, 10) * tg * (v.z() * v.z()) + vec_type3(0,      0, 2 * t)) * fortran_vars.normalization_factor);
          hessian_values(2 * (ii + 5) + 1, point) = vec_type3(((vxxy * vyzz) * 4 * th * (v.z() * v.z()) - vec_type3(0,  4, 4) * tg * (vxxy * vyzz)                                 ) * fortran_vars.normalization_factor);
        }
        ii += 6;
      }
    }
  }
}

template class PointGroup<double>;
template class PointGroup<float>;
}
