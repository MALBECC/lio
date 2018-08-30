#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../common.h"
#include "../cuda_includes.h"
#include "../init.h"
#include "../matrix.h"
#include "../partition.h"
using namespace std;

namespace G2G {

template <class scalar_type>
void PointGroupCPU<scalar_type>::compute_functions(bool forces, bool gga) {
#if !CPU_RECOMPUTE && GPU_KERNELS
  if (this->inGlobal) return;
  this->inGlobal = true;
  forces = gga = true;  // Vamos a cachear asi que guardemos todo y listo
#endif
  /* Load group functions */
  uint group_m = this->total_functions();
  int numpoints = this->number_of_points;

  function_values.resize(group_m, numpoints);
  if (forces || gga) {
    gX.resize(group_m, numpoints);
    gX.zero();
    gY.resize(group_m, numpoints);
    gY.zero();
    gZ.resize(group_m, numpoints);
    gZ.zero();
  }
  if (gga) {
    hPX.resize(group_m, numpoints);
    hPX.zero();
    hPY.resize(group_m, numpoints);
    hPY.zero();
    hPZ.resize(group_m, numpoints);
    hPZ.zero();
    hIX.resize(group_m, numpoints);
    hIX.zero();
    hIY.resize(group_m, numpoints);
    hIY.zero();
    hIZ.resize(group_m, numpoints);
    hIZ.zero();
  }

#pragma omp parallel for schedule(static)
  for (int point = 0; point < this->points.size(); point++) {
    vec_type3 point_position = vec_type3(this->points[point].position.x,
                                         this->points[point].position.y,
                                         this->points[point].position.z);

    for (uint i = 0, ii = 0; i < this->total_functions_simple(); i++) {
      // compute exponential
      uint nuc = this->func2global_nuc(i);
      const vec_type3 v(point_position -
                        vec_type3(fortran_vars.atom_positions(nuc)));
      scalar_type dist = v.length2();

      scalar_type t = 0, tg = 0, th = 0;
      uint global_func = this->local2global_func[i];
      uint contractions = fortran_vars.contractions(global_func);
      for (uint contraction = 0; contraction < contractions; contraction++) {
        scalar_type a = fortran_vars.a_values(global_func, contraction);
        scalar_type c = fortran_vars.c_values(global_func, contraction);
        scalar_type t0 = exp(-(a * dist)) * c;
        t += t0;
        if (forces || gga) tg += t0 * a;
        if (gga) th += t0 * (a * a);
      }

      const vec_type3 vxxy(v.x, v.x, v.y);
      const vec_type3 vyzz(v.y, v.z, v.z);

      // compute s, p, d
      if (i < this->s_functions) {
        function_values(ii, point) = t;
        if (forces || gga) {
          // gradient_values(ii, point) = vec_type3(v * (-2 * tg));
          gX(ii, point) = v.x * -2 * tg;
          gY(ii, point) = v.y * -2 * tg;
          gZ(ii, point) = v.z * -2 * tg;
        }
        if (gga) {
          // hessian_values(2 * ii + 0, point) = vec_type3((v * v) * 4 * th - 2
          // * tg); // Fxx, Fxy, Fxz
          hPX(ii, point) = v.x * v.x * 4 * th - 2 * tg;
          hPY(ii, point) = v.y * v.y * 4 * th - 2 * tg;
          hPZ(ii, point) = v.z * v.z * 4 * th - 2 * tg;
          // hessian_values(2 * ii + 1, point) = vec_type3(vxxy * vyzz * 4 *
          // th); // Fxy, Fxz, Fyz
          hIX(ii, point) = vxxy.x * vyzz.x * 4 * th;
          hIY(ii, point) = vxxy.y * vyzz.y * 4 * th;
          hIZ(ii, point) = vxxy.z * vyzz.z * 4 * th;
        }

        ii++;
      } else if (i < this->s_functions + this->p_functions) {
        function_values(ii + 0, point) = v.x * t;
        function_values(ii + 1, point) = v.y * t;
        function_values(ii + 2, point) = v.z * t;

        if (forces || gga) {
          // gradient_values(ii + 0, point) = vec_type3(vec_type3(t, 0, 0) - v *
          // 2 * tg * v.x);
          gX(ii, point) = t - v.x * 2 * tg * v.x;
          gY(ii, point) = 0 - v.y * 2 * tg * v.x;
          gZ(ii, point) = 0 - v.z * 2 * tg * v.x;
          // gradient_values(ii + 1, point) = vec_type3(vec_type3(0, t, 0) - v *
          // 2 * tg * v.y);
          gX(ii + 1, point) = 0 - v.x * 2 * tg * v.y;
          gY(ii + 1, point) = t - v.y * 2 * tg * v.y;
          gZ(ii + 1, point) = 0 - v.z * 2 * tg * v.y;
          // gradient_values(ii + 2, point) = vec_type3(vec_type3(0, 0, t) - v *
          // 2 * tg * v.z);
          gX(ii + 2, point) = 0 - v.x * 2 * tg * v.z;
          gY(ii + 2, point) = 0 - v.y * 2 * tg * v.z;
          gZ(ii + 2, point) = t - v.z * 2 * tg * v.z;
        }
        if (gga) {
          // hessian_values(2 * (ii + 0) + 0, point) = vec_type3((v * v) *
          // 4 * th * v.x - vec_type3(6, 2, 2)     * tg * v.x);
          hPX(ii, point) = v.x * v.x * 4 * th * v.x - 6 * tg * v.x;
          hPY(ii, point) = v.y * v.y * 4 * th * v.x - 2 * tg * v.x;
          hPZ(ii, point) = v.z * v.z * 4 * th * v.x - 2 * tg * v.x;
          // hessian_values(2 * (ii + 0) + 1, point) = vec_type3((vxxy * vyzz) *
          // 4 * th * v.x - vec_type3(v.y, v.z, 0) * 2 * tg);
          hIX(ii, point) = vxxy.x * vyzz.x * 4 * th * v.x - v.y * 2 * tg;
          hIY(ii, point) = vxxy.y * vyzz.y * 4 * th * v.x - v.z * 2 * tg;
          hIZ(ii, point) = vxxy.z * vyzz.z * 4 * th * v.x - 0 * 2 * tg;
          // hessian_values(2 * (ii + 1) + 0, point) = vec_type3((v * v) *
          // 4 * th * v.y - vec_type3(2, 6, 2)     * tg * v.y);
          hPX(ii + 1, point) = v.x * v.x * 4 * th * v.y - 2 * tg * v.y;
          hPY(ii + 1, point) = v.y * v.y * 4 * th * v.y - 6 * tg * v.y;
          hPZ(ii + 1, point) = v.z * v.z * 4 * th * v.y - 2 * tg * v.y;
          // hessian_values(2 * (ii + 1) + 1, point) = vec_type3((vxxy * vyzz) *
          // 4 * th * v.y - vec_type3(v.x, 0, v.z) * 2 * tg);
          hIX(ii + 1, point) = vxxy.x * vyzz.x * 4 * th * v.y - 2 * tg * v.x;
          hIY(ii + 1, point) = vxxy.y * vyzz.y * 4 * th * v.y - 2 * tg * 0;
          hIZ(ii + 1, point) = vxxy.z * vyzz.z * 4 * th * v.y - 2 * tg * v.z;
          // hessian_values(2 * (ii + 2) + 0, point) = vec_type3((v * v)       *
          // 4 * th * v.z - vec_type3(2, 2, 6)     * tg * v.z);
          hPX(ii + 2, point) = v.x * v.x * 4 * th * v.z - 2 * tg * v.z;
          hPY(ii + 2, point) = v.y * v.y * 4 * th * v.z - 2 * tg * v.z;
          hPZ(ii + 2, point) = v.z * v.z * 4 * th * v.z - 6 * tg * v.z;
          //          hessian_values(2 * (ii + 2) + 1, point) = vec_type3((vxxy
          //          * vyzz) * 4 * th * v.z - vec_type3(0, v.x, v.y) * 2 * tg);
          hIX(ii + 2, point) = vxxy.x * vyzz.x * 4 * th * v.z - 2 * tg * 0;
          hIY(ii + 2, point) = vxxy.y * vyzz.y * 4 * th * v.z - 2 * tg * v.x;
          hIZ(ii + 2, point) = vxxy.z * vyzz.z * 4 * th * v.z - 2 * tg * v.y;
        }

        ii += 3;
      } else {
        function_values(ii + 0, point) =
            t * v.x * v.x * fortran_vars.normalization_factor;
        function_values(ii + 1, point) = t * v.y * v.x;
        function_values(ii + 2, point) =
            t * v.y * v.y * fortran_vars.normalization_factor;
        function_values(ii + 3, point) = t * v.z * v.x;
        function_values(ii + 4, point) = t * v.z * v.y;
        function_values(ii + 5, point) =
            t * v.z * v.z * fortran_vars.normalization_factor;

        if (forces || gga) {
          // gradient_values(ii + 0, point) = vec_type3((vec_type3(2 * v.x, 0
          // , 0      ) * t - v * 2 * tg * v.x * v.x) *
          // fortran_vars.normalization_factor);
          gX(ii, point) = (2 * v.x * t - v.x * 2 * tg * v.x * v.x) *
                          fortran_vars.normalization_factor;
          gY(ii, point) = (0 * t - v.y * 2 * tg * v.x * v.x) *
                          fortran_vars.normalization_factor;
          gZ(ii, point) = (0 * t - v.z * 2 * tg * v.x * v.x) *
                          fortran_vars.normalization_factor;
          // gradient_values(ii + 1, point) = vec_type3(vec_type3(v.y     , v.x
          // , 0      ) * t - v * 2 * tg * v.y * v.x);
          gX(ii + 1, point) = v.y * t - v.x * 2 * tg * v.y * v.x;
          gY(ii + 1, point) = v.x * t - v.y * 2 * tg * v.y * v.x;
          gZ(ii + 1, point) = 0 * t - v.z * 2 * tg * v.y * v.x;
          // gradient_values(ii + 2, point) = vec_type3((vec_type3(0      , 2 *
          // v.y, 0      ) * t - v * 2 * tg * v.y * v.y) *
          // fortran_vars.normalization_factor);
          gX(ii + 2, point) = (0 * t - v.x * 2 * tg * v.y * v.y) *
                              fortran_vars.normalization_factor;
          gY(ii + 2, point) = (2 * v.y * t - v.y * 2 * tg * v.y * v.y) *
                              fortran_vars.normalization_factor;
          gZ(ii + 2, point) = (0 * t - v.z * 2 * tg * v.y * v.y) *
                              fortran_vars.normalization_factor;
          // gradient_values(ii + 3, point) = vec_type3( vec_type3(v.z    , 0
          // , v.x    ) * t - v * 2 * tg * v.z * v.x);
          gX(ii + 3, point) = v.z * t - v.x * 2 * tg * v.z * v.x;
          gY(ii + 3, point) = 0 * t - v.y * 2 * tg * v.z * v.x;
          gZ(ii + 3, point) = v.x * t - v.z * 2 * tg * v.z * v.x;
          // gradient_values(ii + 4, point) = vec_type3(vec_type3(0       , v.z
          // , v.y    ) * t - v * 2 * tg * v.z * v.y);
          gX(ii + 4, point) = 0 * t - v.x * 2 * tg * v.z * v.y;
          gY(ii + 4, point) = v.z * t - v.y * 2 * tg * v.z * v.y;
          gZ(ii + 4, point) = v.y * t - v.z * 2 * tg * v.z * v.y;
          // gradient_values(ii + 5, point) = vec_type3((vec_type3(0      , 0
          // , 2 * v.z) * t - v * 2 * tg * v.z * v.z) *
          // fortran_vars.normalization_factor);
          gX(ii + 5, point) = (0 * t - v.x * 2 * tg * v.z * v.z) *
                              fortran_vars.normalization_factor;
          gY(ii + 5, point) = (0 * t - v.y * 2 * tg * v.z * v.z) *
                              fortran_vars.normalization_factor;
          gZ(ii + 5, point) = (2 * v.z * t - v.z * 2 * tg * v.z * v.z) *
                              fortran_vars.normalization_factor;
        }

        if (gga) {
          //          hessian_values(2 * (ii + 0) + 0, point) = vec_type3(((v *
          //          v)       * 4 * th * (v.x * v.x) - vec_type3(10, 2, 2) * tg
          //          * (v.x * v.x)    + vec_type3(2 * t, 0    , 0)) *
          //          fortran_vars.normalization_factor);
          hPX(ii, point) =
              (v.x * v.x * 4 * th * v.x * v.x - 10 * tg * v.x * v.x + 2 * t) *
              fortran_vars.normalization_factor;
          hPY(ii, point) =
              (v.y * v.y * 4 * th * v.x * v.x - 2 * tg * v.x * v.x + 0) *
              fortran_vars.normalization_factor;
          hPZ(ii, point) =
              (v.z * v.z * 4 * th * v.x * v.x - 2 * tg * v.x * v.x + 0) *
              fortran_vars.normalization_factor;
          //          hessian_values(2 * (ii + 0) + 1, point) = vec_type3(((vxxy
          //          * vyzz) * 4 * th * (v.x * v.x) - vec_type3(4,  4, 0) * tg
          //          * (vxxy * vyzz)                                 ) *
          //          fortran_vars.normalization_factor);
          hIX(ii, point) = (vxxy.x * vyzz.x * 4 * th * v.x * v.x -
                            4 * tg * vxxy.x * vyzz.x) *
                           fortran_vars.normalization_factor;
          hIY(ii, point) = (vxxy.y * vyzz.y * 4 * th * v.x * v.x -
                            4 * tg * vxxy.y * vyzz.y) *
                           fortran_vars.normalization_factor;
          hIZ(ii, point) = (vxxy.z * vyzz.z * 4 * th * v.x * v.x -
                            0 * tg * vxxy.z * vyzz.z) *
                           fortran_vars.normalization_factor;
          //          hessian_values(2 * (ii + 1) + 0, point) = vec_type3(((v *
          //          v)       * 4 * th * (v.x * v.y) - vec_type3(6,  6, 2) * tg
          //          * (v.x * v.y)                                   ));
          hPX(ii + 1, point) =
              v.x * v.x * 4 * th * v.x * v.y - 6 * tg * v.x * v.y;
          hPY(ii + 1, point) =
              v.y * v.y * 4 * th * v.x * v.y - 6 * tg * v.x * v.y;
          hPZ(ii + 1, point) =
              v.z * v.z * 4 * th * v.x * v.y - 2 * tg * v.x * v.y;
          //          hessian_values(2 * (ii + 1) + 1, point) = vec_type3(((vxxy
          //          * vyzz) * 4 * th * (v.x * v.y) - vec_type3(2 * (v.x * v.x
          //          + v.y * v.y), 2 * v.y * v.z, 2 * v.x * v.z) * tg +
          //          vec_type3(t     , 0    , 0)));
          hIX(ii + 1, point) = vxxy.x * vyzz.x * 4 * th * v.x * v.y -
                               2 * (v.x * v.x + v.y * v.y) * tg + t;
          hIY(ii + 1, point) =
              vxxy.y * vyzz.y * 4 * th * v.x * v.y - 2 * v.y * v.z * tg + 0;
          hIZ(ii + 1, point) =
              vxxy.z * vyzz.z * 4 * th * v.x * v.y - 2 * v.x * v.z * tg + 0;
          //          hessian_values(2 * (ii + 2) + 0, point) = vec_type3(((v *
          //          v)       * 4 * th * (v.y * v.y) - vec_type3(2, 10, 2) * tg
          //          * (v.y * v.y) + vec_type3(0    , 2 * t, 0)) *
          //          fortran_vars.normalization_factor);
          hPX(ii + 2, point) =
              (v.x * v.x * 4 * th * v.y * v.y - 2 * tg * v.y * v.y + 0) *
              fortran_vars.normalization_factor;
          hPY(ii + 2, point) =
              (v.y * v.y * 4 * th * v.y * v.y - 10 * tg * v.y * v.y + 2 * t) *
              fortran_vars.normalization_factor;
          hPZ(ii + 2, point) =
              (v.z * v.z * 4 * th * v.y * v.y - 2 * tg * v.y * v.y + 0) *
              fortran_vars.normalization_factor;
          //          hessian_values(2 * (ii + 2) + 1, point) = vec_type3(((vxxy
          //          * vyzz) * 4 * th * (v.y * v.y) - vec_type3(4,  0, 4) * tg
          //          * (vxxy * vyzz)                                ) *
          //          fortran_vars.normalization_factor);
          hIX(ii + 2, point) = (vxxy.x * vyzz.x * 4 * th * (v.y * v.y) -
                                4 * tg * vxxy.x * vyzz.x) *
                               fortran_vars.normalization_factor;
          hIY(ii + 2, point) = (vxxy.y * vyzz.y * 4 * th * (v.y * v.y) -
                                0 * tg * vxxy.y * vyzz.y) *
                               fortran_vars.normalization_factor;
          hIZ(ii + 2, point) = (vxxy.z * vyzz.z * 4 * th * (v.y * v.y) -
                                4 * tg * vxxy.z * vyzz.z) *
                               fortran_vars.normalization_factor;
          //          hessian_values(2 * (ii + 3) + 0, point) = vec_type3(((v *
          //          v)       * 4 * th * (v.x * v.z) - vec_type3(6,  2, 6) * tg
          //          * (v.x * v.z)                                  ));
          hPX(ii + 3, point) =
              v.x * v.x * 4 * th * v.x * v.z - 6 * tg * v.x * v.z;
          hPY(ii + 3, point) =
              v.y * v.y * 4 * th * v.x * v.z - 2 * tg * v.x * v.z;
          hPZ(ii + 3, point) =
              v.z * v.z * 4 * th * v.x * v.z - 6 * tg * v.x * v.z;
          //          hessian_values(2 * (ii + 3) + 1, point) = vec_type3(((vxxy
          //          * vyzz) * 4 * th * (v.x * v.z) - vec_type3(2 * v.y * v.z,
          //          2 * (v.x * v.x + v.z * v.z), 2 * v.x * v.y) * tg +
          //          vec_type3(0,      t,     0)));
          hIX(ii + 3, point) =
              vxxy.x * vyzz.x * 4 * th * v.x * v.z - 2 * v.y * v.z * tg + 0;
          hIY(ii + 3, point) = vxxy.y * vyzz.y * 4 * th * v.x * v.z -
                               2 * (v.x * v.x + v.z * v.z) * tg + t;
          hIZ(ii + 3, point) =
              vxxy.z * vyzz.z * 4 * th * v.x * v.z - 2 * v.x * v.y * tg + 0;
          //           hessian_values(2 * (ii + 4) + 0, point) = vec_type3(((v *
          //           v)       * 4 * th * (v.y * v.z) - vec_type3(2,  6, 6) *
          //           tg * (v.y * v.z)                                ));
          hPX(ii + 4, point) =
              v.x * v.x * 4 * th * v.y * v.z - 2 * tg * v.y * v.z;
          hPY(ii + 4, point) =
              v.y * v.y * 4 * th * v.y * v.z - 6 * tg * v.y * v.z;
          hPZ(ii + 4, point) =
              v.z * v.z * 4 * th * v.y * v.z - 6 * tg * v.y * v.z;
          //          hessian_values(2 * (ii + 4) + 1, point) = vec_type3(((vxxy
          //          * vyzz) * 4 * th * (v.y * v.z) - vec_type3(2 * v.x * v.z,
          //          2 * v.x * v.y, 2 * (v.y * v.y + v.z * v.z)) * tg +
          //          vec_type3(0,      0,     t)));
          hIX(ii + 4, point) =
              vxxy.x * vyzz.x * 4 * th * v.y * v.z - 2 * v.x * v.z * tg + 0;
          hIY(ii + 4, point) =
              vxxy.y * vyzz.y * 4 * th * v.y * v.z - 2 * v.x * v.y * tg + 0;
          hIZ(ii + 4, point) = vxxy.z * vyzz.z * 4 * th * v.y * v.z -
                               2 * (v.y * v.y + v.z * v.z) * tg + t;
          //          hessian_values(2 * (ii + 5) + 0, point) = vec_type3(((v *
          //          v)       * 4 * th * (v.z * v.z) - vec_type3(2,  2, 10) *
          //          tg * (v.z * v.z) + vec_type3(0,      0, 2 * t)) *
          //          fortran_vars.normalization_factor);
          hPX(ii + 5, point) =
              (v.x * v.x * 4 * th * v.z * v.z - 2 * tg * v.z * v.z + 0) *
              fortran_vars.normalization_factor;
          hPY(ii + 5, point) =
              (v.y * v.y * 4 * th * v.z * v.z - 2 * tg * v.z * v.z + 0) *
              fortran_vars.normalization_factor;
          hPZ(ii + 5, point) =
              (v.z * v.z * 4 * th * v.z * v.z - 10 * tg * v.z * v.z + 2 * t) *
              fortran_vars.normalization_factor;
          // hessian_values(2 * (ii + 5) + 1, point) = vec_type3(((vxxy * vyzz)
          // * 4 * th * (v.z * v.z) - vec_type3(0,  4, 4) * tg * (vxxy * vyzz)
          // ) * fortran_vars.normalization_factor);
          hIX(ii + 5, point) = (vxxy.x * vyzz.x * 4 * th * v.z * v.z -
                                0 * tg * vxxy.x * vyzz.x) *
                               fortran_vars.normalization_factor;
          hIY(ii + 5, point) = (vxxy.y * vyzz.y * 4 * th * v.z * v.z -
                                4 * tg * vxxy.y * vyzz.y) *
                               fortran_vars.normalization_factor;
          hIZ(ii + 5, point) = (vxxy.z * vyzz.z * 4 * th * v.z * v.z -
                                4 * tg * vxxy.z * vyzz.z) *
                               fortran_vars.normalization_factor;
        }
        ii += 6;
      }
    }
  }

  function_values.transpose(function_values_transposed);
}

#if FULL_DOUBLE
template class PointGroupCPU<double>;
#else
template class PointGroupCPU<float>;
#endif
}
