#include <iostream>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include "../common.h"
#include "../cuda/cuda_extra.h"
#include "../init.h"
#include "../matrix.h"
#include "../partition.h"
using namespace std;
using namespace G2G;

template <bool forces, bool gga> void g2g_compute_functions(void)
{
	for (list<PointGroup>::iterator it = final_partition.begin(); it != final_partition.end(); ++it) {
		PointGroup& group = *it;

		/* Load group functions */
		uint group_m = group.total_functions();
    group.function_values.resize(group_m, group.number_of_points);
    if (forces || gga) group.gradient_values.resize(group_m, group.number_of_points);
    if (gga) group.hessian_values.resize(group_m * 2, group.number_of_points);

    uint point = 0;
    for (list<Point>::const_iterator point_it = group.points.begin(); point_it != group.points.end(); ++point_it, ++point) {
      float3 point_position = make_float3(point_it->position.x, point_it->position.y, point_it->position.z);

      for (uint i = 0, ii = 0; i < group.functions.size(); i++) {
        // compute exponential
        uint nuc = fortran_vars.nucleii(group.functions[i]) - 1;
        float3 v = point_position - to_float3(fortran_vars.atom_positions(nuc));
        float dist = length2(v);

        float t = 0, tg = 0, th = 0;
        uint contractions = fortran_vars.contractions(group.functions[i]);
        for (uint contraction = 0; contraction < contractions; contraction++) {
          float a = fortran_vars.a_values(group.functions[i], contraction);
          float c = fortran_vars.c_values(group.functions[i], contraction);
          float t0 = expf(-(a * dist)) * c;
          t += t0;
          if (forces || gga) tg += t0 * a;
          if (gga) th += t0 * (a * a);
        }

				float3 vxxy, vyzz;
				if (gga) { vxxy = make_float3(v.x, v.x, v.y); vyzz = make_float3(v.y, v.z, v.z); }

        // compute s, p, d
        if (i < group.s_functions) {
          group.function_values(ii, point) = t;
          if (forces || gga) {
            group.gradient_values(ii, point) = cfloat3(v * (-2 * tg));
          }
          if (gga) {
            group.hessian_values(2 * ii + 0, point) = cfloat3((v * v) * 4 * th - 2 * tg); // Fxx, Fxy, Fxz
            group.hessian_values(2 * ii + 1, point) = cfloat3(vxxy * vyzz * 4 * th); // Fxy, Fxz, Fyz
          }
          
          ii++;
        }
        else if (i < group.s_functions + group.p_functions) {
          group.function_values(ii + 0, point) = v.x * t;
          group.function_values(ii + 1, point) = v.y * t;
          group.function_values(ii + 2, point) = v.z * t;

          if (forces || gga) {
            group.gradient_values(ii + 0, point) = cfloat3(make_float3(t, 0, 0) - v * 2 * tg * v.x);
            group.gradient_values(ii + 1, point) = cfloat3(make_float3(0, t, 0) - v * 2 * tg * v.y);
            group.gradient_values(ii + 2, point) = cfloat3(make_float3(0, 0, t) - v * 2 * tg * v.z);
          }
          if (gga) {
            group.hessian_values(2 * (ii + 0) + 0, point) = cfloat3(4 * th * v.x * (v * v) - tg * v.x * make_float3(6, 2, 2));
            group.hessian_values(2 * (ii + 0) + 1, point) = cfloat3(4 * th * v.x * (vxxy * vyzz) - 2 * tg * make_float3(v.y, v.z, 0));
            group.hessian_values(2 * (ii + 1) + 0, point) = cfloat3(4 * th * v.y * (v * v) - tg * v.y * make_float3(2, 6, 2));
            group.hessian_values(2 * (ii + 1) + 1, point) = cfloat3(4 * th * v.y * (vxxy * vyzz) - 2 * tg * make_float3(v.x, 0, v.z));
            group.hessian_values(2 * (ii + 2) + 0, point) = cfloat3(4 * th * v.z * (v * v) - tg * v.z * make_float3(2, 2, 6));
            group.hessian_values(2 * (ii + 2) + 1, point) = cfloat3(4 * th * v.z * (vxxy * vyzz) - 2 * tg * make_float3(0, v.x, v.y));
          }

          ii += 3;
        }
        else {
          group.function_values(ii + 0, point) = t * v.x * v.x * fortran_vars.normalization_factor;
          group.function_values(ii + 1, point) = t * v.y * v.x;
          group.function_values(ii + 2, point) = t * v.y * v.y * fortran_vars.normalization_factor;
          group.function_values(ii + 3, point) = t * v.z * v.x;
          group.function_values(ii + 4, point) = t * v.z * v.y;
          group.function_values(ii + 5, point) = t * v.z * v.z * fortran_vars.normalization_factor;

          if (gga) {
						group.hessian_values(2 * (ii + 0) + 0, point) = cfloat3((4 * th * (v.x * v.x) * (v * v)       - tg * (v.x * v.x)   * make_float3(10, 2, 2) + make_float3(2 * t, 0    , 0)) * fortran_vars.normalization_factor);
            group.hessian_values(2 * (ii + 0) + 1, point) = cfloat3((4 * th * (v.x * v.x) * (vxxy * vyzz) - tg * (vxxy * vyzz) * make_float3(4,  4, 0)                               ) * fortran_vars.normalization_factor);
            group.hessian_values(2 * (ii + 1) + 0, point) = cfloat3((4 * th * (v.x * v.y) * (v * v)       - tg * (v.x * v.y)   * make_float3(6,  6, 2)                               ));
            group.hessian_values(2 * (ii + 1) + 1, point) = cfloat3((4 * th * (v.x * v.y) * (vxxy * vyzz) - tg * make_float3(2 * (v.x * v.x + v.y * v.y), 2 * v.y * v.z, 2 * v.x * v.z) + make_float3(t     , 0    , 0)));
            group.hessian_values(2 * (ii + 2) + 0, point) = cfloat3((4 * th * (v.y * v.y) * (v * v)       - tg * (v.y * v.y)   * make_float3(2, 10, 2) + make_float3(0    , 2 * t, 0)) * fortran_vars.normalization_factor);
            group.hessian_values(2 * (ii + 2) + 1, point) = cfloat3((4 * th * (v.y * v.y) * (vxxy * vyzz) - tg * (vxxy * vyzz) * make_float3(4,  0, 4)                                ) * fortran_vars.normalization_factor);
						group.hessian_values(2 * (ii + 3) + 0, point) = cfloat3((4 * th * (v.x * v.z) * (v * v)       - tg * (v.x * v.z)   * make_float3(6,  2, 6)                                ));
            group.hessian_values(2 * (ii + 3) + 1, point) = cfloat3((4 * th * (v.x * v.z) * (vxxy * vyzz) - tg * make_float3(2 * v.y * v.z, 2 * (v.x * v.x + v.z * v.z), 2 * v.x * v.y) + make_float3(0,      t,     0)));
            group.hessian_values(2 * (ii + 4) + 0, point) = cfloat3((4 * th * (v.y * v.z) * (v * v)       - tg * (v.y * v.z)   * make_float3(2,  6, 6)                                ));
            group.hessian_values(2 * (ii + 4) + 1, point) = cfloat3((4 * th * (v.y * v.z) * (vxxy * vyzz) - tg * make_float3(2 * v.x * v.z, 2 * v.x * v.y, 2 * (v.y * v.y + v.z * v.z)) + make_float3(0,      0,     t)));
            group.hessian_values(2 * (ii + 5) + 0, point) = cfloat3((4 * th * (v.z * v.z) * (v * v)       - tg * (v.z * v.z)   * make_float3(2,  2, 10) + make_float3(0,      0, 2 * t)) * fortran_vars.normalization_factor);
            group.hessian_values(2 * (ii + 5) + 1, point) = cfloat3((4 * th * (v.z * v.z) * (vxxy * vyzz) - tg * (vxxy * vyzz) * make_float3(0,  4, 4)                                 ) * fortran_vars.normalization_factor);
          }
          if (forces || gga) {
            group.gradient_values(ii + 0, point) = cfloat3((make_float3(2 * v.x, 0      , 0      ) * t - 2 * tg * v * v.x * v.x) * fortran_vars.normalization_factor);
            group.gradient_values(ii + 1, point) =  cfloat3(make_float3(v.y    , v.x    , 0      ) * t - 2 * tg * v * v.y * v.x);
            group.gradient_values(ii + 2, point) = cfloat3((make_float3(0      , 2 * v.y, 0      ) * t - 2 * tg * v * v.y * v.y) * fortran_vars.normalization_factor);
            group.gradient_values(ii + 3, point) = cfloat3( make_float3(v.z    , 0      , v.x    ) * t - 2 * tg * v * v.z * v.x);
            group.gradient_values(ii + 4, point) =  cfloat3(make_float3(0      , v.z    , v.y    ) * t - 2 * tg * v * v.z * v.y);
            group.gradient_values(ii + 5, point) = cfloat3((make_float3(0      , 0      , 2 * v.z) * t - 2 * tg * v * v.z * v.z) * fortran_vars.normalization_factor);
          }
          ii += 6;
        }
      }
    }
	}
}

template void g2g_compute_functions<true, false>(void);
template void g2g_compute_functions<true, true>(void);
template void g2g_compute_functions<false, false>(void);
template void g2g_compute_functions<false, true>(void);
