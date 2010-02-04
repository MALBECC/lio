#include <iostream>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include "../common.h"
#include "../init.h"
#include "../matrix.h"
#include "../partition.h"
using namespace std;
using namespace G2G;

#if CPU_KERNELS
template <bool forces> void g2g_compute_functions(void)
{
	for (list<PointGroup>::iterator it = final_partition.begin(); it != final_partition.end(); ++it) {
		PointGroup& group = *it;

		/* Load group functions */
		uint group_m = group.s_functions + group.p_functions * 3 + group.d_functions * 6;
    uint group_spd = group.s_functions + group.p_functions + group.d_functions;

    group.function_values.resize(group_m, group.number_of_points);
    if (forces) group.gradient_values.resize(group_m, group.number_of_points);

    uint point = 0;
    for (list<Point>::const_iterator point_it = group.points.begin(); point_it != group.points.end(); ++point_it, ++point) {
      float3 point_position = make_float3(point_it->position.x, point_it->position.y, point_it->position.z);

      uint ii = 0, i = 0;
      for (uint i = 0; i < group.functions.size(); i++) {
        // compute exponential
        uint nuc = fortran_vars.nucleii.get(group.functions[i]) - 1;
        float3 v = point_position - fortran_vars.atom_positions.get(nuc).to_float3();
        float dist = length2(v);

        float t = 0, tg = 0;
        uint contractions = fortran_vars.contractions.get(group.functions[i]);
        for (uint contraction = 0; contraction < contractions; contraction++) {
          float a = fortran_vars.a_values.get(group.functions[i], contraction);
          float c = fortran_vars.c_values.get(group.functions[i], contraction);
          float t0 = exp(-(a * dist)) * c;
          t += t0;
          if (forces) tg += t0 * a;
        }

        float3 c;
        if (forces) c = v * (2.0f * tg);

        // compute s, p, d
        if (i < group.s_functions) {
          group.function_values.get(ii, point) = t;
          if (forces) { group.gradient_values.get(ii, point) = c; }
          ii++;
        }
        else if (i < group.s_functions + group.p_functions) {
          group.function_values.get(ii + 0, point) = v.x * t;
          group.function_values.get(ii + 1, point) = v.y * t;
          group.function_values.get(ii + 2, point) = v.z * t;

          if (forces) {
            group.gradient_values.get(ii + 0, point) = c * v.x - make_float3(t, 0, 0);
            group.gradient_values.get(ii + 1, point) = c * v.y - make_float3(0, t, 0);
            group.gradient_values.get(ii + 2, point) = c * v.z - make_float3(0, 0, t);
          }
          ii += 3;
        }
        else {
          group.function_values.get(ii + 0, point) = t * v.x * v.x * fortran_vars.normalization_factor;
          group.function_values.get(ii + 1, point) = t * v.y * v.x;
          group.function_values.get(ii + 2, point) = t * v.y * v.y * fortran_vars.normalization_factor;
          group.function_values.get(ii + 3, point) = t * v.z * v.x;
          group.function_values.get(ii + 4, point) = t * v.z * v.y;
          group.function_values.get(ii + 5, point) = t * v.z * v.z * fortran_vars.normalization_factor;

          if (forces) {
            group.gradient_values.get(ii + 0, point) = (c * v.x * v.x - make_float3(2 * t * v.x, 0          , 0          )) * fortran_vars.normalization_factor;
            group.gradient_values.get(ii + 1, point) =  c * v.y * v.x - make_float3(t * v.y    , t * v.x    , 0          );
            group.gradient_values.get(ii + 2, point) = (c * v.y * v.y - make_float3(0          , 2 * t * v.y, 0          )) * fortran_vars.normalization_factor;
            group.gradient_values.get(ii + 3, point) =  c * v.z * v.x - make_float3(t * v.z    , 0          , t * v.x    );
            group.gradient_values.get(ii + 4, point) =  c * v.z * v.y - make_float3(0          , t * v.z    , t * v.y    );
            group.gradient_values.get(ii + 5, point) = (c * v.z * v.z - make_float3(0          , 0          , 2 * t * v.z)) * fortran_vars.normalization_factor;
          }
          ii += 6;
        }
      }
    }
	}
}

template void g2g_compute_functions<true>(void);
template void g2g_compute_functions<false>(void);

#endif
