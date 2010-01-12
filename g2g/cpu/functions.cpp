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

    group.function_values.resize(group.number_of_points, group_m);
    if (forces) group.gradient_values.resize(group.number_of_points, group_m);

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
          t += expf(-(a * dist)) * c;
          if (forces) tg += t * a;
        }

        // compute s, p, d
        if (i < group.s_functions) {
          group.function_values.get(point, ii) = t;
          if (forces) group.gradient_values.get(point, ii) = v * tg;
          ii++;
        }
        else if (i < group.s_functions + group.p_functions) {
          group.function_values.get(point, ii + 0) = v.x * t;
          group.function_values.get(point, ii + 1) = v.y * t;
          group.function_values.get(point, ii + 2) = v.z * t;

          if (forces) {
            float3 c = v * (2.0f * tg);
            group.gradient_values.get(point, ii + 0) = v * c.x - make_float3(t, 0, 0);
            group.gradient_values.get(point, ii + 1) = v * c.y - make_float3(0, t, 0);
            group.gradient_values.get(point, ii + 2) = v * c.z - make_float3(0, 0, t);
          }
          ii += 3;
        }
        else {
          group.function_values.get(point, ii + 0) = t * v.x * v.x * fortran_vars.normalization_factor;
          group.function_values.get(point, ii + 1) = t * v.y * v.x;
          group.function_values.get(point, ii + 2) = t * v.y * v.y * fortran_vars.normalization_factor;
          group.function_values.get(point, ii + 3) = t * v.z * v.x;
          group.function_values.get(point, ii + 4) = t * v.z * v.y;
          group.function_values.get(point, ii + 5) = t * v.z * v.z * fortran_vars.normalization_factor;

          if (forces) {
            group.gradient_values.get(point, ii + 0) = v * 2.0f * tg * v.x * v.x * fortran_vars.normalization_factor  - make_float3(2 * t * v.x * fortran_vars.normalization_factor, 0, 0);
            group.gradient_values.get(point, ii + 1) = v * 2.0f * tg * v.y * v.x                                      - make_float3(t * v.y, t * v.x, 0);
            group.gradient_values.get(point, ii + 2) = v * 2.0f * tg * v.y * v.y * fortran_vars.normalization_factor  - make_float3(0, 2 * t * v.y * fortran_vars.normalization_factor , 0);
            group.gradient_values.get(point, ii + 3) = v * 2.0f * tg * v.z * v.x                                      - make_float3(t * v.z, 0, t * v.x);
            group.gradient_values.get(point, ii + 4) = v * 2.0f * tg * v.z * v.y                                      - make_float3(0, t * v.z, t * v.y);
            group.gradient_values.get(point, ii + 5) = v * 2.0f * tg * v.z * v.z * fortran_vars.normalization_factor  - make_float3(0, 0, 2 * t * v.z * fortran_vars.normalization_factor);
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
