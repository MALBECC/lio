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
void g2g_compute_functions(void)
{
	for (list<PointGroup>::iterator it = final_partition.begin(); it != final_partition.end(); ++it) {
		PointGroup& group = *it;

		/* Load group functions */
		uint group_m = group.s_functions + group.p_functions * 3 + group.d_functions * 6;
    uint group_spd = group.s_functions + group.p_functions + group.d_functions;

    group.function_values.resize(group.number_of_points, group_m);

    uint point = 0;
    for (list<Point>::const_iterator point_it = group.points.begin(); point_it != group.points.end(); ++point_it, ++point) {
      float3 point_position = make_float3(point_it->position.x, point_it->position.y, point_it->position.z);

      uint ii = 0, i = 0;
      for (uint i = 0; i < group.functions.size(); i++) {
        // compute exponential
        uint nuc = fortran_vars.nucleii.get(group.functions[i]) - 1;
        float3 v = point_position - fortran_vars.atom_positions.get(nuc).to_float3();
        float dist = length2(v);

        float t = 0;
        uint contractions = fortran_vars.contractions.get(group.functions[i]);
        for (uint contraction = 0; contraction < contractions; contraction++) {
          float a = fortran_vars.a_values.get(group.functions[i], contraction);
          float c = fortran_vars.c_values.get(group.functions[i], contraction);
          t += expf(-(a * dist)) * c;
        }

        // compute s, p, d
        if (i < group.s_functions) { group.function_values.get(point, ii) = t; ii++; }
        else if (i < group.s_functions + group.p_functions) {
          group.function_values.get(point, ii) = v.x * t; ii++;
          group.function_values.get(point, ii) = v.y * t; ii++;
          group.function_values.get(point, ii) = v.z * t; ii++;
        }
        else {
          group.function_values.get(point, ii) = t * v.x * v.x * fortran_vars.normalization_factor; ii++;
          group.function_values.get(point, ii) = t * v.y * v.x;                                     ii++;
          group.function_values.get(point, ii) = t * v.y * v.y * fortran_vars.normalization_factor; ii++;
          group.function_values.get(point, ii) = t * v.z * v.x;                                     ii++;
          group.function_values.get(point, ii) = t * v.z * v.y;                                     ii++;
          group.function_values.get(point, ii) = t * v.z * v.z * fortran_vars.normalization_factor; ii++;
        }
      }
    }
	}
}
#endif
