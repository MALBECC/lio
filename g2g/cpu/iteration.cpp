/* -*- mode: c -*- */
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include "../common.h"
#include "../init.h"
#include "../cuda/cuda_extra.h"
#include "../matrix.h"
#include "../timer.h"
#include "../partition.h"

#include "cpu/pot.h"

using namespace std;
using namespace G2G;

#if CPU_KERNELS
extern "C" void g2g_solve_groups_(const uint& computation_type, double* fort_energy_ptr, double* fort_forces_ptr)
{
  Timer timer_total;
  timer_total.start();

 	cout << "<================ iteracion [";
	switch(computation_type) {
		case COMPUTE_ENERGY_ONLY: cout << "energia"; break;
		case COMPUTE_RMM: cout << "rmm"; break;
		case COMPUTE_FORCE_ONLY: cout << "fuerzas"; break;
		case COMPUTE_ENERGY_FORCE: cout << "energia+fuerzas"; break;
	}
	cout << "] ==========>" << endl;

  bool compute_energy = (computation_type == COMPUTE_ENERGY_ONLY || computation_type == COMPUTE_ENERGY_FORCE);
  bool compute_forces = (computation_type == COMPUTE_FORCE_ONLY || computation_type == COMPUTE_ENERGY_FORCE);

  double total_energy = 0;

  HostMatrixFloat rmm_output;
  HostMatrixFloat w(fortran_vars.nco);
  HostMatrixFloat rmm_input;

	for (list<PointGroup>::const_iterator group_it = final_partition.begin(); group_it != final_partition.end(); ++group_it)
  {
		const PointGroup& group = *group_it;
    uint group_m = group.s_functions + group.p_functions * 3 + group.d_functions * 6;

    if (computation_type == COMPUTE_RMM) {
      rmm_output.resize(group_m, group_m); rmm_output.fill(0);
    }

    /******** density *******/
    // each group point
    uint point = 0;
    for (list<Point>::const_iterator point_it = group.points.begin(); point_it != group.points.end(); ++point_it, ++point)
    {
      w.fill(0);
      float partial_density = 0;

      // each group function
      uint ii = 0;
      for (uint i = 0; i < group.functions.size(); i++)
      {
        uint inc;
        if (i < group.s_functions) inc = 1;
        else if (i < group.s_functions + group.p_functions) inc = 3;
        else inc = 6;

        for (uint j = 0; j < inc; j++, ii++) {
          float f = group.function_values.get(point, ii);

          for (uint k = 0; k < fortran_vars.nco; k++) {
            float r = fortran_vars.rmm_input.get(group.functions[i] + j, k); // TODO: hacer una copia para mejor acceso
          }
        }
      }

      for (uint k = 0; k < fortran_vars.nco; k++) { partial_density += w.get(k) * w.get(k); }
      partial_density *= 2;

      /******** energy *******/
      float exc, corr, y2a;
      if (compute_energy) {
        cpu_pot<true, false>(partial_density, exc, corr, y2a);
        total_energy += (partial_density * point_it->weight) * (exc + corr);
      }

      /******** RMM *******/
      if (computation_type == COMPUTE_RMM) {
        cpu_pot<false, true>(partial_density, exc, corr, y2a);
        float factor = point_it->weight * y2a;

        uint small_fi = 0;
        for (vector<uint>::const_iterator it_fi = group.functions.begin(); it_fi != group.functions.end(); ++it_fi) {
          uint fi_advance;
          if (*it_fi < fortran_vars.s_funcs) fi_advance = 1;
          else if (*it_fi < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) fi_advance = 3;
          else fi_advance = 6;

          for (uint i = 0; i < fi_advance; i++) {

            uint small_fj = 0;
            for (vector<uint>::const_iterator it_fj = group.functions.begin(); it_fj != group.functions.end(); ++it_fj) {
              uint fj_advance;
              if (*it_fj < fortran_vars.s_funcs) fj_advance = 1;
              else if (*it_fj < fortran_vars.s_funcs + fortran_vars.p_funcs * 3) fj_advance = 3;
              else fj_advance = 6;

              for (uint j = 0; j < fj_advance; j++) {
                uint fi = *it_fi + i; uint fj = *it_fj + j;
                if (fi > fj) continue;
                uint big_index = (fi * fortran_vars.m - (fi * (fi - 1)) / 2) + (fj - fi);
                //cout << small_fi << " " << small_fj << " " << small_fj + small_fi << endl;
                float Fi = group.function_values.get(point, small_fi);
                float Fj = group.function_values.get(point, small_fj + small_fi);
                fortran_vars.rmm_output.get(big_index) += Fi * Fj * factor;
                small_fj++;
              }
            }
            small_fi++;
          }
        }
      }
    }
  }

  /***** send results to fortran ****/
  if (compute_energy) { *fort_energy_ptr = total_energy; }

  timer_total.stop();
  cout << "TIME: " << timer_total << endl;
}
#endif
