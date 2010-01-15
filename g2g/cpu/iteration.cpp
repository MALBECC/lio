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

  Timer t_ciclos;

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

  HostMatrixFloat w(fortran_vars.nco);

  HostMatrixFloat3 density_derivs, forces;
  if (compute_forces) { density_derivs.resize(fortran_vars.atoms, 1); forces.resize(fortran_vars.atoms, 1); }

  Timer t_rmm, t_resto;

  /********** iterate all groups ***********/
	for (list<PointGroup>::const_iterator group_it = final_partition.begin(); group_it != final_partition.end(); ++group_it)
  {
		const PointGroup& group = *group_it;
    uint group_m = group.s_functions + group.p_functions * 3 + group.d_functions * 6;

    /******** each point *******/
    uint point = 0;
    for (list<Point>::const_iterator point_it = group.points.begin(); point_it != group.points.end(); ++point_it, ++point)
    {
      t_resto.start();
      /** density **/
      float partial_density = 0;
      w.fill(0);

      uint ii = 0;
      for (uint i = 0; i < group.functions.size(); i++)
      {
        uint inc;
        if (i < group.s_functions) inc = 1;
        else if (i < group.s_functions + group.p_functions) inc = 3;
        else inc = 6;

        uint big_i = group.functions[i];

        for (uint j = 0; j < inc; j++, ii++) {
          float f = group.function_values.get(point, ii);

          for (uint k = 0; k < fortran_vars.nco; k++) {
            float r = fortran_vars.rmm_input.get(big_i + j, k);
            w.get(k) += f * r;
          }
        }
      }

      for (uint k = 0; k < fortran_vars.nco; k++) { partial_density += w.get(k) * w.get(k); }
      partial_density *= 2;

      /** energy / potential **/
      float exc = 0, corr = 0, y2a = 0;
      if (compute_energy) {
        if (compute_forces) cpu_pot<true, true>(partial_density, exc, corr, y2a);
        else cpu_pot<true, false>(partial_density, exc, corr, y2a);
        total_energy += (partial_density * point_it->weight) * (exc + corr);
      }
      else cpu_pot<false, true>(partial_density, exc, corr, y2a);

      /** density derivatives / forces **/
      if (compute_forces) {
        float factor = point_it->weight * y2a;
        density_derivs.fill(make_float3(0,0,0));

        ii = 0;
        for (uint i = 0; i < group.functions.size(); i++)
        {
          uint inc;
          if (i < group.s_functions) inc = 1;
          else if (i < group.s_functions + group.p_functions) inc = 3;
          else inc = 6;

          uint big_i = group.functions[i];

          for (uint j = 0; j < inc; j++, ii++) {
            float wrdm = 0;
            for (uint k = 0; k < fortran_vars.nco; k++) {
              float r = fortran_vars.rmm_input.get(big_i + j, k);
              wrdm += r * w.get(k);
              uint nuc = fortran_vars.nucleii.get(big_i + j) - 1;
              density_derivs.get(nuc) += group.gradient_values.get(point, ii) * wrdm;
            }
          }
        }

        for (set<uint>::const_iterator it = group.nucleii.begin(); it != group.nucleii.end(); ++it)
          forces.get(*it) += density_derivs.get(*it) * factor;
      }
      t_resto.pause();
      t_rmm.start();

      /******** RMM *******/
      if (computation_type == COMPUTE_RMM) {
        float factor = point_it->weight * y2a;

        ii = 0;
        for (uint i = 0; i < group.functions.size(); i++) {
          uint inc_i;
          if (i < group.s_functions) inc_i = 1;
          else if (i < group.s_functions + group.p_functions) inc_i = 3;
          else inc_i = 6;

          for (uint k = 0; k < inc_i; k++, ii++) {
            float Fi = group.function_values.get(point, ii);
            uint big_i = group.functions[i] + k;

            uint jj = 0;
            for (uint j = 0; j < group.functions.size(); j++) {
              uint inc_j;
              if (j < group.s_functions) inc_j = 1;
              else if (j < group.s_functions + group.p_functions) inc_j = 3;
              else inc_j = 6;

              for (uint l = 0; l < inc_j; l++, jj++) {
                uint big_j = group.functions[j] + l;
                if (big_i > big_j) continue;

                float Fj = group.function_values.get(point, jj);
                uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
                fortran_vars.rmm_output.get(big_index) += Fi * Fj * factor;
              }
            }
          }
        }
      }
      t_rmm.pause();
    }
  }

  /***** send results to fortran ****/
  if (compute_energy) *fort_energy_ptr = total_energy;
  
  if (compute_forces) {
    FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, FORTRAN_MAX_ATOMS);
    for (uint i = 0; i < fortran_vars.atoms; i++)  {
      float3 this_force = forces.get(i);
      fort_forces.get(i,0) += this_force.x;
      fort_forces.get(i,1) += this_force.y;
      fort_forces.get(i,2) += this_force.z;
    }
  }

  timer_total.stop();
  cout << "iteration: " << timer_total << endl;
  cout << "rmm: " << t_rmm << " t_resto: " << t_resto << endl;
}
#endif
