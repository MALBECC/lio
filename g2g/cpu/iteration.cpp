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

/**
 * Nota: tener presente que el get() puede llegar a ser muy costoso
 */


#if CPU_KERNELS
extern "C" void g2g_solve_groups_(const uint& computation_type, double* fort_energy_ptr, double* fort_forces_ptr)
{
  Timer timer_total;
  timer_total.start();

  Timer t_ciclos;

 	cout << "<================ iteracion [";
	switch(computation_type) {
    case COMPUTE_RMM: cout << "rmm"; break;
		case COMPUTE_ENERGY_ONLY: cout << "energia"; break;
		case COMPUTE_ENERGY_FORCE: cout << "energia+fuerzas"; break;
		case COMPUTE_FORCE_ONLY: cout << "fuerzas"; break;		
	}
	cout << "] ==========>" << endl;

  bool compute_energy = (computation_type == COMPUTE_ENERGY_ONLY || computation_type == COMPUTE_ENERGY_FORCE);
  bool compute_forces = (computation_type == COMPUTE_FORCE_ONLY || computation_type == COMPUTE_ENERGY_FORCE);
  bool compute_rmm = (computation_type != COMPUTE_ENERGY_ONLY);

  double total_energy = 0;

  HostMatrixFloat3 density_derivs, forces;
  if (compute_forces) { density_derivs.resize(fortran_vars.atoms, 1); forces.resize(fortran_vars.atoms, 1); forces.zero(); }

  HostMatrixFloat rmm_output;

  Timer t_rmm, t_density, t_resto;

  /********** iterate all groups ***********/
	for (list<PointGroup>::const_iterator group_it = final_partition.begin(); group_it != final_partition.end(); ++group_it)
  {
		const PointGroup& group = *group_it;
    uint group_m = group.s_functions + group.p_functions * 3 + group.d_functions * 6;
    if (compute_rmm) { rmm_output.resize(group_m, group_m); rmm_output.zero(); }

    // prepare rmm_input for this group
    t_density.start();
    HostMatrixFloat rmm_input(group_m, group_m);
    rmm_input.zero();
    uint ii = 0;
    for (uint i = 0; i < group.functions.size(); i++) {
      uint inc_i;
      if (i < group.s_functions) inc_i = 1;
      else if (i < group.s_functions + group.p_functions) inc_i = 3;
      else inc_i = 6;

      for (uint k = 0; k < inc_i; k++, ii++) {
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
            uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
            rmm_input.get(ii, jj) = fortran_vars.rmm_input_ndens1.data[big_index];
          }
        }
      }
    }
    t_density.pause();

    /******** each point *******/
    uint point = 0;
    for (list<Point>::const_iterator point_it = group.points.begin(); point_it != group.points.end(); ++point_it, ++point)
    {
      t_density.start();
      /** density **/
      float partial_density = 0;

      for (uint i = 0; i < group_m; i++) {
        float w = 0.0f;
        float Fi = group.function_values.data[point * group_m + i];
        for (uint j = i; j < group_m; j++) {
          float Fj = group.function_values.data[point * group_m + j];
          w += rmm_input.data[j * group_m + i] * Fj;
        }
        partial_density += Fi * w;
      }

      /** energy / potential **/
      float exc = 0, corr = 0, y2a = 0;
      if (compute_energy) {
        if (compute_forces) cpu_pot<true, true>(partial_density, exc, corr, y2a);
        else cpu_pot<true, false>(partial_density, exc, corr, y2a);
        total_energy += (partial_density * point_it->weight) * (exc + corr);
      }
      else cpu_pot<false, true>(partial_density, exc, corr, y2a);

      t_density.pause();
      t_resto.start();

      /** compute forces **/
      if (compute_forces) {
        float factor = point_it->weight * y2a;
        
        // compute w
        HostMatrixFloat w(fortran_vars.nco);
        for (uint k = 0; k < fortran_vars.nco; k++)
        {
          float w_local = 0;
          ii = 0;
          for (uint i = 0; i < group.functions.size(); i++)
          {
            uint inc;
            if (i < group.s_functions) inc = 1;
            else if (i < group.s_functions + group.p_functions) inc = 3;
            else inc = 6;
            uint big_i = group.functions[i];
            for (uint j = 0; j < inc; j++, ii++) w_local += group.function_values.data[point * group_m + ii] * fortran_vars.rmm_input.get(big_i + j, k);
          }
          w.data[k] = w_local;
        }
        
        // compute density derivative
        HostMatrixFloat3 density_derivs(fortran_vars.atoms);
        density_derivs.zero();
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
              wrdm += r * w.data[k];
            }
            uint nuc = fortran_vars.nucleii.get(big_i + j) - 1;
            density_derivs.get(nuc) += group.gradient_values.get(ii, point) * wrdm;
          }
        }
        for (set<uint>::const_iterator it = group.nucleii.begin(); it != group.nucleii.end(); ++it)
          forces.get(*it) += density_derivs.get(*it) * factor * 4.0f;
      }

      t_resto.pause();

      t_rmm.start();
      /******** RMM *******/
      if (compute_rmm) {
        float factor = point_it->weight * y2a;
        for (uint i = 0; i < group_m; i++) {
          float Fi = group.function_values.get(i, point);
          for (uint j = i; j < group_m; j++) {
            float Fj = group.function_values.get(j, point);
            rmm_output.get(i, j) += Fi * Fj * factor;
          }
        }
      }

      t_rmm.pause();
    }

    t_rmm.start();
    /* accumulate RMM results for this group */
    if (compute_rmm) {
      uint ii = 0;
      for (uint i = 0; i < group.functions.size(); i++) {
        uint inc_i;
        if (i < group.s_functions) inc_i = 1;
        else if (i < group.s_functions + group.p_functions) inc_i = 3;
        else inc_i = 6;

        for (uint k = 0; k < inc_i; k++, ii++) {
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

              uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
              fortran_vars.rmm_output.get(big_index) += rmm_output.get(ii, jj);
            }
          }
        }
      }
    }
    t_rmm.pause();

    if (compute_forces) {
      FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, FORTRAN_MAX_ATOMS); // TODO: mover esto a init.cpp
      for (uint i = 0; i < fortran_vars.atoms; i++) {
        float3 this_force = forces.get(i);
        //cout << "F: " << forces.get(i).x << " " << forces.get(i).y << " " << forces.get(i).z << endl;
        fort_forces.get(i,0) += this_force.x;
        fort_forces.get(i,1) += this_force.y;
        fort_forces.get(i,2) += this_force.z;        
      }
      forces.zero();
    }
  }

  /***** send results to fortran ****/
  if (compute_energy) *fort_energy_ptr = total_energy;

  timer_total.stop();
  cout << "iteration: " << timer_total << endl;
  cout << "rmm: " << t_rmm << " density: " << t_density << " resto: " << t_resto << endl;
}
#endif
