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
  bool compute_rmm = (computation_type == COMPUTE_RMM);

  double total_energy = 0;

  HostMatrixFloat3 dd, forces;
  if (compute_forces) { dd.resize(fortran_vars.atoms, 1); forces.resize(fortran_vars.atoms, 1); forces.zero(); }

  HostMatrixFloat rmm_output;

  Timer t_rmm, t_density, t_resto;

  /********** iterate all groups ***********/
	for (list<PointGroup>::const_iterator group_it = final_partition.begin(); group_it != final_partition.end(); ++group_it)
  {
		const PointGroup& group = *group_it;
    uint group_m = group.total_functions();
    if (compute_rmm) { rmm_output.resize(group_m, group_m); rmm_output.zero(); }

    // prepare rmm_input for this group
    t_density.start();
    HostMatrixFloat rmm_input(group_m, group_m);
    rmm_input.zero();
    for (uint i = 0, ii = 0; i < group.functions.size(); i++) {
      uint inc_i = group.small_function_type(i);

      for (uint k = 0; k < inc_i; k++, ii++) {
        uint big_i = group.functions[i] + k;
        for (uint j = 0, jj = 0; j < group.functions.size(); j++) {
          uint inc_j = group.small_function_type(j);

          for (uint l = 0; l < inc_j; l++, jj++) {
            uint big_j = group.functions[j] + l;
            if (big_i > big_j) continue;
            uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
            if (isnan(fortran_vars.rmm_input_ndens1.data[big_index])) { cout << big_i << " " << big_j << endl; throw runtime_error("bleh"); }
            rmm_input.get(ii, jj) = fortran_vars.rmm_input_ndens1.data[big_index];
            rmm_input.get(jj, ii) = rmm_input.get(ii, jj);
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
      float3 dxyz = make_float3(0,0,0);
      float3 dd1 = make_float3(0,0,0);
      float3 dd2 = make_float3(0,0,0);

      if (fortran_vars.lda) {
        for (uint i = 0; i < group_m; i++) {
          float w = 0.0f;
          float Fi = group.function_values.data[point * group_m + i];
          for (uint j = i; j < group_m; j++) {
            float Fj = group.function_values.data[point * group_m + j];
            w += rmm_input.data[i * group_m + j] * Fj;
          }
          partial_density += Fi * w;
        }
      }
      else {
        for (uint i = 0; i < group_m; i++) {
          float w = 0.0f;
          float3 w3 = make_float3(0,0,0);
          float3 ww1 = make_float3(0,0,0);
          float3 ww2 = make_float3(0,0,0);

          float Fi = group.function_values.data[point * group_m + i];
          float3 Fgi = group.gradient_values.data[point * group_m + i];
          float3 Fhi1 = group.hessian_values.get(2 * (i + 0) + 0, point);
          float3 Fhi2 = group.hessian_values.get(2 * (i + 0) + 1, point);
          if (isnan(Fhi1.x) || isnan(Fhi1.y) || isnan(Fhi1.z)) throw runtime_error("kaspuf");
          if (isnan(Fhi2.x) || isnan(Fhi2.y) || isnan(Fhi2.z)) throw runtime_error("kaspuf");

          for (uint j = i; j < group_m; j++) {
            float rmm = rmm_input.data[i * group_m + j];
            float Fj = group.function_values.data[point * group_m + j];
            w += Fj * rmm;

            float3 Fgj = group.gradient_values.data[point * group_m + j];
            w3 += Fgj * rmm;

            float3 Fhj1 = group.hessian_values.get(2 * (j + 0) + 0, point);
            float3 Fhj2 = group.hessian_values.get(2 * (j + 0) + 1, point);
            ww1 += Fhj1 * rmm;
            ww2 += Fhj2 * rmm;
          }
          partial_density += Fi * w;
          //cout << w << endl;
          if (isnan(Fgi.x) || isnan(Fgi.y) || isnan(Fgi.z)) throw runtime_error("kaspuf");
          if (isnan(w3.x) || isnan(w3.y) || isnan(w3.z)) throw runtime_error("kaspuf");
          if (isnan(w) || isnan(Fi)) throw runtime_error("kaspuf");

          dxyz += Fgi * w + w3 * Fi;
          dd1 += Fgi * w3 * 2 + Fhi1 * w + ww1 * Fi;
          dd2.x += Fgi.x * w3.y + Fgi.y * w3.x + Fhi2.x * w + ww2.x * Fi;
          dd2.y += Fgi.x * w3.z + Fgi.z * w3.x + Fhi2.y * w + ww2.y * Fi;
          dd2.z += Fgi.y * w3.z + Fgi.z * w3.y + Fhi2.z * w + ww2.z * Fi;

          if (isnan(dd1.x) || isnan(dd1.y) || isnan(dd1.z)) throw runtime_error("kaspuf");
          if (isnan(dd2.x) || isnan(dd2.y) || isnan(dd2.z)) throw runtime_error("kaspuf");
          if (isnan(dxyz.x) || isnan(dxyz.y) || isnan(dxyz.z)) throw runtime_error("kaspuf");
        }
      }

      //cout << partial_density << endl;//" " << exc << " " << corr << " " << y2a << endl;

      /** density derivatives **/
      if (compute_forces) {
        dd.zero();
        for (uint i = 0, ii = 0; i < group.functions.size(); i++) {
          uint nuc = group.nuc_map[i];
          uint inc_i = group.small_function_type(i);
          float3 this_dd = make_float3(0,0,0);
          for (uint k = 0; k < inc_i; k++, ii++) {
            float w = 0.0f;
            for (uint j = 0; j < group_m; j++) {
              float Fj = group.function_values.data[point * group_m + j];
              w += rmm_input.data[ii * group_m + j] * Fj * (ii == j ? 2 : 1);
            }
            this_dd -= group.gradient_values.data[point * group_m + ii] * w;
          }
          dd.get(nuc) += this_dd;
        }
      }

      /** energy / potential **/
      float exc = 0, corr = 0, y2a = 0;
      if (fortran_vars.lda)
        cpu_pot(partial_density, exc, corr, y2a);
      else {
        //cout << "antes: " << partial_density << " " << dxyz << " " << dd1 << " " << dd2 << endl;
        cpu_potg(partial_density, dxyz, dd1, dd2, exc, corr, y2a);
        cout << exc << " " << corr << " " << y2a << endl;
      }

      if (compute_energy)
        total_energy += (partial_density * point_it->weight) * (exc + corr);

      t_density.pause();
      t_resto.start();

      if (compute_forces) {
        float factor = point_it->weight * y2a;
        for (set<uint>::const_iterator it = group.nucleii.begin(); it != group.nucleii.end(); ++it)
          forces.get(*it) += dd.data[*it] * factor;
      }        

      t_resto.pause();

      t_rmm.start();
      /******** RMM *******/
      if (compute_rmm) {
        float factor = point_it->weight * y2a;
        for (uint i = 0; i < group_m; i++) {
          float Fi = group.function_values.data[point * group_m + i];
          for (uint j = i; j < group_m; j++) {
            float Fj = group.function_values.data[point * group_m + j];
            rmm_output.data[j * group_m + i] += Fi * Fj * factor;
          }
        }
      }

      t_rmm.pause();
    }

    t_rmm.start();
    /* accumulate RMM results for this group */
    if (compute_rmm) {
      for (uint i = 0, ii = 0; i < group.functions.size(); i++) {
        uint inc_i = group.small_function_type(i);
        for (uint k = 0; k < inc_i; k++, ii++) {
          uint big_i = group.functions[i] + k;

          for (uint j = 0, jj = 0; j < group.functions.size(); j++) {
            uint inc_j = group.small_function_type(j);
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
  }

  if (compute_forces) {
    FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, FORTRAN_MAX_ATOMS); // TODO: mover esto a init.cpp
    for (uint i = 0; i < fortran_vars.atoms; i++) {
      float3 this_force = forces.get(i);
      fort_forces.get(i,0) += this_force.x;
      fort_forces.get(i,1) += this_force.y;
      fort_forces.get(i,2) += this_force.z;
    }
  }

  /***** send results to fortran ****/
  if (compute_energy) *fort_energy_ptr = total_energy;

  timer_total.stop();
  cout << "iteration: " << timer_total << endl;
  cout << "rmm: " << t_rmm << " density: " << t_density << " resto: " << t_resto << endl;
}
#endif
