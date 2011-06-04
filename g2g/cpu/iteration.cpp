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
#include "scalar_vector_types.h"

#include "cpu/pot.h"

using std::cout;
using std::endl;
using std::list;
using namespace G2G;


template <bool compute_rmm, bool lda, bool compute_forces>
void g2g_iteration(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr) {  
  Timers timers;
  timers.total.start();

  partition.solve(timers, compute_rmm, lda, compute_forces, compute_energy, fort_energy_ptr, fort_forces_ptr);

  timers.total.stop();
  cout << timers << endl;  
}

template<class scalar_type>
void PointGroup<scalar_type>::solve(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, bool compute_energy,
                                    double* fort_energy_ptr, double* fort_forces_ptr)
{
  HostMatrix<scalar_type> rmm_output;
  uint group_m = total_functions();
  if (compute_rmm) { rmm_output.resize(group_m, group_m); rmm_output.zero(); }

  #if CPU_RECOMPUTE
  /** Compute functions **/
  timers.functions.start();
  compute_functions<compute_forces, !lda>();
  timers.functions.pause();
  #endif

  // prepare rmm_input for this group
  timers.density.start();
  HostMatrix<scalar_type> rmm_input(group_m, group_m);
  get_rmm_input(rmm_input);
  timers.density.pause();

  HostMatrix<vec_type3> forces(total_nucleii(), 1); forces.zero();
  HostMatrix<vec_type3> dd;

  /******** each point *******/
  uint point = 0;
  for (list<Point>::const_iterator point_it = points.begin(); point_it != points.end(); ++point_it, ++point)
  {
    timers.density.start();

    /** density **/
    scalar_type partial_density = 0;
    vec_type3 dxyz(0,0,0);
    vec_type3 dd1(0,0,0);
    vec_type3 dd2(0,0,0);

    if (lda) {
      for (uint i = 0; i < group_m; i++) {
        float w = 0.0;
        float Fi = function_values(i, point);
        for (uint j = i; j < group_m; j++) {
          scalar_type Fj = function_values(j, point);
          w += rmm_input(j, i) * Fj;
        }
        partial_density += Fi * w;
      }
    }
    else {
      for (uint i = 0; i < group_m; i++) {
        float w = 0.0;
        vec_type3 w3(0,0,0);
        vec_type3 ww1(0,0,0);
        vec_type3 ww2(0,0,0);

        scalar_type Fi = function_values(i, point);
        vec_type3 Fgi(gradient_values(i, point));
        vec_type3 Fhi1(hessian_values(2 * (i + 0) + 0, point));
        vec_type3 Fhi2(hessian_values(2 * (i + 0) + 1, point));

        for (uint j = 0; j <= i; j++) {
          scalar_type rmm = rmm_input(j,i);
          scalar_type Fj = function_values(j, point);
          w += Fj * rmm;

          vec_type3 Fgj(gradient_values(j, point));
          w3 += Fgj * rmm;

          vec_type3 Fhj1(hessian_values(2 * (j + 0) + 0, point));
          vec_type3 Fhj2(hessian_values(2 * (j + 0) + 1, point));
          ww1 += Fhj1 * rmm;
          ww2 += Fhj2 * rmm;
        }
        partial_density += Fi * w;

        dxyz += Fgi * w + w3 * Fi;
        dd1 += Fgi * w3 * 2 + Fhi1 * w + ww1 * Fi;

        vec_type3 FgXXY(Fgi.x(), Fgi.x(), Fgi.y());
        vec_type3 w3YZZ(w3.y(), w3.z(), w3.z());
        vec_type3 FgiYZZ(Fgi.y(), Fgi.z(), Fgi.z());
        vec_type3 w3XXY(w3.x(), w3.x(), w3.y());

        dd2 += FgXXY * w3YZZ + FgiYZZ * w3XXY + Fhi2 * w + ww2 * Fi;
      }
    }
    timers.density.pause();

    timers.forces.start();
    /** density derivatives **/
    if (compute_forces) {
      dd.resize(total_nucleii(), 1); dd.zero();
      for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
        uint nuc = func2local_nuc(ii);
        uint inc_i = small_function_type(i);
        vec_type3 this_dd = vec_type3(0,0,0);
        for (uint k = 0; k < inc_i; k++, ii++) {
          scalar_type w = 0.0;
          for (uint j = 0; j < group_m; j++) {
            scalar_type Fj = function_values(j, point);
            w += rmm_input(j, ii) * Fj * (ii == j ? 2 : 1);
          }
          this_dd -= gradient_values(ii, point) * w;
        }
        dd(nuc) += this_dd;
      }
    }
    timers.forces.pause();

    timers.pot.start();

    timers.density.start();
    /** energy / potential **/
    scalar_type exc = 0, corr = 0, y2a = 0;
    if (lda)
      cpu_pot(partial_density, exc, corr, y2a);
    else {
      cpu_potg(partial_density, dxyz, dd1, dd2, exc, corr, y2a);
    }

    timers.pot.pause();

    if (compute_energy)
      *fort_energy_ptr += (partial_density * point_it->weight) * (exc + corr);
    timers.density.pause();


    /** forces **/
    timers.forces.start();
    if (compute_forces) {
      scalar_type factor = point_it->weight * y2a;
      for (uint i = 0; i < total_nucleii(); i++)
        forces(i) += dd(i) * factor;
    }
    timers.forces.pause();

    /** RMM **/
    timers.rmm.start();
    if (compute_rmm) {
      scalar_type factor = point_it->weight * y2a;
      HostMatrix<scalar_type>::blas_ssyr(LowerTriangle, factor, function_values, rmm_output, point);
    }
    timers.rmm.pause();
  }

  timers.forces.start();
  /* accumulate force results for this group */
  if (compute_forces) {
    FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, FORTRAN_MAX_ATOMS); // TODO: mover esto a init.cpp
    for (uint i = 0; i < total_nucleii(); i++) {
      uint global_atom = local2global_nuc[i];
      vec_type3 this_force = forces(i);
      fort_forces(global_atom,0) += this_force.x();
      fort_forces(global_atom,1) += this_force.y();
      fort_forces(global_atom,2) += this_force.z();
    }
  }
  timers.forces.pause();

  timers.rmm.start();
  /* accumulate RMM results for this group */
  if (compute_rmm) {
    for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
      uint inc_i = small_function_type(i);
      for (uint k = 0; k < inc_i; k++, ii++) {
        uint big_i = local2global_func[i] + k;

        for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
          uint inc_j = small_function_type(j);
          for (uint l = 0; l < inc_j; l++, jj++) {
            uint big_j = local2global_func[j] + l;
            if (big_i > big_j) continue;

            uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
            fortran_vars.rmm_output(big_index) += rmm_output(ii, jj);
          }
        }
      }
    }
  }
  timers.rmm.pause();

  #if CPU_RECOMPUTE
  /* clear functions */
  function_values.deallocate();
  gradient_values.deallocate();
  hessian_values.deallocate();
  #endif
}

template void g2g_iteration<true, true, true>(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<true, true, false>(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<true, false, true>(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<true, false, false>(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<false, true, true>(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<false, true, false>(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<false, false, true>(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr);
template void g2g_iteration<false, false, false>(bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr);
