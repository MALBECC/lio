#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cstring>
#include <omp.h>
#include <cstdio>
#include "../common.h"
#include "../init.h"
#include "../cuda_includes.h"
#include "../matrix.h"
#include "../timer.h"
#include "../partition.h"

#include <stdlib.h>
#include "../pointxc/calc_ggaCS.h"
#include "../pointxc/calc_ggaOS.h"
#include "../pointxc/calc_ldaCS.h"

#if USE_LIBXC
#include <xc.h>
#include "../libxc/libxcproxy.h"
#endif

using std::cout;
using std::endl;
using std::vector;

namespace G2G {

template <class scalar_type>
void PointGroupCPU<scalar_type>::solve(
    Timers& timers, bool compute_rmm, bool lda, bool compute_forces,
    bool compute_energy, double& energy, double& energy_i, double& energy_c,
    double& energy_c1, double& energy_c2, HostMatrix<double>& fort_forces,
    int inner_threads, HostMatrix<double>& rmm_global_output, bool OPEN) {}

template <class scalar_type>
void PointGroupCPU<scalar_type>::solve_closed(
    Timers& timers, bool compute_rmm, bool lda, bool compute_forces,
    bool compute_energy, double& energy, HostMatrix<double>& fort_forces,
    int inner_threads, HostMatrix<double>& rmm_global_output) {
  const uint group_m = this->total_functions();
  const int npoints = this->points.size();

  //printf("solve_closed(...)\n");

#if CPU_RECOMPUTE or !GPU_KERNELS
  /** Compute functions **/
  timers.functions.start();
  compute_functions(compute_forces, !lda);
  timers.functions.pause();
#endif

#if USE_LIBXC
#if LIBXC_CPU
  /** Libxc CPU - version **/
  const int nspin = XC_UNPOLARIZED;
  const int functionalExchange = fortran_vars.ex_functional_id; //101;
  const int functionalCorrelation = fortran_vars.ec_functional_id; // 130;
  LibxcProxy<scalar_type,3> libxcProxy(functionalExchange, functionalCorrelation, nspin);
#endif
#endif

  double localenergy = 0.0;

  // prepare rmm_input for this group
  timers.density.start();

  HostMatrix<scalar_type> rmm_input(group_m, group_m);
  get_rmm_input(rmm_input);

  vector<vec_type3> forces;
  vector<std::vector<vec_type3> > forces_mat;
  HostMatrix<scalar_type> factors_rmm;

  if (compute_rmm || compute_forces) factors_rmm.resize(this->points.size(), 1);

  if (compute_forces) {
    forces.resize(this->total_nucleii(), vec_type3(0.f, 0.f, 0.f));
    forces_mat.resize(
        this->points.size(),
        vector<vec_type3>(this->total_nucleii(), vec_type3(0.f, 0.f, 0.f)));
  }

  const int iexch = fortran_vars.iexch;

  /** density **/
  if (lda) {
    for (int point = 0; point < this->points.size(); point++) {
      scalar_type partial_density = 0;
      for (int i = 0; i < group_m; i++) {
        scalar_type w = 0.0;
        scalar_type Fi = function_values(i, point);
        for (int j = i; j < group_m; j++) {
          scalar_type Fj = function_values(j, point);
          w += rmm_input(j, i) * Fj;
        }
        partial_density += Fi * w;
      }
      scalar_type exc = 0, corr = 0, y2a = 0;

      calc_ldaCS_in(partial_density, exc, corr, y2a, iexch);

      if (compute_energy) {
        localenergy +=
            (partial_density * this->points[point].weight) * (exc + corr);
      }

      /** RMM **/
      if (compute_rmm || compute_forces) {
        factors_rmm(point) = this->points[point].weight * y2a;
      }
    }
  } else {
#pragma omp parallel for num_threads(inner_threads) \
    reduction(+ : localenergy) schedule(static)
    for (int point = 0; point < npoints; point++) {
      scalar_type pd, tdx, tdy, tdz, tdd1x, tdd1y, tdd1z, tdd2x, tdd2y, tdd2z;
      pd = tdx = tdy = tdz = tdd1x = tdd1y = tdd1z = tdd2x = tdd2y = tdd2z =0.0;

      const scalar_type* fv = function_values.row(point);
      const scalar_type* gxv = gX.row(point);
      const scalar_type* gyv = gY.row(point);
      const scalar_type* gzv = gZ.row(point);
      const scalar_type* hpxv = hPX.row(point);
      const scalar_type* hpyv = hPY.row(point);
      const scalar_type* hpzv = hPZ.row(point);
      const scalar_type* hixv = hIX.row(point);
      const scalar_type* hiyv = hIY.row(point);
      const scalar_type* hizv = hIZ.row(point);

      for (int i = 0; i < group_m; i++) {
        scalar_type w = 0;
        scalar_type w3xc = 0, w3yc = 0, w3zc = 0;
        scalar_type ww1xc = 0, ww1yc = 0, ww1zc = 0;
        scalar_type ww2xc = 0, ww2yc = 0, ww2zc = 0;

        const scalar_type* rm = rmm_input.row(i);
#if INTEL_COMP
#pragma vector always
#endif
        for (int j = 0; j <= i; j++) {
          const scalar_type rmj = rm[j];
          w += fv[j] * rmj;
          w3xc += gxv[j] * rmj;
          w3yc += gyv[j] * rmj;
          w3zc += gzv[j] * rmj;
          ww1xc += hpxv[j] * rmj;
          ww1yc += hpyv[j] * rmj;
          ww1zc += hpzv[j] * rmj;
          ww2xc += hixv[j] * rmj;
          ww2yc += hiyv[j] * rmj;
          ww2zc += hizv[j] * rmj;
        }

        const scalar_type Fi = fv[i];
        const scalar_type gx = gxv[i], gy = gyv[i], gz = gzv[i];
        const scalar_type hpx = hpxv[i], hpy = hpyv[i], hpz = hpzv[i];
        const scalar_type hix = hixv[i], hiy = hiyv[i], hiz = hizv[i];

        pd += Fi * w;

        tdx += gx * w + w3xc * Fi;
        tdd1x += gx * w3xc * 2 + hpx * w + ww1xc * Fi;
        tdd2x += gx * w3yc + gy * w3xc + hix * w + ww2xc * Fi;

        tdy += gy * w + w3yc * Fi;
        tdd1y += gy * w3yc * 2 + hpy * w + ww1yc * Fi;
        tdd2y += gx * w3zc + gz * w3xc + hiy * w + ww2yc * Fi;

        tdz += gz * w + w3zc * Fi;
        tdd1z += gz * w3zc * 2 + hpz * w + ww1zc * Fi;
        tdd2z += gy * w3zc + gz * w3yc + hiz * w + ww2zc * Fi;
      }

      /** energy / potential **/
      scalar_type exc = 0, corr = 0, y2a = 0;
      const vec_type3 dxyz(tdx, tdy, tdz);
      const vec_type3 dd1(tdd1x, tdd1y, tdd1z);
      const vec_type3 dd2(tdd2x, tdd2y, tdd2z);

#if USE_LIBXC
    /** Libxc CPU - version **/
#if LIBXC_CPU
    if (fortran_vars.use_libxc) {
	libxcProxy.doGGA(pd, dxyz, dd1, dd2, exc, corr, y2a);
    } else {
	calc_ggaCS_in<scalar_type, 3>(pd, dxyz, dd1, dd2, exc, corr, y2a, iexch);
    }
#else
    calc_ggaCS_in<scalar_type, 3>(pd, dxyz, dd1, dd2, exc, corr, y2a, iexch);
#endif
#else
      calc_ggaCS_in<scalar_type, 3>(pd, dxyz, dd1, dd2, exc, corr, y2a, iexch);
#endif

      const scalar_type wp = this->points[point].weight;

      if (compute_energy) {
        localenergy += (pd * wp) * (exc + corr);
      }

      /** RMM **/
      if (compute_rmm || compute_forces) {
        factors_rmm(point) = wp * y2a;
      }
    }
  }
  timers.density.pause();

  timers.forces.start();
  if (compute_forces) {
    HostMatrix<scalar_type> ddx, ddy, ddz;
    ddx.resize(this->total_nucleii(), 1);
    ddy.resize(this->total_nucleii(), 1);
    ddz.resize(this->total_nucleii(), 1);
#pragma omp parallel for num_threads(inner_threads)
    for (int point = 0; point < this->points.size(); point++) {
      ddx.zero();
      ddy.zero();
      ddz.zero();
      for (int i = 0, ii = 0; i < this->total_functions_simple(); i++) {
        uint nuc = this->func2local_nuc(ii);
        uint inc_i = this->small_function_type(i);
        scalar_type tddx = 0, tddy = 0, tddz = 0;
        for (uint k = 0; k < inc_i; k++, ii++) {
          scalar_type w = 0.0;
          for (uint j = 0; j < group_m; j++) {
            scalar_type Fj = function_values(j, point);
            w += rmm_input(j, ii) * Fj * (ii == j ? 2 : 1);
          }
          tddx -= w * gX(ii, point);
          tddy -= w * gY(ii, point);
          tddz -= w * gZ(ii, point);
        }
        ddx(nuc) += tddx;
        ddy(nuc) += tddy;
        ddz(nuc) += tddz;
      }
      scalar_type factor = factors_rmm(point);
      for (int i = 0; i < this->total_nucleii(); i++) {
        forces_mat[point][i] = vec_type3(ddx(i), ddy(i), ddz(i)) * factor;
      }
    }
    /* accumulate forces for each point */
    if (forces_mat.size() > 0) {
#pragma omp parallel for num_threads(inner_threads) schedule(static)
      for (int j = 0; j < forces_mat[0].size(); j++) {
        vec_type3 acum(0.f, 0.f, 0.f);
        for (int i = 0; i < forces_mat.size(); i++) {
          acum += forces_mat[i][j];
        }
        forces[j] = acum;
      }
    }
/* accumulate force results for this group */
#pragma omp parallel for num_threads(inner_threads)
    for (int i = 0; i < this->total_nucleii(); i++) {
      uint global_atom = this->local2global_nuc[i];
      vec_type3 this_force = forces[i];
      fort_forces(global_atom, 0) += this_force.x;
      fort_forces(global_atom, 1) += this_force.y;
      fort_forces(global_atom, 2) += this_force.z;
    }
  }
  timers.forces.pause();

  timers.rmm.start();
  /* accumulate RMM results for this group */
  if (compute_rmm) {
    const int indexes = this->rmm_bigs.size();
#pragma omp parallel for num_threads(inner_threads) schedule(static)
    for (int i = 0; i < indexes; i++) {
      int bi = this->rmm_bigs[i], row = this->rmm_rows[i],
          col = this->rmm_cols[i];

      double res = 0;
      const scalar_type* fvr = function_values_transposed.row(row);
      const scalar_type* fvc = function_values_transposed.row(col);

#if INTEL_COMP
#pragma vector always
#endif
      for (int point = 0; point < npoints; point++) {
        res += fvr[point] * fvc[point] * factors_rmm(point);
      }
      rmm_global_output(bi) += res;
    }
  }
  timers.rmm.pause();

  energy += localenergy;

#if CPU_RECOMPUTE or !GPU_KERNELS
  /* clear functions */
  gX.deallocate();
  gY.deallocate();
  gZ.deallocate();
  hIX.deallocate();
  hIY.deallocate();
  hIZ.deallocate();
  hPX.deallocate();
  hPY.deallocate();
  hPZ.deallocate();
  function_values_transposed.deallocate();
#endif
}

template <class scalar_type>
void PointGroupCPU<scalar_type>::solve_opened(
    Timers& timers, bool compute_rmm, bool lda, bool compute_forces,
    bool compute_energy, double& energy, double& energy_i, double& energy_c,
    double& energy_c1, double& energy_c2, HostMatrix<double>& fort_forces,
    HostMatrix<double>& rmm_output_local_a,
    HostMatrix<double>& rmm_output_local_b) {
  //   std::exit(0);
  int inner_threads = 1;
  const uint group_m = this->total_functions();
  const int npoints = this->points.size();

#if CPU_RECOMPUTE or !GPU_KERNELS
  /** Compute functions **/
  timers.functions.start();
  compute_functions(compute_forces, !lda);
  timers.functions.pause();
#endif

  double localenergy = 0.0;

  // prepare rmm_input for this group
  timers.density.start();

  HostMatrix<scalar_type> rmm_input_a(group_m, group_m),
      rmm_input_b(group_m, group_m);
  get_rmm_input(rmm_input_a, rmm_input_b);

  vector<vec_type3> forces_a, forces_b;
  vector<std::vector<vec_type3> > forces_mat_a, forces_mat_b;
  HostMatrix<scalar_type> factors_rmm_a, factors_rmm_b;

  if (compute_rmm || compute_forces) {
    factors_rmm_a.resize(this->points.size(), 1);
    factors_rmm_b.resize(this->points.size(), 1);
  }

  if (compute_forces) {
    forces_a.resize(this->total_nucleii(), vec_type3(0.f, 0.f, 0.f));
    forces_b.resize(this->total_nucleii(), vec_type3(0.f, 0.f, 0.f));
    forces_mat_a.resize(
        this->points.size(),
        vector<vec_type3>(this->total_nucleii(), vec_type3(0.f, 0.f, 0.f)));
    forces_mat_b.resize(
        this->points.size(),
        vector<vec_type3>(this->total_nucleii(), vec_type3(0.f, 0.f, 0.f)));
  }

  const int iexch = fortran_vars.iexch;

  /** density **/
  if (lda) {
  } else {
#pragma omp parallel for num_threads(inner_threads) \
    reduction(+ : localenergy) schedule(static)
    for (int point = 0; point < npoints; point++) {
      // Partial density and space derivatives values in point (accumulates
      // functions)
      scalar_type pd_a, pd_b;
      scalar_type tdx_a, tdy_a, tdz_a, tdx_b, tdy_b, tdz_b;  // dx, dy, dz
      scalar_type tdd1x_a, tdd1y_a, tdd1z_a, tdd1x_b, tdd1y_b,
          tdd1z_b;  // dxx, dyy, dzz
      scalar_type tdd2x_a, tdd2y_a, tdd2z_a, tdd2x_b, tdd2y_b,
          tdd2z_b;  // dxy, dxz, dyz

      pd_a = tdx_a = tdy_a = tdz_a = tdd1x_a = tdd1y_a = tdd1z_a = tdd2x_a =
          tdd2y_a = tdd2z_a = 0;
      pd_b = tdx_b = tdy_b = tdz_b = tdd1x_b = tdd1y_b = tdd1z_b = tdd2x_b =
          tdd2y_b = tdd2z_b = 0;

      // Evaluated basis functions, derivatives and gradient for the point
      const scalar_type* fv = function_values.row(point);
      const scalar_type* gxv = gX.row(point);    // dx
      const scalar_type* gyv = gY.row(point);    // dy
      const scalar_type* gzv = gZ.row(point);    // dz
      const scalar_type* hpxv = hPX.row(point);  // dxx
      const scalar_type* hpyv = hPY.row(point);  // dyy
      const scalar_type* hpzv = hPZ.row(point);  // dzz
      const scalar_type* hixv = hIX.row(point);  // dxy
      const scalar_type* hiyv = hIY.row(point);  // dxz
      const scalar_type* hizv = hIZ.row(point);  // dyz

      for (int i = 0; i < group_m; i++) {
        // Accumulates j functions weighted for i function and point
        scalar_type w_a, w_b;
        scalar_type w3x_a, w3y_a, w3z_a, w3x_b, w3y_b, w3z_b;
        scalar_type ww1x_a, ww1y_a, ww1z_a, ww1x_b, ww1y_b, ww1z_b;
        scalar_type ww2x_a, ww2y_a, ww2z_a, ww2x_b, ww2y_b, ww2z_b;
        w_a = w3x_a = w3y_a = w3z_a = ww1x_a = ww1y_a = ww1z_a = ww2x_a =
            ww2y_a = ww2z_a = 0.0;
        w_b = w3x_b = w3y_b = w3z_b = ww1x_b = ww1y_b = ww1z_b = ww2x_b =
            ww2y_b = ww2z_b = 0.0;

        const scalar_type* rmm_a = rmm_input_a.row(i);
        const scalar_type* rmm_b = rmm_input_b.row(i);
#if INTEL_COMP
#pragma vector always
#endif
        for (int j = 0; j <= i; j++) {
          const scalar_type rmm_aj = rmm_a[j];
          const scalar_type rmm_bj = rmm_b[j];
          w_a += fv[j] * rmm_aj;
          w3x_a += gxv[j] * rmm_aj;
          w3y_a += gyv[j] * rmm_aj;
          w3z_a += gzv[j] * rmm_aj;
          ww1x_a += hpxv[j] * rmm_aj;
          ww1y_a += hpyv[j] * rmm_aj;
          ww1z_a += hpzv[j] * rmm_aj;
          ww2x_a += hixv[j] * rmm_aj;
          ww2y_a += hiyv[j] * rmm_aj;
          ww2z_a += hizv[j] * rmm_aj;

          w_b += fv[j] * rmm_bj;
          w3x_b += gxv[j] * rmm_bj;
          w3y_b += gyv[j] * rmm_bj;
          w3z_b += gzv[j] * rmm_bj;
          ww1x_b += hpxv[j] * rmm_bj;
          ww1y_b += hpyv[j] * rmm_bj;
          ww1z_b += hpzv[j] * rmm_bj;
          ww2x_b += hixv[j] * rmm_bj;
          ww2y_b += hiyv[j] * rmm_bj;
          ww2z_b += hizv[j] * rmm_bj;
        }

        const scalar_type Fi = fv[i];
        const scalar_type gx = gxv[i], gy = gyv[i], gz = gzv[i];
        const scalar_type hpx = hpxv[i], hpy = hpyv[i], hpz = hpzv[i];
        const scalar_type hix = hixv[i], hiy = hiyv[i], hiz = hizv[i];

        pd_a += Fi * w_a;

        tdx_a += gx * w_a + w3x_a * Fi;
        tdd1x_a += gx * w3x_a * 2 + hpx * w_a + ww1x_a * Fi;
        tdd2x_a += gx * w3y_a + gy * w3x_a + hix * w_a + ww2x_a * Fi;

        tdy_a += gy * w_a + w3y_a * Fi;
        tdd1y_a += gy * w3y_a * 2 + hpy * w_a + ww1y_a * Fi;
        tdd2y_a += gx * w3z_a + gz * w3x_a + hiy * w_a + ww2y_a * Fi;

        tdz_a += gz * w_a + w3z_a * Fi;
        tdd1z_a += gz * w3z_a * 2 + hpz * w_a + ww1z_a * Fi;
        tdd2z_a += gy * w3z_a + gz * w3y_a + hiz * w_a + ww2z_a * Fi;

        pd_b += Fi * w_b;

        tdx_b += gx * w_b + w3x_b * Fi;
        tdd1x_b += gx * w3x_b * 2 + hpx * w_b + ww1x_b * Fi;
        tdd2x_b += gx * w3y_b + gy * w3x_b + hix * w_b + ww2x_b * Fi;

        tdy_b += gy * w_b + w3y_b * Fi;
        tdd1y_b += gy * w3y_b * 2 + hpy * w_b + ww1y_b * Fi;
        tdd2y_b += gx * w3z_b + gz * w3x_b + hiy * w_b + ww2y_b * Fi;

        tdz_b += gz * w_b + w3z_b * Fi;
        tdd1z_b += gz * w3z_b * 2 + hpz * w_b + ww1z_b * Fi;
        tdd2z_b += gy * w3z_b + gz * w3y_b + hiz * w_b + ww2z_b * Fi;
      }

      /** energy / potential **/
      scalar_type exc_corr = 0, corr1 = 0, corr2 = 0;
      scalar_type exc = 0, corr = 0, y2a = 0, y2b = 0;
      const vec_type3 dxyz_a(tdx_a, tdy_a, tdz_a), dxyz_b(tdx_b, tdy_b, tdz_b);
      const vec_type3 dd1_a(tdd1x_a, tdd1y_a, tdd1z_a),
          dd1_b(tdd1x_b, tdd1y_b, tdd1z_b);
      const vec_type3 dd2_a(tdd2x_a, tdd2y_a, tdd2z_a),
          dd2_b(tdd2x_b, tdd2y_b, tdd2z_b);

      calc_ggaOS<scalar_type, 3>(pd_a, pd_b, dxyz_a, dxyz_b, dd1_a, dd1_b,
                                 dd2_a, dd2_b, exc_corr, exc, corr, corr1,
                                 corr2, y2a, y2b, 9);

      const scalar_type wp = this->points[point].weight;

      if (compute_energy) {
        localenergy += ((pd_a + pd_b) * wp) * (exc + corr);
      }

      /** RMM **/
      if (compute_rmm || compute_forces) {
        factors_rmm_a(point) = wp * y2a;
        factors_rmm_b(point) = wp * y2b;
      }
    }
  }
  timers.density.pause();

  timers.forces.start();
  if (compute_forces) {
    HostMatrix<scalar_type> ddx_a, ddy_a, ddz_a;
    HostMatrix<scalar_type> ddx_b, ddy_b, ddz_b;
    ddx_a.resize(this->total_nucleii(), 1);
    ddx_b.resize(this->total_nucleii(), 1);
    ddy_a.resize(this->total_nucleii(), 1);
    ddy_b.resize(this->total_nucleii(), 1);
    ddz_a.resize(this->total_nucleii(), 1);
    ddz_b.resize(this->total_nucleii(), 1);

#pragma omp parallel for num_threads(inner_threads)
    for (int point = 0; point < this->points.size(); point++) {
      ddx_a.zero();
      ddy_a.zero();
      ddz_a.zero();
      ddx_b.zero();
      ddy_b.zero();
      ddz_b.zero();

      for (int i = 0, ii = 0; i < this->total_functions_simple(); i++) {
        uint nuc = this->func2local_nuc(ii);
        uint inc_i = this->small_function_type(i);
        scalar_type tddx_a = 0, tddy_a = 0, tddz_a = 0;
        scalar_type tddx_b = 0, tddy_b = 0, tddz_b = 0;
        for (uint k = 0; k < inc_i; k++, ii++) {
          scalar_type w_a = 0.0, w_b = 0.0;
          for (uint j = 0; j < group_m; j++) {
            scalar_type Fj = function_values(j, point);
            w_a += rmm_input_a(j, ii) * Fj * (ii == j ? 2 : 1);
            w_b += rmm_input_b(j, ii) * Fj * (ii == j ? 2 : 1);
          }
          tddx_a -= w_a * gX(ii, point);
          tddx_b -= w_b * gX(ii, point);
          tddy_a -= w_a * gY(ii, point);
          tddy_b -= w_b * gY(ii, point);
          tddz_a -= w_a * gZ(ii, point);
          tddz_b -= w_b * gZ(ii, point);
        }
        ddx_a(nuc) += tddx_a;
        ddy_a(nuc) += tddy_a;
        ddz_a(nuc) += tddz_a;
        ddx_b(nuc) += tddx_b;
        ddy_b(nuc) += tddy_b;
        ddz_b(nuc) += tddz_b;
      }

      scalar_type factor_a = factors_rmm_a(point);
      scalar_type factor_b = factors_rmm_b(point);
      for (int i = 0; i < this->total_nucleii(); i++) {
        forces_mat_a[point][i] =
            vec_type3(ddx_a(i), ddy_a(i), ddz_a(i)) * factor_a;
        forces_mat_b[point][i] =
            vec_type3(ddx_b(i), ddy_b(i), ddz_b(i)) * factor_b;
      }
    }

    /* accumulate forces for each point */
    if ((forces_mat_a.size() > 0) && (forces_mat_b.size() > 0)) {
#pragma omp parallel for num_threads(inner_threads) schedule(static)
      for (int j = 0; j < forces_mat_a[0].size(); j++) {
        vec_type3 acum_a(0.f, 0.f, 0.f);
        vec_type3 acum_b(0.f, 0.f, 0.f);
        for (int i = 0; i < forces_mat_a.size(); i++) {
          acum_a += forces_mat_a[i][j];
          acum_b += forces_mat_b[i][j];
        }
        forces_a[j] = acum_a;
        forces_b[j] = acum_b;
      }
    }

/* accumulate force results for this group */
#pragma omp parallel for num_threads(inner_threads)
    for (int i = 0; i < this->total_nucleii(); i++) {
      uint global_atom = this->local2global_nuc[i];
      vec_type3 this_force = forces_a[i] + forces_b[i];
      fort_forces(global_atom, 0) += this_force.x;
      fort_forces(global_atom, 1) += this_force.y;
      fort_forces(global_atom, 2) += this_force.z;
    }
  }
  timers.forces.pause();

  timers.rmm.start();
  /* accumulate RMM results for this group */
  if (compute_rmm) {
    const int indexes = this->rmm_bigs.size();
#pragma omp parallel for num_threads(inner_threads) schedule(static)
    for (int i = 0; i < indexes; i++) {
      int bi = this->rmm_bigs[i], row = this->rmm_rows[i],
          col = this->rmm_cols[i];

      double res_a = 0, res_b = 0;
      const scalar_type* fvr = function_values_transposed.row(row);
      const scalar_type* fvc = function_values_transposed.row(col);

#if INTEL_COMP
#pragma vector always
#endif
      for (int point = 0; point < npoints; point++) {
        res_a += fvr[point] * fvc[point] * factors_rmm_a(point);
        res_b += fvr[point] * fvc[point] * factors_rmm_b(point);
      }
      rmm_output_local_a(bi) += res_a;
      rmm_output_local_b(bi) += res_b;
    }
  }
  timers.rmm.pause();

  energy += localenergy;

#if CPU_RECOMPUTE or !GPU_KERNELS
  /* clear functions */
  gX.deallocate();
  gY.deallocate();
  gZ.deallocate();
  hIX.deallocate();
  hIY.deallocate();
  hIZ.deallocate();
  hPX.deallocate();
  hPY.deallocate();
  hPZ.deallocate();
  function_values_transposed.deallocate();
#endif
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}
