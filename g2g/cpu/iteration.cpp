#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cstring>
#include <omp.h>
#include <cstdio>
#include "mkl.h"
#include "../buffer_pool.h"
#include "../common.h"
#include "../init.h"
#include "../cuda_includes.h"
#include "../matrix.h"
#include "../timer.h"
#include "../partition.h"

#include "cpu/pot.h"

using std::cout;
using std::endl;
using std::vector;

namespace G2G { 

void trmm(const MKL_INT M, const MKL_INT N, const float *A, float *B) {
  return cblas_strmm(CblasRowMajor,CblasRight,CblasLower,CblasNoTrans, CblasNonUnit, 
      M, N, 1.0, A, N, B, N);
}

void trmm(const MKL_INT M, const MKL_INT N, const double *A, double *B) {
  return cblas_dtrmm(CblasRowMajor,CblasRight,CblasLower,CblasNoTrans, CblasNonUnit, 
      M, N, 1.0, A, N, B, N);
}

template<class scalar_type>
void do_trmm(Timers & ts, const HostMatrix<scalar_type> & triagmat, const HostMatrix<scalar_type> & genmat, scalar_type * res) {
  ts.memcpy.start();
  genmat.copy_to_tmp(res);
  ts.memcpy.pause();
  ts.trmms.start();
  trmm(genmat.height, genmat.width, triagmat.asArray(), res);
  ts.trmms.pause();
}

template<class scalar_type>
void PointGroup<scalar_type>::do_trmms(Timers & ts, ThreadBufferPool<scalar_type> & pool, const HostMatrix<scalar_type> & rmm_input, int threads) const {
  pool.reset();
  ts.memcpy.start();
  scalar_type * pointers[10];
  for(int i = 0; i < 10; i++) {
    pointers[i] = pool.get_pool(pool_elements());
  }
  const HostMatrix<scalar_type> * matrices[] = { 
    &function_values, &gX, &gY, &gZ, &hPX, &hPY, &hPZ, &hIX, &hIY, &hIZ 
  };

  int movethreads = std::min(threads, 10);
  #pragma omp parallel for num_threads(movethreads)
  for(int i = 0; i < 10; i++) {
    matrices[i]->copy_to_tmp(pointers[i]);
  }
  ts.memcpy.pause();
  ts.trmms.start();
  mkl_set_num_threads_local(threads);
  trmm(10 * number_of_points, total_functions(), rmm_input.asArray(), pool.pool_start());
  ts.trmms.pause();
  pool.reset();
}

template<class scalar_type> void PointGroup<scalar_type>::solve(Timers& timers,
  bool compute_rmm, bool lda, bool compute_forces, bool compute_energy, 
  double& energy, double& energy_i, double& energy_c, double& energy_c1, double& energy_c2,
  HostMatrix<double> & fort_forces, ThreadBufferPool<scalar_type> & pool, 
  int inner_threads, HostMatrix<scalar_type> & rmm_global_output, bool OPEN)
{

  solve_closed(timers, compute_rmm, lda, compute_forces, compute_energy, 
    energy, fort_forces, pool, inner_threads, rmm_global_output);
}

template<class scalar_type> void PointGroup<scalar_type>::solve_closed(Timers& timers,
  bool compute_rmm, bool lda, bool compute_forces, bool compute_energy, 
  double& energy, HostMatrix<double> & fort_forces, ThreadBufferPool<scalar_type> & pool, 
  int inner_threads, HostMatrix<scalar_type> & rmm_global_output)
{
  if(points.size() < 32) inner_threads = 1;
  const uint group_m = total_functions();

  #if CPU_RECOMPUTE
  /** Compute functions **/
  timers.functions.start();
  compute_functions(compute_forces, !lda);
  timers.functions.pause();
  #endif

  double localenergy = 0.0;

  // prepare rmm_input for this group
  timers.density.start();

  timers.rmm_input.start();
  HostMatrix<scalar_type> rmm_input(group_m, group_m);
  get_rmm_input(rmm_input);
  timers.rmm_input.pause();

  timers.density.pause();

  vector<vec_type3> forces;
  vector< std::vector<vec_type3> > forces_mat;
  vector<scalar_type> factors_rmm, factors_y2a;
      
  if(compute_rmm) factors_rmm.resize(points.size(), 0.0);

  if(compute_forces) {
    factors_y2a.resize(points.size());
    forces.resize(total_nucleii(), vec_type3(0.f,0.f,0.f)); 
    forces_mat.resize(points.size(), vector<vec_type3>(total_nucleii(), vec_type3(0.f,0.f,0.f)));
  }

  timers.density.start();

  const int iexch = fortran_vars.iexch;
  timers.density_calcs.start();

  /** density **/
  if (lda) {
    for(int point = 0; point < points.size(); point++) {
      scalar_type partial_density = 0;
      for (uint i = 0; i < group_m; i++) {
        scalar_type w = 0.0;
        scalar_type Fi = function_values(i, point);
        for (uint j = i; j < group_m; j++) {
          scalar_type Fj = function_values(j, point);
          w += rmm_input(j, i) * Fj;
        }
        partial_density += Fi * w;
      }
      scalar_type exc = 0, corr = 0, y2a = 0;
      cpu_pot(partial_density, exc, corr, y2a, iexch);

      if (compute_energy) {
          localenergy +=  (partial_density * points[point].weight) * (exc + corr);
      }

      /** forces **/
      if (compute_forces) {
        factors_y2a[point] = y2a;
      }

      /** RMM **/
      if (compute_rmm) {
        factors_rmm[point] = points[point].weight * y2a;
      }
    }
  } else {
    timers.density.start();

    do_trmms(timers, pool, rmm_input, inner_threads);

    int elements = pool_elements();
    scalar_type * wv = pool.get_pool(elements);
    scalar_type * w3x = pool.get_pool(elements); 
    scalar_type * w3y = pool.get_pool(elements); 
    scalar_type * w3z = pool.get_pool(elements);
    scalar_type * ww1x = pool.get_pool(elements); 
    scalar_type * ww1y = pool.get_pool(elements); 
    scalar_type * ww1z = pool.get_pool(elements);
    scalar_type * ww2x = pool.get_pool(elements); 
    scalar_type * ww2y = pool.get_pool(elements); 
    scalar_type * ww2z = pool.get_pool(elements);

    timers.density.pause();

    #pragma omp parallel for num_threads(inner_threads)
    for(int point = 0; point < points.size(); point++) {
      scalar_type pd, tdx, tdy, tdz, tdd1x, tdd1y, tdd1z,tdd2x, tdd2y, tdd2z;
      pd = tdx = tdy = tdz = tdd1x = tdd1y = tdd1z = tdd2x = tdd2y = tdd2z = 0;

      const scalar_type * fv = function_values.row(point);
      const scalar_type * gxv = gX.row(point);
      const scalar_type * gyv = gY.row(point);
      const scalar_type * gzv = gZ.row(point);
      const scalar_type * hpxv = hPX.row(point);
      const scalar_type * hpyv = hPY.row(point);
      const scalar_type * hpzv = hPZ.row(point);
      const scalar_type * hixv = hIX.row(point);
      const scalar_type * hiyv = hIY.row(point);
      const scalar_type * hizv = hIZ.row(point);

      #pragma vector aligned always nontemporal
      for (int i = 0; i < group_m; i++) {
        const int ai = point * group_m + i;
        const scalar_type w = wv[ai];
        const scalar_type Fi = fv[i];
        const scalar_type gx = gxv[i], gy = gyv[i], gz = gzv[i];
        const scalar_type hpx = hpxv[i], hpy = hpyv[i], hpz = hpzv[i];
        const scalar_type hix = hixv[i], hiy = hiyv[i], hiz = hizv[i];
        const scalar_type w3xc = w3x[ai], w3yc = w3y[ai], w3zc = w3z[ai];
        const scalar_type ww1xc = ww1x[ai], ww1yc = ww1y[ai], ww1zc = ww1z[ai];
        const scalar_type ww2xc = ww2x[ai], ww2yc = ww2y[ai], ww2zc = ww2z[ai];

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
      vec_type3 dxyz(tdx,tdy,tdz);
      vec_type3 dd1(tdd1x,tdd1y,tdd1z);
      vec_type3 dd2(tdd2x,tdd2y,tdd2z);

      cpu_potg(pd, dxyz, dd1, dd2, exc, corr, y2a, iexch);

      if (compute_energy) {
        localenergy +=  (pd * points[point].weight) * (exc + corr);
      }

      /** forces **/
      if (compute_forces) {
        factors_y2a[point] = y2a;
      }

      /** RMM **/
      if (compute_rmm) {
        factors_rmm[point] = points[point].weight * y2a;
      }
    }
  }
  timers.density_calcs.pause();

  timers.forces.start();
  if (compute_forces) {
    HostMatrix<vec_type3> dd; 
    dd.resize(total_nucleii(), 1); 
    #pragma omp parallel for num_threads(inner_threads)
    for(int point = 0; point < points.size(); point++) {
      dd.zero();
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
          this_dd -= vec_type3(w*gX(ii, point),w*gY(ii,point),w*gZ(ii,point));
        }
        dd(nuc) += this_dd;
      }
      scalar_type factor = points[point].weight * factors_y2a[point];
      for (uint i = 0; i < total_nucleii(); i++) {
        forces_mat[point][i] = dd(i) * factor;
      }
    }
    /* accumulate forces for each point */
    if(forces_mat.size() > 0) {
      for (int j = 0; j < forces_mat[0].size(); j++) {
        vec_type3 acum(0.f,0.f,0.f);
        for (int i = 0; i < forces_mat.size(); i++) {
          acum += forces_mat[i][j];
        }
        forces[j] = acum;
      }
    }
    /* accumulate force results for this group */
    for (uint i = 0; i < total_nucleii(); i++) {
      uint global_atom = local2global_nuc[i];
      vec_type3 this_force = forces[i];
      fort_forces(global_atom,0) += this_force.x();
      fort_forces(global_atom,1) += this_force.y();
      fort_forces(global_atom,2) += this_force.z();
    }
  }
  timers.forces.pause();

  timers.rmm.start();
  /* accumulate RMM results for this group */
  if(compute_rmm) {
    HostMatrix<scalar_type> rmm_output[inner_threads];

    const int ps = points.size();
    const int chunk = (ps + inner_threads - 1) / inner_threads;

    timers.rmm_calcs.start();
    #pragma omp parallel for num_threads(inner_threads)
    for(int i = 0; i < inner_threads; i++) {
      rmm_output[i].resize(group_m, group_m);
      rmm_output[i].zero();
      mkl_set_num_threads_local(1);
      for(int p = chunk*i;p < chunk*(i+1) && p < ps; p++) {
        HostMatrix<scalar_type>::blas_ssyr(LowerTriangle, factors_rmm[p],function_values, rmm_output[i], p);
      }
    }
    timers.rmm_calcs.pause();

    for(int step = 1; step < inner_threads; step *= 2) {
      const int incr = step * 2, stop = inner_threads - step;
      const int elems = group_m * group_m;
      #pragma omp parallel for num_threads(inner_threads / step)
      for(int k = 0; k < stop; k += incr) {
        scalar_type * out = rmm_output[k].data;
        const scalar_type * in = rmm_output[k+step].data;
        #pragma vector aligned always nontemporal
        for(int i = 0; i < elems; i++) {
            out[i] += in[i];
        }
      }
    }
  
    timers.rmm_update.start();
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
            rmm_global_output(big_index) += rmm_output[0](ii,jj);
          }
        }
      }
    }
    timers.rmm_update.pause();
  }

  timers.rmm.pause();
  energy+=localenergy;

#if CPU_RECOMPUTE
  /* clear functions */
  function_values.deallocate();
  gradient_values.deallocate();
  hessian_values.deallocate();
#endif
}

template class PointGroup<double>;
template class PointGroup<float>;

}
