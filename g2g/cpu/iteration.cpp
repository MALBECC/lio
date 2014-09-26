/* -*- mode: c -*- */
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

void trmm(const CBLAS_ORDER Order, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, 
        const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const MKL_INT M, 
        const MKL_INT N, const float alpha, const float *A, const MKL_INT lda, 
        float *B, const MKL_INT ldb) {
    return cblas_strmm(Order,Side,Uplo,TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

void trmm(const CBLAS_ORDER Order, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, 
        const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const MKL_INT M, 
        const MKL_INT N, const double alpha, const double *A, const MKL_INT lda, 
        double *B, const MKL_INT ldb) {
    return cblas_dtrmm(Order,Side,Uplo,TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

template<class scalar_type>
void do_trmm(Timers & ts, const HostMatrix<scalar_type> & triagmat, const HostMatrix<scalar_type> & genmat, scalar_type * res) {
    ts.memcpy.start();
    genmat.copy_to_tmp(res);
    ts.memcpy.pause();
    ts.trmms.start();
    trmm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 
        genmat.height, genmat.width, 1.0, triagmat.asArray(), genmat.width, res, triagmat.height);
    ts.trmms.pause();
}

template<class scalar_type> void PointGroup<scalar_type>::solve(Timers& timers, 
    bool compute_rmm, bool lda, bool compute_forces, bool compute_energy, 
    double& energy, HostMatrix<double> & fort_forces, ThreadBufferPool<scalar_type> & pool, int pieces, HostMatrix<scalar_type> & rmm_global_output)
{
  uint group_m = total_functions();

  #if CPU_RECOMPUTE
  /** Compute functions **/
  timers.functions.start();
  compute_functions(compute_forces, !lda);
  timers.functions.pause();
  #endif

  double localenergy = 0.0;

  // prepare rmm_input for this group
  timers.density.start();
  HostMatrix<scalar_type> rmm_input(ALIGN(group_m), ALIGN(group_m));
  get_rmm_input(rmm_input);
  timers.density.pause();

  vector<vec_type3> forces(total_nucleii(), vec_type3(0.f,0.f,0.f)); 
  vector<std::vector<vec_type3> > forces_mat(
      points.size(), vector<vec_type3>(total_nucleii(), vec_type3(0.f,0.f,0.f)));
  vector<scalar_type> factors_rmm(points.size(),0);

  scalar_type * wv,* w3x,* w3y,* w3z,* ww1x,*ww1y,*ww1z,*ww2x,*ww2y,*ww2z;
  wv = w3x = w3y = w3z = ww1x = ww1y = ww1z = ww2x = ww2y = ww2z = NULL;
 
  if(!lda){
    timers.density.start();

    wv = pool.get_pool(); 
    w3x = pool.get_pool(); w3y = pool.get_pool(); w3z = pool.get_pool();
    ww1x = pool.get_pool(); ww1y = pool.get_pool(); ww1z = pool.get_pool();
    ww2x = pool.get_pool(); ww2y = pool.get_pool(); ww2z = pool.get_pool();
    
    int assign[10];
    memset(assign, 0, sizeof(assign));
    for(int i = 0; i < 10; i++) {
        assign[i] = i % pieces;
    }
    #pragma omp parallel for
    for(int i = 0; i < pieces; i++) {
        if (assign[0] == i) do_trmm(timers, rmm_input, function_values, wv);
        if (assign[1] == i) do_trmm(timers, rmm_input,  gX, w3x);
        if (assign[2] == i) do_trmm(timers, rmm_input,  gY, w3y);
        if (assign[3] == i) do_trmm(timers, rmm_input,  gZ, w3z);
        if (assign[4] == i) do_trmm(timers, rmm_input, hPX, ww1x);
        if (assign[5] == i) do_trmm(timers, rmm_input, hPY, ww1y);
        if (assign[6] == i) do_trmm(timers, rmm_input, hPZ, ww1z);
        if (assign[7] == i) do_trmm(timers, rmm_input, hIX, ww2x);
        if (assign[8] == i) do_trmm(timers, rmm_input, hIY, ww2y);
        if (assign[9] == i) do_trmm(timers, rmm_input, hIZ, ww2z);
    }

    timers.density.pause();
  }

  #pragma omp parallel for reduction(+:localenergy)
  for(int point = 0; point< points.size(); point++) {
    HostMatrix<vec_type3> dd;
    /** density **/
    scalar_type partial_density = 0;
    vec_type3 dxyz(0,0,0);
    vec_type3 dd1(0,0,0);
    vec_type3 dd2(0,0,0);

    timers.density.start();
    int skip = ALIGN(group_m);
    if (lda) {
      for (uint i = 0; i < skip; i++) {
        scalar_type w = 0.0;
        scalar_type Fi = function_values(i, point);
        for (uint j = i; j < skip; j++) {
          scalar_type Fj = function_values(j, point);
          w += rmm_input(j, i) * Fj;
        }
        partial_density += Fi * w;
      }
    } else {
      timers.density_calcs.start();
      #pragma ivdep
      #pragma vector aligned always
      for (int i = 0; i < group_m; i++) {
        int ai = point * skip + i;
        scalar_type w = wv[ai];
        vec_type3 w3(w3x[ai],w3y[ai],w3z[ai]);
        vec_type3 ww1(ww1x[ai],ww1y[ai],ww1z[ai]);
        vec_type3 ww2(ww2x[ai],ww2y[ai],ww2z[ai]);

        scalar_type Fi = function_values(i, point);
        vec_type3 Fgi(gX(i,point), gY(i, point), gZ(i, point));
        vec_type3 Fhi1(hPX(i,point), hPY(i, point), hPZ(i, point));
        vec_type3 Fhi2(hIX(i,point), hIY(i, point), hIZ(i, point));

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
    timers.density_calcs.pause();
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
      localenergy += (partial_density * points[point].weight) * (exc + corr);

    timers.density.pause();

    /** forces **/
    timers.forces.start();
    if (compute_forces) {
      scalar_type factor = points[point].weight * y2a;
      for (uint i = 0; i < total_nucleii(); i++) {
        forces_mat[point][i] = dd(i) * factor;
      }
    }
    timers.forces.pause();

    /** RMM **/
    timers.rmm.start();
    if (compute_rmm) {
      scalar_type factor = points[point].weight * y2a;
      factors_rmm[point] = factor;
    }
    timers.rmm.pause();
  } // end for

  timers.forces.start();
  /* accumulate forces for each point */
  if (compute_forces) {
    if(forces_mat.size() > 0) {
      #pragma omp parallel for
      for (int j = 0; j < forces_mat[0].size(); j++) {
        vec_type3 acum(0.f,0.f,0.f);
        for (int i = 0; i < forces_mat.size(); i++) {
          acum += forces_mat[i][j];
        }
        forces[j] = acum;
      }
    }
  }

  /* accumulate force results for this group */
  if (compute_forces) {
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
    HostMatrix<scalar_type> rmm_output_piece[pieces];
    for(int i = 0; i < pieces; i++) {
        rmm_output_piece[i].resize(ALIGN(group_m), ALIGN(group_m));
        rmm_output_piece[i].zero();
    }
    #pragma omp parallel for
    for(int piece = 0; piece < pieces; piece++){
        for(int i = 0;i < points.size();i++) {
          if(i % pieces != piece) continue;
          scalar_type factor = factors_rmm[i];
          HostMatrix<scalar_type>::blas_ssyr(LowerTriangle, factor, function_values, rmm_output_piece[piece], i);
        }
    }

    #pragma ivdep
    #pragma vector aligned always
    for(int p = 1; p < pieces; p++){
        for(int i = 0; i < group_m; i++) {
            for(int j = 0; j < group_m; j++) {
                rmm_output_piece[0](i,j) += rmm_output_piece[p](i,j);
            }
        }
    }
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
            rmm_global_output(big_index) += rmm_output_piece[0](ii, jj);
          }
        }
      }
    }
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
