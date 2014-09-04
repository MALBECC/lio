/* -*- mode: c -*- */
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cstring>
#include "mkl.h"
#include "../common.h"
#include "../init.h"
#include "../cuda_includes.h"
#include "../matrix.h"
#include "../timer.h"
#include "../partition.h"

#include "cpu/pot.h"

using std::cout;
using std::endl;
using std::list;
using std::vector;

namespace G2G { 

float * do_trmm(const HostMatrix<float> & triagmat, const HostMatrix<float> & genmat) {
    float * res = (float *) mkl_malloc(genmat.width * genmat.height * sizeof(float),64);
    memcpy(res, genmat.asArray(), genmat.width * genmat.height * sizeof(float));
    cblas_strmm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 
        genmat.height, genmat.width, 1.0, triagmat.asArray(), triagmat.height, res, triagmat.height);
    return res;
}

double * do_trmm(const HostMatrix<double> & triagmat, const HostMatrix<double> & genmat) {
    return NULL;
}

template< int compo, int skip, int start >
float * do_trmm_proyect(const HostMatrix<float> & triagmat, const HostMatrix< vec_type< float, 3> > & genmat) {
    int width = genmat.width / skip;
    float * res = (float *) mkl_malloc(width * genmat.height * sizeof(float),64);
    for(int row = 0, pos = 0; row < genmat.height; row++){
        for(int col = start; col < genmat.width; col += skip){
            vec_type<float, 3> e = genmat(col, row);
            res[pos++] = (compo == 0) ? e.x() : (compo == 1) ? e.y() : e.z();
        }
    }
    cblas_strmm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 
        genmat.height, width, 1.0, triagmat.asArray(), triagmat.height, res, triagmat.height);
    return res;
}

template< int compo, int skip, int start >
double * do_trmm_proyect(const HostMatrix<double> & triagmat, const HostMatrix< vec_type< double, 3> > & genmat) {
    return NULL;
}

template<class scalar_type> void PointGroup<scalar_type>::solve(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, 
    bool compute_energy, double& energy, double* fort_forces_ptr)
{
  HostMatrix<scalar_type> rmm_output;
  uint group_m = total_functions();
  if (compute_rmm) { rmm_output.resize(group_m, group_m); rmm_output.zero(); }

  #if CPU_RECOMPUTE
  /** Compute functions **/
  timers.functions.start();
  compute_functions(compute_forces, !lda);
  timers.functions.pause();
  #endif
  double localenergy = 0.0;
  // prepare rmm_input for this group
  timers.density.start();
  HostMatrix<scalar_type> rmm_input(group_m, group_m);
  get_rmm_input(rmm_input);
  timers.density.pause();

  HostMatrix<vec_type3> forces(total_nucleii(), 1); forces.zero();
  vector<std::vector<vec_type3> > forces_mat(
      points.size(), vector<vec_type3>(total_nucleii(), vec_type3(0.f,0.f,0.f)));
  vector<scalar_type> factors_rmm(points.size(),0);
  /******** each point *******/
  vector<Point> _points(points.begin(),points.end());

  scalar_type * wv,* w3x,* w3y,* w3z,* ww1x,*ww1y,*ww1z,*ww2x,*ww2y,*ww2z;
  wv = w3x = w3y = w3z = ww1x = ww1y = ww1z = ww2x = ww2y = ww2z = NULL;
 
  if(!lda){
    wv   = do_trmm(rmm_input, function_values);
    w3x  = do_trmm_proyect<0,1,0>(rmm_input, gradient_values);
    w3y  = do_trmm_proyect<1,1,0>(rmm_input, gradient_values);
    w3z  = do_trmm_proyect<2,1,0>(rmm_input, gradient_values);
    ww1x = do_trmm_proyect<0,2,0>(rmm_input, hessian_values);
    ww1y = do_trmm_proyect<1,2,0>(rmm_input, hessian_values);
    ww1z = do_trmm_proyect<2,2,0>(rmm_input, hessian_values);
    ww2x = do_trmm_proyect<0,2,1>(rmm_input, hessian_values);
    ww2y = do_trmm_proyect<1,2,1>(rmm_input, hessian_values);
    ww2z = do_trmm_proyect<2,2,1>(rmm_input, hessian_values);
  }

#pragma omp parallel for reduction(+:localenergy)
  for(int point = 0; point<_points.size(); point++)
  {
    HostMatrix<vec_type3> dd;
    /** density **/
    scalar_type partial_density = 0;
    vec_type3 dxyz(0,0,0);
    vec_type3 dd1(0,0,0);
    vec_type3 dd2(0,0,0);

    timers.density.start();
    if (lda) {
      for (uint i = 0; i < group_m; i++) {
        scalar_type w = 0.0;
        scalar_type Fi = function_values(i, point);
        for (uint j = i; j < group_m; j++) {
          scalar_type Fj = function_values(j, point);
          w += rmm_input(j, i) * Fj;
        }
        partial_density += Fi * w;
      }
    }
    else {
      for (int i = 0; i < group_m; i++) {
        int ai = point * group_m + i;
        scalar_type w = wv[ai];
        vec_type3 w3(w3x[ai],w3y[ai],w3z[ai]);
        vec_type3 ww1(ww1x[ai],ww1y[ai],ww1z[ai]);
        vec_type3 ww2(ww2x[ai],ww2y[ai],ww2z[ai]);

        scalar_type Fi = function_values(i, point);
        vec_type3 Fgi(gradient_values(i, point));
        vec_type3 Fhi1(hessian_values(2 * (i + 0) + 0, point));
        vec_type3 Fhi2(hessian_values(2 * (i + 0) + 1, point));

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
      localenergy += (partial_density * _points[point].weight) * (exc + corr);

    timers.density.pause();

    /** forces **/
    timers.forces.start();
    if (compute_forces) {
      scalar_type factor = _points[point].weight * y2a;
      for (uint i = 0; i < total_nucleii(); i++) {
        forces_mat[point][i] = dd(i) * factor;
      }
    }
    timers.forces.pause();

    /** RMM **/
    timers.rmm.start();
    if (compute_rmm) {
      scalar_type factor = _points[point].weight * y2a;
      factors_rmm[point] = factor;
    }
    timers.rmm.pause();
  } // end for

  if (compute_rmm) {
    for(int i=0; i<_points.size(); i++) {
      scalar_type factor = factors_rmm[i];
      HostMatrix<scalar_type>::blas_ssyr(LowerTriangle, factor, function_values, rmm_output, i);
    }
  }

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
        forces(j) = acum;
      }
    }
  }

  /* accumulate force results for this group */
  if (compute_forces) {
    FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, fortran_vars.max_atoms); // TODO: mover esto a init.cpp
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
  energy+=localenergy;

  mkl_free(wv);
  mkl_free(w3x);
  mkl_free(w3y);
  mkl_free(w3z);
  mkl_free(ww1x);
  mkl_free(ww1y);
  mkl_free(ww1z);
  mkl_free(ww2x);
  mkl_free(ww2y);
  mkl_free(ww2z);

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
