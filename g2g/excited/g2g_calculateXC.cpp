#include <iostream>
#include <omp.h>

#include <stdio.h>
#include <string.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"
#include "../libxc/libxcproxy.h"

using namespace G2G;
extern Partition partition;

extern "C" void g2g_calculatexc_(double* Tmat, double* Fv) {
  partition.solve_lr(Tmat, Fv);
}

namespace G2G {

void Partition::solve_lr(double* T, double* F) {
  int M = fortran_vars.m;
  std::vector<HostMatrix<double> > Fock_output(G2G::cpu_threads +
                                               G2G::gpu_threads);

#pragma omp parallel for num_threads(cpu_threads + gpu_threads) schedule(static)
  for (uint i = 0; i < work.size(); i++) {
#if GPU_KERNELS
    bool gpu_thread = false;
    if (i >= cpu_threads) {
      gpu_thread = true;
      cudaSetDevice(i - cpu_threads);
    }
#endif

    Fock_output[i].resize(fortran_vars.rmm_output.width,
                          fortran_vars.rmm_output.height);
    // height es siempre 1: es un vector en realidad
    Fock_output[i].zero();
    for (uint j = 0; j < work[i].size(); j++) {
      Timer element;
      element.start_and_sync();
      int ind = work[i][j];
      if (ind >= cubes.size()) {
        spheres[ind - cubes.size()]->solve_closed_lr(T, Fock_output[i]);
      } else {
        cubes[ind]->solve_closed_lr(T, Fock_output[i]);
      }
#if GPU_KERNELS
      if (gpu_thread) cudaDeviceSynchronize();
#endif
      element.stop_and_sync();
    }
  }

  for (uint i = 0; i < work.size(); i++) {
    int index = 0;
    for (uint j = 0; j < M; j++) {
      F[j * M + j] += Fock_output[i](index, 0);
      index += 1;
      for (uint k = j + 1; k < M; k++) {
        F[j * M + k] += Fock_output[i](index, 0);
        F[k * M + j] = F[j * M + k];
        index += 1;
      }
    }
  }
  vector<HostMatrix<double> >().swap(Fock_output);
}

template <class scalar_type>
void PointGroupCPU<scalar_type>::solve_closed_lr(double* T,
                                                 HostMatrix<double>& Fock) {
  // aqui T no es la transicion density, sino es B * C**T
  const uint group_m = this->total_functions();
  const int npoints = this->points.size();
  bool lda = false;
  bool compute_forces = false;
  compute_functions(compute_forces, !lda);
  HostMatrix<scalar_type> rmm_input(group_m, group_m);

  int M = fortran_vars.m;
  get_rmm_input(rmm_input);

  double* smallFock = (double*)malloc(group_m * group_m * sizeof(double));
  memset(smallFock, 0.0f, group_m * group_m * sizeof(double));

  // FORMAMOS LA TRANSITION DENSITY REDUCIDA
  HostMatrix<double> tred(group_m, group_m);
  HostMatrix<double> Tbig(M * (M + 1) / 2);
  int index = 0;
  for (int row = 0; row < M; row++) {
    Tbig(index) = T[row * M + row];
    index += 1;
    for (int col = row + 1; col < M; col++) {
      Tbig(index) = T[row * M + col] + T[col * M + row];
      index += 1;
    }
  }
  get_tred_input(tred, Tbig);
  Tbig.deallocate();

// LIBXC INITIALIZATION
#define libxc_init_param                                              \
  fortran_vars.func_id, fortran_vars.func_coef, fortran_vars.nx_func, \
      fortran_vars.nc_func, fortran_vars.nsr_id, fortran_vars.screen, \
      XC_UNPOLARIZED
  LibxcProxy<scalar_type, 3> libxcProxy(libxc_init_param);
#undef libxc_init_param

  double* lrCoef = new double[3];
  double* tot_term = new double[group_m];

  for (int point = 0; point < npoints; point++) {
    scalar_type pd, fxc, tdx, tdy, tdz;
    pd = fxc = tdx = tdy = tdz = 0.0f;
    scalar_type red, redx, redy, redz;
    red = redx = redy = redz = 0.0f;
    const scalar_type* fv = function_values.row(point);
    const scalar_type* gxv = gX.row(point);
    const scalar_type* gyv = gY.row(point);
    const scalar_type* gzv = gZ.row(point);
    for (int i = 0; i < group_m; i++) {
      double z3xc, z3yc, z3zc, z;
      z3xc = z3yc = z3zc = z = 0.0f;
      for (int j = 0; j <= i; j++) {
        // Transition Density
        z += fv[j] * tred(i, j);
        z3xc += gxv[j] * tred(i, j);
        z3yc += gyv[j] * tred(i, j);
        z3zc += gzv[j] * tred(i, j);
      }
      const double Fi = fv[i];
      const double gx = gxv[i], gy = gyv[i], gz = gzv[i];
      // Transition Density
      red += Fi * z;
      redx += gx * z + z3xc * Fi;
      redy += gy * z + z3yc * Fi;
      redz += gz * z + z3zc * Fi;
    }

    // Ground State Density
    pd = rho_values(0, point);
    tdx = rho_values(1, point);
    tdy = rho_values(2, point);
    tdz = rho_values(3, point);

    double sigma = tdx * tdx + tdy * tdy + tdz * tdz;
    double cruz = redx * tdx + redy * tdy + redz * tdz;
    cruz *= 0.50f;
    tdx *= 0.5f;
    tdy *= 0.5f;
    tdz *= 0.5f;

    libxcProxy.coefLR(&pd, &sigma, red, cruz, lrCoef);
    const scalar_type wp = this->points[point].weight;
    double term1, term2, term3, term4, tot_term_ii, result;

    for (int i = 0; i < group_m; i++) {
      term1 = lrCoef[0] * 0.5f * fv[i] + lrCoef[1] * tdx * gxv[i];
      term2 = lrCoef[1] * tdy * gyv[i] + lrCoef[1] * tdz * gzv[i];
      term3 = lrCoef[2] * redx * gxv[i] + lrCoef[2] * redy * gyv[i];
      term4 = lrCoef[2] * redz * gzv[i];
      tot_term[i] = (term1 + term2 + term3 + term4) * wp;
      tot_term_ii = tot_term[i];
      for (int j = 0; j <= i; j++) {
        result = fv[i] * tot_term[j] + tot_term_ii * fv[j];
        smallFock[i * group_m + j] += result;
      }
    }

  }  // END points loop

  const int indexes = this->rmm_bigs.size();
  for (int i = 0; i < indexes; i++) {
    int bi = this->rmm_bigs[i], row = this->rmm_rows[i],
        col = this->rmm_cols[i];
    Fock(bi) += smallFock[col * group_m + row];
  }
  // Free Memory
  free(smallFock);
  smallFock = NULL;
  delete[] tot_term;
  tot_term = NULL;
  delete[] lrCoef;
  lrCoef = NULL;
}
template <class scalar_type>
void PointGroupCPU<scalar_type>::get_tred_input(
    HostMatrix<scalar_type>& tred_input, HostMatrix<double>& source) const {
  tred_input.zero();
  const int indexes = this->rmm_bigs.size();
  for (int i = 0; i < indexes; i++) {
    int ii = this->rmm_rows[i], jj = this->rmm_cols[i], bi = this->rmm_bigs[i];
    tred_input(ii, jj) = tred_input(jj, ii) = (scalar_type)source(bi);
  }
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}  // namespace G2G
