/* -*- mode: c -*- */
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math_constants.h>
#include <string>
#include <vector>
#include <stdio.h>
#include "../../init.h"
#include "../../common.h"
#include "../../partition.h"
#include "../../timer.h"

#if USE_LIBXC
#include "../../libxc/libxc_accumulate_point.h"
#endif

#include "accumulate_values.h"

using namespace std;

namespace G2G {
#if FULL_DOUBLE
texture<int2, 2, cudaReadModeElementType> tred_gpu_3rd_tex;
#else
texture<float, 2, cudaReadModeElementType> tred_gpu_3rd_tex;
#endif

#include "../../cuda/kernels/transpose.h"
#include "obtain_fock_cuda.h"
#include "obtain_terms.h"

// comienzan las nuevas
#include "ES_compute_3rd_partial.h"

template <class scalar_type>
void PointGroupGPU<scalar_type>::solve_3rd_der(double* T,
                                               HostMatrix<double>& Fock,
                                               int DER) {
  uint group_m = this->total_functions();
  bool lda = false;
  bool compute_forces = false;

  compute_functions(compute_forces, !lda);

  CudaMatrix<scalar_type> point_weights_gpu;
  HostMatrix<scalar_type> point_weights_cpu(this->number_of_points, 1);
  uint i = 0;
  for (vector<Point>::const_iterator p = this->points.begin();
       p != this->points.end(); ++p, ++i) {
    point_weights_cpu(i) = p->weight;
  }
  point_weights_gpu = point_weights_cpu;
  point_weights_cpu.deallocate();

  // Variables to kernels: gpu_compute
  dim3 threadBlock, threadGrid;
  const int block_height = divUp(group_m, 2 * DENSITY_BLOCK_SIZE);
  threadBlock =
      dim3(DENSITY_BLOCK_SIZE, 1,
           1);  // Hay que asegurarse que la cantidad de funciones este en rango
  threadGrid = dim3(this->number_of_points, block_height, 1);

  // Partial Transition Density GPU
  CudaMatrix<scalar_type> partial_tred_gpu;
  CudaMatrix<vec_type<scalar_type, 4>> tredxyz_gpu;
  partial_tred_gpu.resize(COALESCED_DIMENSION(this->number_of_points),
                          block_height);
  tredxyz_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);

  // Accumulate Values
  CudaMatrix<scalar_type> tred_accum_gpu;
  CudaMatrix<vec_type<scalar_type, 4>> tredxyz_accum_gpu;
  tred_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
  tredxyz_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));

  // Variables to kernels: accumulate density
  const dim3 threadGrid_accumulate(
      divUp(this->number_of_points, DENSITY_ACCUM_BLOCK_SIZE), 1, 1);
  const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE, 1, 1);

  // Transpose functions and gradients
#define BLOCK_DIM 16
  int transposed_width = COALESCED_DIMENSION(this->number_of_points);
  dim3 transpose_grid(transposed_width / BLOCK_DIM, divUp((group_m), BLOCK_DIM),
                      1);
  dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);

  CudaMatrix<scalar_type> function_values_transposed;
  CudaMatrix<vec_type<scalar_type, 4>> gradient_values_transposed;
  function_values_transposed.resize(
      group_m, COALESCED_DIMENSION(this->number_of_points));
  gradient_values_transposed.resize(
      group_m, COALESCED_DIMENSION(this->number_of_points));

  transpose<<<transpose_grid, transpose_threads>>>(
      function_values_transposed.data, function_values.data,
      COALESCED_DIMENSION(this->number_of_points), group_m);

  transpose<<<transpose_grid, transpose_threads>>>(
      gradient_values_transposed.data, gradient_values.data,
      COALESCED_DIMENSION(this->number_of_points), group_m);

  // FORM reduce Transition density
  HostMatrix<scalar_type> tred_cpu(COALESCED_DIMENSION(group_m),
                                   group_m + DENSITY_BLOCK_SIZE);
  int M = fortran_vars.m;

  HostMatrix<double> Tbig(M * (M + 1) / 2);
  int index = 0;
  int row, col;
  for (row = 0; row < M; row++) {
    Tbig(index) = T[row * M + row];
    index += 1;
    for (col = row + 1; col < M; col++) {
      Tbig(index) = T[row * M + col] + T[col * M + row];
      index += 1;
    }
  }
  get_tred_input(tred_cpu, Tbig);
  Tbig.deallocate();

  // ponemos ceros fuera de group_m
  for (uint i = 0; i < (group_m + DENSITY_BLOCK_SIZE); i++) {
    for (uint j = 0; j < COALESCED_DIMENSION(group_m); j++) {
      if ((i >= group_m) || (j >= group_m) || (j > i)) {
        tred_cpu.data[COALESCED_DIMENSION(group_m) * i + j] = 0.0f;
      }
    }
  }

  // Form Bind Textures
  cudaArray* cuArraytred;
  cudaMallocArray(&cuArraytred, &tred_gpu_3rd_tex.channelDesc, tred_cpu.width,
                  tred_cpu.height);
  cudaMemcpyToArray(cuArraytred, 0, 0, tred_cpu.data,
                    sizeof(scalar_type) * tred_cpu.width * tred_cpu.height,
                    cudaMemcpyHostToDevice);
  cudaBindTextureToArray(tred_gpu_3rd_tex, cuArraytred);
  tred_cpu.deallocate();

// CALCULATE PARTIAL DENSITIES
#define compden_parameter                                           \
  this->number_of_points, function_values_transposed.data, group_m, \
      gradient_values_transposed.data, partial_tred_gpu.data, tredxyz_gpu.data
  ES_compute_3rd_partial<scalar_type, true, true, false>
      <<<threadGrid, threadBlock>>>(compden_parameter);

// ACCUMULATE DENSITIES
#define accumden_parameter                                     \
  this->number_of_points, block_height, partial_tred_gpu.data, \
      tredxyz_gpu.data, tred_accum_gpu.data, tredxyz_accum_gpu.data
  accumulate_values<scalar_type, true, true, false>
      <<<threadGrid_accumulate, threadBlock_accumulate>>>(accumden_parameter);

#undef compute_parameter
#undef accumden_parameter

  // LIBXC INITIALIZATION
  fortran_vars.fexc = fortran_vars.func_coef[0];
#define libxc_init_param                                              \
  fortran_vars.func_id, fortran_vars.func_coef, fortran_vars.nx_func, \
      fortran_vars.nc_func, fortran_vars.nsr_id, fortran_vars.screen, \
      XC_UNPOLARIZED
  LibxcProxy_cuda<scalar_type, 4> libxcProxy_cuda(libxc_init_param);
#undef libxc_init_param

  CudaMatrix<scalar_type> lrCoef_gpu;
  lrCoef_gpu.resize(COALESCED_DIMENSION(this->number_of_points));

  // DEFINE OUTPUTS
  CudaMatrix<vec_type<scalar_type, 4>> Txyz;
  CudaMatrix<vec_type<scalar_type, 4>> Dxyz;
  Txyz.resize(COALESCED_DIMENSION(this->number_of_points));
  Dxyz.resize(COALESCED_DIMENSION(this->number_of_points));

  if (DER == 2) {
    libxc_gpu_coefLR<scalar_type, true, true, false>(
        &libxcProxy_cuda, this->number_of_points, rmm_accum_gpu.data,
        tred_accum_gpu.data, dxyz_accum_gpu.data, tredxyz_accum_gpu.data,
        // Outputs
        Dxyz.data, Txyz.data, lrCoef_gpu.data);
  } else {
    libxc_gpu_coefZv<scalar_type, true, true, false>(
        &libxcProxy_cuda, this->number_of_points, rmm_accum_gpu.data,
        tred_accum_gpu.data, dxyz_accum_gpu.data, tredxyz_accum_gpu.data,
        // Outputs
        Dxyz.data, Txyz.data, lrCoef_gpu.data);
  }

  // CALCULATE TERMS
  CudaMatrix<scalar_type> terms_lr;
  terms_lr.resize(group_m, COALESCED_DIMENSION(this->number_of_points));
  threadGrid = dim3(this->number_of_points);
  threadBlock = dim3(group_m);

#define compute_parameters                                              \
  this->number_of_points, group_m, point_weights_gpu.data,              \
      function_values_transposed.data, gradient_values_transposed.data, \
      lrCoef_gpu.data, Dxyz.data, Txyz.data, terms_lr.data

  gpu_obtain_term<scalar_type, false, true, false>
      <<<threadGrid, threadBlock>>>(compute_parameters);
#undef compute_parameters

  // CALCULATE FOCK
  const int block_size = divUp(group_m, RMM_BLOCK_SIZE_XY);
  const int MM = block_size * (block_size + 1) / 2;
  threadGrid = dim3(MM, 1, 1);
  threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);

  CudaMatrix<scalar_type> smallFock_gpu(group_m, group_m);
  smallFock_gpu.zero();
#define compute_parameters                                 \
  this->number_of_points, group_m, point_weights_gpu.data, \
      function_values_transposed.data, terms_lr.data, smallFock_gpu.data

  gpu_obtain_fock<scalar_type, false, true, false>
      <<<threadGrid, threadBlock>>>(compute_parameters);
#undef compute_parameters

  // Obtain the global fock
  HostMatrix<scalar_type> smallFock(smallFock_gpu);
  smallFock_gpu.deallocate();
  this->add_rmm_output(smallFock, Fock);

  // Free Memory
  smallFock.deallocate();
  cudaUnbindTexture(tred_gpu_3rd_tex);
  cudaFreeArray(cuArraytred);
  Txyz.deallocate();
  Dxyz.deallocate();
  partial_tred_gpu.deallocate();
  tredxyz_gpu.deallocate();
  tred_accum_gpu.deallocate();
  tredxyz_accum_gpu.deallocate();
  lrCoef_gpu.deallocate();
  point_weights_gpu.deallocate();
  function_values_transposed.deallocate();
  gradient_values_transposed.deallocate();
  terms_lr.deallocate();
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupGPU<double>;
#else
template class PointGroup<float>;
template class PointGroupGPU<float>;
#endif
}  // namespace G2G
