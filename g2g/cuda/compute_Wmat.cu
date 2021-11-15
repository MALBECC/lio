/* -*- mode: c -*- */
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math_constants.h>
#include <string>
#include <vector>

#include "../common.h"
#include "../init.h"
#include "cuda_extra.h"
#include "../matrix.h"
#include "../partition.h"
#include "../scalar_vector_types.h"
#include "../global_memory_pool.h"

#include "gpu_variables.h"
#include "kernels/rmm.h"
#include "kernels/becke.h"

/*******************************************************************
** This subroutine is similar to solve_x, and it is used only     **
** for CDFT calculations. It only computes the W matrix required  **
** for Hab calculations in Mixed CDFT.                            **
** Since there is no difference between open and closed shells,   **
** at least at this level, we only need one subroutine.           **
*******************************************************************/
namespace G2G {

template <class scalar_type>
void PointGroupGPU<scalar_type>::calc_W_mat(HostMatrix<double>& W_output_local,
                                            CDFTVars& my_cdft_vars) {
  int device;
  cudaGetDevice(&device);
  current_device = device;

  /** Compute this group's functions **/
  compute_functions(false, false);
  uint group_m = this->total_functions();
  uint n_points = this->number_of_points;

  /* CUDA-related dimensions */
  const int block_height = divUp(group_m, 2 * DENSITY_BLOCK_SIZE);
  dim3 threadBlock(DENSITY_BLOCK_SIZE, 1, 1);
  dim3 threadGrid(n_points, block_height, 1);

  /* These are used for w accumulation in a given point of the grid */
  const dim3 threadGrid_accumulate(divUp(n_points, DENSITY_ACCUM_BLOCK_SIZE), 1,
                                   1);
  const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE, 1, 1);

  /** Load points from group **/
  HostMatrix<scalar_type> point_weights_cpu(n_points, 1);
  uint i = 0;
  for (vector<Point>::const_iterator p = this->points.begin();
       p != this->points.end(); ++p, ++i) {
    point_weights_cpu(i) = p->weight;
  }
  CudaMatrix<scalar_type> point_weights_gpu;
  point_weights_gpu = point_weights_cpu;

  // Becke Weights setup.
  CudaMatrix<scalar_type> becke_w_gpu;

  becke_w_gpu.resize(fortran_vars.atoms * n_points);
  HostMatrix<scalar_type> becke_w_cpu(fortran_vars.atoms * n_points);

  for (unsigned int jpoint = 0; jpoint < n_points; jpoint++) {
    for (unsigned int iatom = 0; iatom < fortran_vars.atoms; iatom++) {
      becke_w_cpu(jpoint * fortran_vars.atoms + iatom) =
          (scalar_type)this->points[jpoint].atom_weights(iatom);
    }
  }
  becke_w_gpu = becke_w_cpu;

  CudaMatrix<scalar_type> cdft_factors_gpu;
  cdft_factors_gpu.resize(my_cdft_vars.regions * n_points);
  cdft_factors_gpu.zero();

  CudaMatrix<uint> cdft_atoms(my_cdft_vars.atoms);
  CudaMatrix<uint> cdft_natom(my_cdft_vars.natom);

  // Here we accumulate the w for each point.
  gpu_cdft_factors<scalar_type>
      <<<threadGrid_accumulate, threadBlock_accumulate>>>(
          cdft_factors_gpu.data, cdft_natom.data, cdft_atoms.data,
          point_weights_gpu.data, becke_w_gpu.data, n_points,
          fortran_vars.atoms, my_cdft_vars.regions, my_cdft_vars.max_nat);

  cudaAssertNoError("compute_cdft_weights");

  // We now accumulate the constraints in a single W matrix.
  CudaMatrix<scalar_type> W_factors_gpu;
  W_factors_gpu.resize(n_points);
  W_factors_gpu.zero();

  CudaMatrix<scalar_type> cdft_Vc;
  if (my_cdft_vars.do_chrg) {
    HostMatrix<scalar_type> cdft_Vc_cpu(my_cdft_vars.regions);
    cdft_Vc.resize(my_cdft_vars.regions);
    for (unsigned int i = 0; i < my_cdft_vars.regions; i++) {
      cdft_Vc_cpu(i) = (scalar_type)my_cdft_vars.Vc(i);
    }
    cdft_Vc = cdft_Vc_cpu;

    gpu_cdft_factors_accum<scalar_type>
        <<<threadGrid_accumulate, threadBlock_accumulate>>>(
            cdft_factors_gpu.data, this->number_of_points, my_cdft_vars.regions,
            cdft_Vc.data, W_factors_gpu.data);
  }

  CudaMatrix<scalar_type> cdft_Vs;
  if (my_cdft_vars.do_spin) {
    HostMatrix<scalar_type> cdft_Vs_cpu(my_cdft_vars.regions);
    cdft_Vs.resize(my_cdft_vars.regions);
    for (unsigned int i = 0; i < my_cdft_vars.regions; i++) {
      cdft_Vs_cpu(i) = (scalar_type)-my_cdft_vars.Vs(i);
    }
    cdft_Vs = cdft_Vs_cpu;

    gpu_cdft_factors_accum<scalar_type>
        <<<threadGrid_accumulate, threadBlock_accumulate>>>(
            cdft_factors_gpu.data, this->number_of_points, my_cdft_vars.regions,
            cdft_Vs.data, W_factors_gpu.data);
  }

  // This part calculates Fi * w * Fj, with w being cdft_factors.
  // Only use enough blocks for lower triangle
  uint blocksPerRow = divUp(group_m, RMM_BLOCK_SIZE_XY);
  threadGrid = dim3(blocksPerRow * (blocksPerRow + 1) / 2);
  threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);

  CudaMatrix<scalar_type> w_output_gpu(COALESCED_DIMENSION(group_m), group_m);
  w_output_gpu.zero();

  // For calls with a single block (pretty common with cubes) don't bother doing
  // the arithmetic to get block position in the matrix
  if (blocksPerRow > 1) {
    gpu_update_rmm<scalar_type, true><<<threadGrid, threadBlock>>>(
        W_factors_gpu.data, n_points, w_output_gpu.data, function_values.data,
        group_m);
  } else {
    gpu_update_rmm<scalar_type, false><<<threadGrid, threadBlock>>>(
        W_factors_gpu.data, n_points, w_output_gpu.data, function_values.data,
        group_m);
  }

  cudaAssertNoError("update_wmat");

  /*** Contribute this RMM to the total RMM ***/
  HostMatrix<scalar_type> w_output_cpu(w_output_gpu);
  this->add_rmm_output(w_output_cpu, W_output_local);

  /* clear functions */
  if (!(this->inGlobal)) {
    function_values.deallocate();
  }
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupGPU<double>;
#else
template class PointGroup<float>;
template class PointGroupGPU<float>;
#endif
}  // namespace G2G