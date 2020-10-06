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

template<class scalar_type>
void PointGroupGPU<scalar_type>::calc_W_mat(HostMatrix<double>& W_output_local){

   int device;
   cudaGetDevice(&device);
   current_device = device;

   /** Compute this group's functions **/
   compute_functions(false, true);
   uint group_m  = this->total_functions();
   uint n_points = this->number_of_points

   /* CUDA-related dimensions */
   const int block_height = divUp(group_m, 2*DENSITY_BLOCK_SIZE);
   dim3 threadBlock(DENSITY_BLOCK_SIZE,1,1);
   dim3 threadGrid(n_points,block_height,1);

   /* These are used for w accumulation in a given point of the grid */
   const dim3 threadGrid_accumulate(divUp(n_points,DENSITY_ACCUM_BLOCK_SIZE),1,1);
   const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE,1,1);

   /** Load points from group **/
   HostMatrix<scalar_type> point_weights_cpu(n_points, 1);
   uint i = 0;
   for (vector<Point>::const_iterator p = this->points.begin(); p != this->points.end(); ++p, ++i) {
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
                           (scalar_type) this->points[jpoint].atom_weights(iatom);
   }
   }
   becke_w_gpu = becke_w_cpu;

   CudaMatrix<scalar_type> cdft_factors_gpu;
   cdft_factors_gpu.resize(cdft_vars.regions * n_points);
   cdft_factors_gpu.zero();

   CudaMatrix<uint> cdft_atoms(cdft_vars.atoms);
   CudaMatrix<uint> cdft_natom(cdft_vars.natom);

   // Here we accumulate the w for each point.
   gpu_cdft_factors<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                             cdft_factors_gpu.data, cdft_natom.data, 
                                             cdft_atoms.data,  point_weights_gpu.data,
                                             becke_w_gpu.data, n_points,
                                             fortran_vars.atoms, cdft_vars.regions,
                                             cdft_vars.max_nat);
   
   cudaAssertNoError("compute_cdft_weights");
   

   // This part calculates Fi * w * Fj, with w being cdft_factors.
   // Only use enough blocks for lower triangle
   threadGrid = dim3(blocksPerRow*(blocksPerRow+1)/2);
   threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);

   CudaMatrix<scalar_type> w_output_gpu(COALESCED_DIMENSION(group_m), group_m);
   w_output_gpu.zero();

   // For calls with a single block (pretty common with cubes) don't bother doing 
   // the arithmetic to get block position in the matrix
   uint blocksPerRow = divUp(group_m, RMM_BLOCK_SIZE_XY);
   if (blocksPerRow > 1) {
      gpu_update_rmm<scalar_type,true><<<threadGrid, threadBlock>>>(cdft_factors_gpu.data, n_points,
                                                                    w_output_gpu.data, function_values.data,
                                                                    group_m);
   } else {
      gpu_update_rmm<scalar_type,false><<<threadGrid, threadBlock>>>(cdft_factors_gpu.data, n_points,
                                                                     w_output_gpu.data, function_values.data,
                                                                     group_m);
   }

   cudaAssertNoError("update_wmat");

   /*** Contribute this RMM to the total RMM ***/
   HostMatrix<scalar_type> rmm_output_cpu(rmm_output_gpu);
   this->add_rmm_output(rmm_output_cpu, rmm_output_local);

   /* clear functions */
   if(!(this->inGlobal)) {
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
}