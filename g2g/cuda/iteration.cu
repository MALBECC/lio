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
#include "../timer.h"
#include "../partition.h"
#include "../scalar_vector_types.h"
#include "../global_memory_pool.h"

#include "../pointxc/calc_ggaCS.h"
#include "../pointxc/calc_ggaOS.h"
#include "../pointxc/calc_ldaCS.h"

#if USE_LIBXC
#include "../libxc/libxc_accumulate_point.h"
#endif

namespace G2G {
/*#if FULL_DOUBLE
texture<int2, 2, cudaReadModeElementType> rmm_input_gpu_tex;
texture<int2, 2, cudaReadModeElementType> rmm_input_gpu_tex2;
#else
texture<float, 2, cudaReadModeElementType> rmm_input_gpu_tex;
texture<float, 2, cudaReadModeElementType> rmm_input_gpu_tex2;
#endif*/
/** KERNELS **/
// Including CUDA kernel header files
#include "gpu_variables.h"
#include "kernels/accumulate_point.h"
#include "kernels/energy.h"
#include "kernels/energy_open.h"
#include "kernels/energy_derivs.h"
#include "kernels/rmm.h"
#include "kernels/weight.h"
#include "kernels/functions.h"
#include "kernels/force.h"
#include "kernels/transpose.h"
#include "kernels/becke.h"

using std::cout;
using std::vector;
using std::endl;

//extern "C" void g2g_timer_sum_start_(const char* timer_name, unsigned int length_arg);
//extern "C" void g2g_timer_sum_stop_(const char* timer_name, unsigned int length_arg);
//extern "C" void g2g_timer_sum_pause_(const char* timer_name, unsigned int length_arg);

/**
 * Sets up GPU global variables by transferring parameters from host to device
 * memory across all available GPUs.
 */
void gpu_set_variables(void) {
  int previous_device; cudaGetDevice(&previous_device);
  int gpu_devices = cudaGetGPUCount();
  
  // Copy essential parameters to constant memory on every available GPU
  for(int i = 0; i < gpu_devices; i++) {
    if(cudaSetDevice(i) != cudaSuccess)
      std::cout << "Error: can't set the device " << i << std::endl;
    cudaMemcpyToSymbol(gpu_normalization_factor, &fortran_vars.normalization_factor, sizeof(fortran_vars.normalization_factor), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(gpu_atoms, &fortran_vars.atoms, sizeof(fortran_vars.atoms), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(gpu_Iexch, &fortran_vars.iexch, sizeof(fortran_vars.iexch), 0, cudaMemcpyHostToDevice);
  }
  
  // Restore original CUDA device
  cudaSetDevice(previous_device);
  cudaAssertNoError("set_gpu_variables");
}

/**
 * Copies atomic positions from host to device memory for all available GPUs.
 * 
 * @param m Host matrix containing atomic position data
 */
template<class T> void gpu_set_atom_positions(const HostMatrix<T>& m) {
  int previous_device; cudaGetDevice(&previous_device);
  int gpu_devices = cudaGetGPUCount();
  
  // Copy positions to each available GPU
  for(int i = 0; i < gpu_devices; i++) {
    if(cudaSetDevice(i) != cudaSuccess)
      std::cout << "Error: can't set the device " << i << std::endl;
    cudaMemcpyToSymbol(gpu_atom_positions, m.data, m.bytes(), 0, cudaMemcpyHostToDevice);
  }
  
  // Restore original CUDA device
  cudaSetDevice(previous_device);
}

// Explicit template instantiation depending on precision setting
#if FULL_DOUBLE
template void gpu_set_atom_positions<double3>(const HostMatrix<double3>& m);
#else
template void gpu_set_atom_positions<float3>(const HostMatrix<float3>& m);
#endif

/**
 * Main entry point for DFT calculations on a point group.
 * This is a wrapper that routes calculations to either closed-shell or open-shell
 * implementations depending on the 'open' parameter.
 */
template<class scalar_type>
void PointGroupGPU<scalar_type>::solve(
    Timers& timers, bool compute_rmm, bool lda, bool compute_forces,
    bool compute_energy, double& energy,double& energy_i, double& energy_c,
    double& energy_c1, double& energy_c2,  HostMatrix<double>& fort_forces_ms,
    int inner_threads, HostMatrix<double>& rmm_output_local, bool open){
/*
  if ( open ) {
      solve_opened( timers, compute_rmm, lda, compute_forces, compute_energy,
                    energy, energy_i, energy_c, energy_c1, energy_c2,
                    fort_forces_ms );
  }
  else {
      solve_closed( timers, compute_rmm, lda, compute_forces, compute_energy,
                    energy, fort_forces_ms, inner_threads, rmm_output_local );
  }
*/
//  counter_iter++;                                                            // For Debug FF
//  std::cout << "Grupo " << counter_iter << " Energia : " << energy << " \n"; // For Debug FF
}

/**
 * Implementation of closed-shell DFT calculations for a point group.
 * 
 * This is the main computational workhorse for closed-shell systems, handling:
 * - Function evaluation
 * - Density calculation
 * - Energy calculation
 * - Force calculation
 * - RMM updates
 * 
 * @param timers Timing objects to measure performance
 * @param compute_rmm Whether to update reduced density matrices
 * @param lda Whether to use LDA (vs. GGA) functionals 
 * @param compute_forces Whether to compute atomic forces
 * @param compute_energy Whether to compute energy
 * @param energy Output energy value (accumulates result)
 * @param fort_forces_ms Output forces matrix
 * @param inner_threads Thread count for inner loops
 * @param rmm_output_local Output matrix for RMM updates
 * @param becke_dens Output densities for Becke partitioning
 * @param my_cdft_vars Variables for constrained DFT
 */
template<class scalar_type>
void PointGroupGPU<scalar_type>::solve_closed(
    Timers& timers,
    bool compute_rmm, bool lda, bool compute_forces, bool compute_energy,
    double& energy,    HostMatrix<double>& fort_forces_ms,
    int inner_threads, HostMatrix<double>& rmm_output_local,
    HostMatrix<double>& becke_dens, CDFTVars& my_cdft_vars){

  // Get current CUDA device ID
  int device;
  cudaGetDevice(&device);
  current_device = device;

  /*** Computo sobre cada cubo ****/
  CudaMatrix<scalar_type> point_weights_gpu;

  /** Compute this group's functions **/
  // Start timer and calculate basis functions at all points
  timers.functions.start_and_sync();
  compute_functions(compute_forces, !lda);
  timers.functions.pause_and_sync();

  // Get total number of basis functions in this group
  uint group_m = this->total_functions();

  // Begin density calculation
  timers.density.start_and_sync();
  
  /** Load points from group **/
  // Transfer point weights from host to device
  HostMatrix<scalar_type> point_weights_cpu(this->number_of_points, 1);

  uint i = 0;
  for (vector<Point>::const_iterator p = this->points.begin(); p != this->points.end(); ++p, ++i) {
    point_weights_cpu(i) = p->weight;
  }

  point_weights_gpu = point_weights_cpu;

  // Setup CUDA thread organization
  dim3 threadBlock, threadGrid;
  /* compute density/factors */

  // Calculate the number of vertical thread blocks needed
  const int block_height= divUp(group_m, 2*DENSITY_BLOCK_SIZE);

  // Configure thread organization for density calculation
  threadBlock = dim3(DENSITY_BLOCK_SIZE,1,1); // Hay que asegurarse que la cantidad de funciones este en rango
  threadGrid = dim3(this->number_of_points,block_height,1);

  // Allocate device matrices for density calculation
  CudaMatrix<scalar_type> partial_densities_gpu;
  CudaMatrix< vec_type<scalar_type,4> > dxyz_gpu;
  CudaMatrix< vec_type<scalar_type,4> > dd1_gpu;
  CudaMatrix< vec_type<scalar_type,4> > dd2_gpu;

  // Resize matrices according to problem dimensions
  partial_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
  dxyz_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height);
  dd1_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height );
  dd2_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height );

  // Allocate additional matrices for LibXC if enabled
#if USE_LIBXC
  CudaMatrix<scalar_type> accumulated_densities_gpu;
  CudaMatrix< vec_type<scalar_type,4> > dxyz_accum_gpu;
  CudaMatrix< vec_type<scalar_type,4> > dd1_accum_gpu;
  CudaMatrix< vec_type<scalar_type,4> > dd2_accum_gpu;

  accumulated_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
  dxyz_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
  dd1_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
  dd2_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
#endif

  // Thread configuration for accumulation kernels
  const dim3 threadGrid_accumulate(divUp(this->number_of_points,DENSITY_ACCUM_BLOCK_SIZE),1,1);
  const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE,1,1);

  // Allocate factors matrix if needed for forces or RMM
  CudaMatrix<scalar_type> factors_gpu;
  if (compute_rmm || compute_forces) {
    factors_gpu.resize(this->number_of_points);
    factors_gpu.zero();
  }

  // Setup for matrix transposition
  int transposed_width = COALESCED_DIMENSION(this->number_of_points);
  #define BLOCK_DIM 16
  dim3 transpose_grid(transposed_width / BLOCK_DIM, divUp((group_m),BLOCK_DIM), 1);
  dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);

  // Transpose function values for better memory coalescing 
  CudaMatrix<scalar_type> function_values_transposed;
  function_values_transposed.resize(group_m, COALESCED_DIMENSION(this->number_of_points));
  transpose<<<transpose_grid, transpose_threads>>> (function_values_transposed.data,
      function_values.data, COALESCED_DIMENSION(this->number_of_points), group_m);

  // Transpose gradient values if needed for forces or GGA
  CudaMatrix<vec_type<scalar_type,4> > gradient_values_transposed;
  if (fortran_vars.do_forces || fortran_vars.gga) {
    gradient_values_transposed.resize( group_m,COALESCED_DIMENSION(this->number_of_points));
    transpose<<<transpose_grid, transpose_threads>>> (gradient_values_transposed.data,
        gradient_values.data, COALESCED_DIMENSION(this->number_of_points), group_m );
  }
  
  // Prepare density matrix input
  HostMatrix<scalar_type> rmm_input_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
  get_rmm_input(rmm_input_cpu); //Achica la matriz densidad a la version reducida del grupo

  // Zero out upper triangle and out-of-bounds elements
  for (uint i=0; i<(group_m+DENSITY_BLOCK_SIZE); i++)
  {
    for(uint j=0; j<COALESCED_DIMENSION(group_m); j++)
    {
      if((i>=group_m) || (j>=group_m) || (j > i))
      {
        rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
      }
    }
  }

  // Transfer density matrix to device
  CudaMatrix<scalar_type> rmm_input_gpu;
  rmm_input_gpu=rmm_input_cpu;
  
  /*
   **********************************************************************
   * Pasando RDM (rmm) a texturas
   **********************************************************************
   */
/*
  cudaArray* cuArray;
  cudaMallocArray(&cuArray, &rmm_input_gpu_tex.channelDesc, rmm_input_cpu.width, rmm_input_cpu.height);
  cudaMemcpyToArray(cuArray, 0, 0, rmm_input_cpu.data, sizeof(scalar_type)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);
  cudaBindTextureToArray(rmm_input_gpu_tex, cuArray);

  rmm_input_gpu_tex.normalized = false;
*/

  // Initialize LibXC proxy if LibXC is enabled
//#if USE_LIBXC
//  if (fortran_vars.use_libxc) fortran_vars.fexc = fortran_vars.func_coef[0];
//#define libxc_init_param \
//  fortran_vars.func_id, fortran_vars.func_coef, fortran_vars.nx_func, \
//  fortran_vars.nc_func, fortran_vars.nsr_id, fortran_vars.screen, \
//  XC_UNPOLARIZED
//  LibxcProxy_cuda<scalar_type,4> libxcProxy_cuda(libxc_init_param);
//#undef libxc_init_param
//#endif

  // Setup for Becke partitioning and CDFT calculations
  CudaMatrix<scalar_type> becke_w_gpu;
  CudaMatrix<scalar_type> cdft_factors_gpu;
  if (((my_cdft_vars.do_chrg || my_cdft_vars.do_spin) && compute_rmm) ||
       (fortran_vars.becke && compute_energy)) {
    becke_w_gpu.resize(fortran_vars.atoms * this->number_of_points);
    HostMatrix<scalar_type> becke_w_cpu(fortran_vars.atoms * this->number_of_points);

    // Transfer Becke weights to GPU
    for (unsigned int jpoint = 0; jpoint < this->number_of_points; jpoint++) {
      for (unsigned int iatom = 0; iatom < fortran_vars.atoms; iatom++) {
        becke_w_cpu(jpoint * fortran_vars.atoms + iatom) =
                          (scalar_type) this->points[jpoint].atom_weights(iatom);
     }
    }
    becke_w_gpu = becke_w_cpu;
  }

  // Energy calculation path
  if (compute_energy) {
    CudaMatrix<scalar_type> energy_gpu(this->number_of_points);

    // Define macro parameters for kernel calls
#define compute_parameters \
    energy_gpu.data, factors_gpu.data, point_weights_gpu.data, this->number_of_points, function_values_transposed.data, \
    gradient_values_transposed.data, hessian_values_transposed.data, group_m, partial_densities_gpu.data, dxyz_gpu.data, \
    dd1_gpu.data,dd2_gpu.data, rmm_input_gpu.data

#define accumulate_parameters \
    energy_gpu.data, factors_gpu.data, point_weights_gpu.data, this->number_of_points, block_height, \
    partial_densities_gpu.data, dxyz_gpu.data, dd1_gpu.data, dd2_gpu.data, fortran_vars.fexc

    // Branch based on whether we need to compute forces/RMM and LDA vs. GGA
    if (compute_forces || compute_rmm) {
      if (lda) {
        // LDA with forces/RMM
        gpu_compute_density<scalar_type, true, true, true><<<threadGrid, threadBlock>>>(compute_parameters);
        gpu_accumulate_point<scalar_type, true, true, true><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
      } else {
        // GGA with forces/RMM
        gpu_compute_density<scalar_type, true, true, false><<<threadGrid, threadBlock>>>(compute_parameters);
#if USE_LIBXC
	      if (fortran_vars.use_libxc) {
          fortran_vars.fexc = fortran_vars.func_coef[0];
#define libxc_init_param \
          fortran_vars.func_id, fortran_vars.func_coef, fortran_vars.nx_func, \
          fortran_vars.nc_func, fortran_vars.nsr_id, fortran_vars.screen, \
          XC_UNPOLARIZED
          LibxcProxy_cuda<scalar_type,4> libxcProxy_cuda(libxc_init_param);
#undef libxc_init_param
	        // Accumulate the data for libxc
	        gpu_accumulate_point_for_libxc<scalar_type, true, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
		          point_weights_gpu.data, this->number_of_points, block_height,
		          partial_densities_gpu.data, dxyz_gpu.data, dd1_gpu.data, dd2_gpu.data,
		          accumulated_densities_gpu.data, dxyz_accum_gpu.data, dd1_accum_gpu.data, dd2_accum_gpu.data);

	        // Compute exc_corr and y2a with libxc GPU version.
	        libxc_exchange_correlation_gpu<scalar_type, true, true, false> (&libxcProxy_cuda,
	          	energy_gpu.data, factors_gpu.data, this->number_of_points,
	        	  accumulated_densities_gpu.data, dxyz_accum_gpu.data, dd1_accum_gpu.data, dd2_accum_gpu.data);

	        // Merge the results.
	        gpu_accumulate_energy_and_forces_from_libxc<scalar_type, true, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
	          	energy_gpu.data, factors_gpu.data, point_weights_gpu.data, this->number_of_points, accumulated_densities_gpu.data);
	      } else {
          gpu_accumulate_point<scalar_type, true, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
	      }
#else
	      gpu_accumulate_point<scalar_type, true, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
#endif
      }

      // Handle constrained DFT calculation if needed
      if (my_cdft_vars.do_chrg) {
        cdft_factors_gpu.resize(my_cdft_vars.regions * this->number_of_points);
        cdft_factors_gpu.zero();

        CudaMatrix<uint> cdft_atoms(my_cdft_vars.atoms);
        CudaMatrix<uint> cdft_natom(my_cdft_vars.natom);

        gpu_cdft_factors<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                            cdft_factors_gpu.data, cdft_natom.data, 
                                            cdft_atoms.data,  point_weights_gpu.data,
                                            becke_w_gpu.data, this->number_of_points,
                                            fortran_vars.atoms, my_cdft_vars.regions, my_cdft_vars.max_nat);
      }
    } else {
      // Energy-only calculation (no forces or RMM)
      if (lda) {
        // LDA, energy only
        gpu_compute_density<scalar_type, true, false, true><<<threadGrid, threadBlock>>>(compute_parameters);
        gpu_accumulate_point<scalar_type, true, false, true><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
      } else {
        // GGA, energy only
        gpu_compute_density<scalar_type, true, false, false><<<threadGrid, threadBlock>>>(compute_parameters);
#if USE_LIBXC
        if (fortran_vars.use_libxc) {
          fortran_vars.fexc = fortran_vars.func_coef[0];
#define libxc_init_param \
          fortran_vars.func_id, fortran_vars.func_coef, fortran_vars.nx_func, \
          fortran_vars.nc_func, fortran_vars.nsr_id, fortran_vars.screen, \
          XC_UNPOLARIZED
          LibxcProxy_cuda<scalar_type,4> libxcProxy_cuda(libxc_init_param);
#undef libxc_init_param
      	  // Accumulate the data.
	        gpu_accumulate_point_for_libxc<scalar_type, true, false, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
              point_weights_gpu.data, this->number_of_points, block_height, partial_densities_gpu.data,
              dxyz_gpu.data, dd1_gpu.data, dd2_gpu.data, accumulated_densities_gpu.data,
              dxyz_accum_gpu.data, dd1_accum_gpu.data, dd2_accum_gpu.data);
              
	        // Compute exc_corr and y2a with libxc GPU version.
	        libxc_exchange_correlation_gpu<scalar_type, true, true, false> (&libxcProxy_cuda,
            energy_gpu.data, factors_gpu.data, this->number_of_points, accumulated_densities_gpu.data,
            dxyz_accum_gpu.data, dd1_accum_gpu.data, dd2_accum_gpu.data);
	        // Merge the results.
	        gpu_accumulate_energy_and_forces_from_libxc<scalar_type, true, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
	          energy_gpu.data, factors_gpu.data, point_weights_gpu.data, this->number_of_points, accumulated_densities_gpu.data);
	      } else {
          gpu_accumulate_point<scalar_type, true, false, false><<<threadGrid_accumulate, threadBlock_accumulate>>>(accumulate_parameters);
        }
#else
        gpu_accumulate_point<scalar_type, true, false, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
#endif
      }
    }
    cudaAssertNoError("compute_density");

    // Transfer energy results back to host and accumulate 
    HostMatrix<scalar_type> energy_cpu(energy_gpu);
    for (uint i = 0; i < this->number_of_points; i++) {
      energy += energy_cpu(i);
    }

    // Handle Becke partitioning if enabled
    if (fortran_vars.becke) {
      CudaMatrix<scalar_type> becke_dens_gpu(fortran_vars.atoms * this->number_of_points);
      becke_dens_gpu.zero();
      gpu_compute_becke_cs<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(becke_dens_gpu.data,
                                          partial_densities_gpu.data, point_weights_gpu.data, becke_w_gpu.data,
                                          this->number_of_points, fortran_vars.atoms, block_height);

      HostMatrix<scalar_type> becke_dens_cpu(becke_dens_gpu);
      for (unsigned int jpoint = 0; jpoint < this->number_of_points; jpoint++) {
        for (unsigned int iatom = 0; iatom < fortran_vars.atoms; iatom++) {
          becke_dens(iatom) += (double) becke_dens_cpu(jpoint * fortran_vars.atoms + iatom);
        }
      }
    }
  } else {
    // No energy calculation, but still compute density for forces or RMM
#undef compute_parameters
#undef accumulate_parameters

#define compute_parameters \
    NULL,factors_gpu.data,point_weights_gpu.data,this->number_of_points,function_values_transposed.data,gradient_values_transposed.data,hessian_values_transposed.data,group_m,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data,rmm_input_gpu.data
#define accumulate_parameters \
    NULL,factors_gpu.data,point_weights_gpu.data,this->number_of_points,block_height,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data, fortran_vars.fexc
    if (lda)
    {
        gpu_compute_density<scalar_type, false, true, true><<<threadGrid, threadBlock>>>(compute_parameters);
        gpu_accumulate_point<scalar_type, false, true, true><<<threadGrid_accumulate, threadBlock_accumulate>>>(accumulate_parameters);
    }
    else
    {
        gpu_compute_density<scalar_type, false, true, false><<<threadGrid, threadBlock>>>(compute_parameters);
#if USE_LIBXC
  if (fortran_vars.use_libxc) {
    fortran_vars.fexc = fortran_vars.func_coef[0];
#define libxc_init_param \
    fortran_vars.func_id, fortran_vars.func_coef, fortran_vars.nx_func, \
    fortran_vars.nc_func, fortran_vars.nsr_id, fortran_vars.screen, \
    XC_UNPOLARIZED
    LibxcProxy_cuda<scalar_type,4> libxcProxy_cuda(libxc_init_param);
#undef libxc_init_param
	  // Accumulate the data.
	  gpu_accumulate_point_for_libxc<scalar_type, false, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (point_weights_gpu.data,
            this->number_of_points, block_height,
	    partial_densities_gpu.data, dxyz_gpu.data, dd1_gpu.data, dd2_gpu.data,
	    accumulated_densities_gpu.data, dxyz_accum_gpu.data, dd1_accum_gpu.data, dd2_accum_gpu.data);

	  // Compute exc_corr and y2a with libxc GPU version.
	  libxc_exchange_correlation_gpu<scalar_type, false, true, false> (&libxcProxy_cuda,
	    NULL, factors_gpu.data, this->number_of_points,
	    accumulated_densities_gpu.data, dxyz_accum_gpu.data, dd1_accum_gpu.data, dd2_accum_gpu.data);

	  // Merge the results.
	  gpu_accumulate_energy_and_forces_from_libxc<scalar_type, false, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
	    NULL,factors_gpu.data, point_weights_gpu.data, this->number_of_points, accumulated_densities_gpu.data);
	} else {
    	  gpu_accumulate_point<scalar_type, false, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>>(accumulate_parameters);
	}
#else
        gpu_accumulate_point<scalar_type, false, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>>(accumulate_parameters);
#endif
    }
    
    // Constrained DFT - charge constraints
    if (my_cdft_vars.do_chrg) {
      cdft_factors_gpu.resize(my_cdft_vars.regions * this->number_of_points);
      cdft_factors_gpu.zero();

      CudaMatrix<uint> cdft_atoms(my_cdft_vars.atoms);
      CudaMatrix<uint> cdft_natom(my_cdft_vars.natom);    

      gpu_cdft_factors<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                          cdft_factors_gpu.data, cdft_natom.data, 
                                          cdft_atoms.data,  point_weights_gpu.data,
                                          becke_w_gpu.data, this->number_of_points,
                                          fortran_vars.atoms, my_cdft_vars.regions, my_cdft_vars.max_nat);
    }    
    cudaAssertNoError("compute_density");
  }
#undef compute_parameters
#undef accumulate_parameters

  timers.density.pause_and_sync();
  
  /* compute forces */
  if (compute_forces) {
    // Restore full density matrix - we need the full matrix for forces
    // (upper triangle elements were zeroed earlier)
    for (uint i=0; i<(group_m); i++) {
      for(uint j=0; j<(group_m); j++) {
        if((i>=group_m) || (j>=group_m) || (j > i))
        {
          rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*j+i] ;
        }
      }
    }

    // Create a density matrix with upper triangle restored
    CudaMatrix<scalar_type> rmm_input_gpu_sinceros;
    rmm_input_gpu_sinceros=rmm_input_cpu; // la versión en la GPU ya existe, puede que haya problemas con eso
    
    timers.density_derivs.start_and_sync();
//    cudaMemcpyToArray(cuArray, 0, 0,rmm_input_cpu.data,
//      sizeof(scalar_type)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);

    // Setup thread configuration for density derivatives calculation
    timers.density_derivs.start_and_sync();
    dim3 threads = dim3(this->number_of_points);
    threadBlock = dim3(DENSITY_DERIV_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);

    // Allocate memory for density derivatives
    CudaMatrix<vec_type4> dd_gpu(COALESCED_DIMENSION(this->number_of_points), this->total_nucleii()); dd_gpu.zero();
    CudaMatrixUInt nuc_gpu(this->func2local_nuc);  // TODO: esto en realidad se podria guardar una sola vez durante su construccion

    // Calculate density derivatives for force computation
    gpu_compute_density_derivs<<<threadGrid, threadBlock>>>(
        function_values.data, gradient_values.data, nuc_gpu.data, dd_gpu.data, this->number_of_points, group_m, this->total_nucleii(),rmm_input_gpu_sinceros.data);
    cudaAssertNoError("density_derivs");
    timers.density_derivs.pause_and_sync();

    // Begin force calculation
    timers.forces.start_and_sync();
    CudaMatrix<vec_type4> forces_gpu(this->total_nucleii());
    forces_gpu.zero();

    // Configure threads for force computation
    threads = dim3(this->total_nucleii());
    threadBlock = dim3(FORCE_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);
    
    // Compute forces from density derivatives
    gpu_compute_forces<<<threadGrid, threadBlock>>>(
        this->number_of_points, factors_gpu.data, dd_gpu.data, forces_gpu.data, this->total_nucleii());
    cudaAssertNoError("forces");

    // Transfer force results back to host and accumulate
    HostMatrix<vec_type4> forces_cpu(forces_gpu);

    for (uint i = 0; i < this->total_nucleii(); ++i) {
      vec_type4 atom_force = forces_cpu(i);
      uint global_nuc = this->local2global_nuc[i];
      fort_forces_ms(global_nuc, 0) += atom_force.x;
      fort_forces_ms(global_nuc, 1) += atom_force.y;
      fort_forces_ms(global_nuc, 2) += atom_force.z;

    }
    timers.forces.pause_and_sync();
  }

  // Begin RMM computation
  timers.rmm.start_and_sync();
  /* compute RMM */
  if (compute_rmm) {
    // Configure thread organization for RMM update
    threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);
    uint blocksPerRow = divUp(group_m, RMM_BLOCK_SIZE_XY);
    
    // Only use enough blocks for lower triangle
    threadGrid = dim3(blocksPerRow*(blocksPerRow+1)/2);

    // Allocate output matrix for RMM updates
    CudaMatrix<scalar_type> rmm_output_gpu(COALESCED_DIMENSION(group_m), group_m);
    rmm_output_gpu.zero();


    // Adds CDFT terms to RMM factors.
    CudaMatrix<scalar_type> cdft_Vc;
    if (my_cdft_vars.do_chrg) {
      // Transfer CDFT constraint potentials to GPU
      HostMatrix<scalar_type> cdft_Vc_cpu(my_cdft_vars.regions);

      cdft_Vc.resize(my_cdft_vars.regions);
      for (unsigned int i = 0; i < my_cdft_vars.regions; i++) {
        cdft_Vc_cpu(i) = (scalar_type) my_cdft_vars.Vc(i);
      }
      cdft_Vc = cdft_Vc_cpu;
   
      // Apply CDFT potentials to integration factors
      gpu_cdft_factors_accum<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                                  cdft_factors_gpu.data, this->number_of_points,
                                                  my_cdft_vars.regions, cdft_Vc.data, factors_gpu.data);
    }


    // Update RMM with optimized kernel selection based on problem size
    // For calls with a single block (pretty common with cubes) don't bother doing the arithmetic to get block position in the matrix
    if (blocksPerRow > 1) {
        gpu_update_rmm<scalar_type,true><<<threadGrid, threadBlock>>>(factors_gpu.data, this->number_of_points,
                                                                      rmm_output_gpu.data, function_values.data,
                                                                      group_m);
    } else {
        gpu_update_rmm<scalar_type,false><<<threadGrid, threadBlock>>>(factors_gpu.data, this->number_of_points,
                                                                       rmm_output_gpu.data, function_values.data,
                                                                       group_m);
    }

    cudaAssertNoError("update_rmm");

    /*** Contribute this RMM to the total RMM ***/
    // Transfer results to host and add to global RMM
    HostMatrix<scalar_type> rmm_output_cpu(rmm_output_gpu);
    this->add_rmm_output(rmm_output_cpu, rmm_output_local);
  }
  timers.rmm.pause_and_sync();

  /* clear functions */
  // Free GPU memory if not needed for later use
  if(!(this->inGlobal)) {
    function_values.deallocate();
    gradient_values.deallocate();
    hessian_values_transposed.deallocate();
  }
  //Deshago el bind de textura de rmm
//  cudaUnbindTexture(rmm_input_gpu_tex); //Enroque el Unbind con el Free, asi parece mas logico. Nano
//  cudaFreeArray(cuArray);
}

//======================
// OPENSHELL
//======================

/**
 * Implementation of open-shell DFT calculations for a point group.
 * 
 * This variant handles systems with unpaired electrons (open-shell),
 * managing separate computations for alpha and beta spin components.
 * 
 * @param timers Timing objects to measure performance
 * @param compute_rmm Whether to update reduced density matrices
 * @param lda Whether to use LDA (vs. GGA) functionals 
 * @param compute_forces Whether to compute atomic forces
 * @param compute_energy Whether to compute energy
 * @param energy Output total energy value
 * @param energy_i Output initial energy value
 * @param energy_c Output correlation energy value
 * @param energy_c1 Output additional correlation energy component
 * @param energy_c2 Output additional correlation energy component
 * @param fort_forces_ms Output forces matrix
 * @param rmm_output_local_a Output matrix for alpha-spin RMM updates
 * @param rmm_output_local_b Output matrix for beta-spin RMM updates
 * @param becke_dens Output densities for Becke partitioning
 * @param becke_spin Output spin densities for Becke partitioning
 * @param my_cdft_vars Variables for constrained DFT
 */
template<class scalar_type>
void PointGroupGPU<scalar_type>::solve_opened(
    Timers& timers, bool compute_rmm, bool lda, bool compute_forces,
    bool compute_energy, double& energy, double& energy_i,
    double& energy_c, double& energy_c1, double& energy_c2,
    HostMatrix<double>& fort_forces_ms, HostMatrix<double>& rmm_output_local_a,
    HostMatrix<double>& rmm_output_local_b, HostMatrix<double>& becke_dens,
    HostMatrix<double>& becke_spin, CDFTVars& my_cdft_vars){

  // Get current CUDA device
  int device;
  cudaGetDevice(&device);
  current_device = device;

  /*** Computo sobre cada cubo ****/
  CudaMatrix<scalar_type> point_weights_gpu;

  /** Compute this group's functions **/
  // Start timer and calculate basis functions
  timers.functions.start_and_sync();
  compute_functions(compute_forces, !lda);
  timers.functions.pause_and_sync();

  // Get total number of basis functions
  uint group_m = this->total_functions();

  // Begin density calculation
  timers.density.start_and_sync();
  
  /** Load points from group **/
  // Transfer point weights from host to device
  HostMatrix<scalar_type> point_weights_cpu(this->number_of_points, 1);

  uint i = 0;
  for (vector<Point>::const_iterator p = this->points.begin(); p != this->points.end(); ++p, ++i) {
    point_weights_cpu(i) = p->weight;
  }
  point_weights_gpu = point_weights_cpu;

  // Setup CUDA thread organization
  dim3 threadBlock, threadGrid;
  const int block_height= divUp(group_m,2*DENSITY_BLOCK_SIZE);

  // This makes sure the amount of functions fits within range.
  threadBlock = dim3(DENSITY_BLOCK_SIZE,1,1);
  threadGrid = dim3(this->number_of_points,block_height,1);

  // Allocate matrices for alpha and beta components
  CudaMatrix<scalar_type> factors_a_gpu;
  CudaMatrix<scalar_type> factors_b_gpu;

  // Gradients (dxyz) and Hessians (dd1,dd2) for alpha/beta.
  CudaMatrix<scalar_type> partial_densities_a_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dxyz_a_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd1_a_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd2_a_gpu;

  CudaMatrix<scalar_type> partial_densities_b_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dxyz_b_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd1_b_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd2_b_gpu;

  // Matrix transpose is needed for better coalescence in density.
  CudaMatrix<scalar_type> function_values_transposed;
  CudaMatrix<vec_type<scalar_type,4> > gradient_values_transposed;

  int transposed_width = COALESCED_DIMENSION(this->number_of_points);

  // Allocate transposed matrices
  function_values_transposed.resize(group_m, COALESCED_DIMENSION(this->number_of_points));
  if (fortran_vars.do_forces || fortran_vars.gga)
      gradient_values_transposed.resize( group_m,COALESCED_DIMENSION(this->number_of_points));

  // Configure thread layout for transposition
  #define BLOCK_DIM 16
  dim3 transpose_grid(transposed_width / BLOCK_DIM, divUp((group_m),BLOCK_DIM));
  dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);

  // Transpose function and gradient values for better memory coalescing
  transpose<<<transpose_grid, transpose_threads>>> (function_values_transposed.data, function_values.data,  COALESCED_DIMENSION(this->number_of_points),group_m   );
  if (fortran_vars.do_forces || fortran_vars.gga)
      transpose<<<transpose_grid, transpose_threads>>> (gradient_values_transposed.data, gradient_values.data, COALESCED_DIMENSION(this->number_of_points), group_m );

  // Allocate matrices for partial results
  partial_densities_a_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
  dxyz_a_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height);
  dd1_a_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height );
  dd2_a_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height );

  partial_densities_b_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
  dxyz_b_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height);
  dd1_b_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height );
  dd2_b_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height );

  // Thread configuration for accumulation
  const dim3 threadGrid_accumulate(divUp(this->number_of_points,DENSITY_ACCUM_BLOCK_SIZE),1,1);
  const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE,1,1);

  // Allocate factors matrices if needed for forces or RMM
  if (compute_rmm || compute_forces) {
    factors_a_gpu.resize(this->number_of_points);
    factors_b_gpu.resize(this->number_of_points);
    factors_a_gpu.zero();
    factors_b_gpu.zero();
  }

  // Prepare density matrices for alpha and beta components
  HostMatrix<scalar_type> rmm_input_a_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
  HostMatrix<scalar_type> rmm_input_b_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
   //Reduces density matrixes (Up,Down) to the reduced group version
  get_rmm_input(rmm_input_a_cpu, rmm_input_b_cpu);

  // Zero out upper triangle and out-of-bounds elements
  for (uint i=0; i<(group_m+DENSITY_BLOCK_SIZE); i++) {
    for(uint j=0; j<COALESCED_DIMENSION(group_m); j++) {
      if((i>=group_m) || (j>=group_m) || (j > i)) {
        rmm_input_a_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
        rmm_input_b_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
      }
    }
  }
  
  // Transfer density matrices to GPU
  CudaMatrix<scalar_type> rmm_input_gpu1;
  CudaMatrix<scalar_type> rmm_input_gpu2;

  rmm_input_gpu1=rmm_input_a_cpu;
  rmm_input_gpu2=rmm_input_b_cpu;

  /*
  **********************************************************************
  * Pasando RDM (rmm) a texturas/
  **********************************************************************
  */
/*
  cudaArray* cuArray1;
  cudaArray* cuArray2;
  cudaMallocArray(&cuArray1, &rmm_input_gpu_tex.channelDesc, rmm_input_a_cpu.width,rmm_input_a_cpu.height);
  cudaMallocArray(&cuArray2, &rmm_input_gpu_tex2.channelDesc, rmm_input_b_cpu.width,rmm_input_b_cpu.height);
  cudaMemcpyToArray(cuArray1, 0, 0,rmm_input_a_cpu.data,sizeof(scalar_type)*rmm_input_a_cpu.width*rmm_input_a_cpu.height, cudaMemcpyHostToDevice);
  cudaMemcpyToArray(cuArray2, 0, 0,rmm_input_b_cpu.data,sizeof(scalar_type)*rmm_input_b_cpu.width*rmm_input_b_cpu.height, cudaMemcpyHostToDevice);
  cudaBindTextureToArray(rmm_input_gpu_tex, cuArray1);
  cudaBindTextureToArray(rmm_input_gpu_tex2, cuArray2);

  rmm_input_gpu_tex.normalized = false;
  rmm_input_gpu_tex2.normalized = false;
*/
  // For CDFT and becke partitioning.
  CudaMatrix<scalar_type> becke_w_gpu;
  CudaMatrix<scalar_type> cdft_factors_gpu;
  if (((my_cdft_vars.do_chrg || my_cdft_vars.do_spin) && compute_rmm) ||
      (fortran_vars.becke && compute_energy)) {
    becke_w_gpu.resize(fortran_vars.atoms * this->number_of_points);
    HostMatrix<scalar_type> becke_w_cpu(fortran_vars.atoms * this->number_of_points);

    // Transfer Becke weights to GPU
    for (unsigned int jpoint = 0; jpoint < this->number_of_points; jpoint++) {
      for (unsigned int iatom = 0; iatom < fortran_vars.atoms; iatom++) {
        becke_w_cpu(jpoint * fortran_vars.atoms + iatom) =
                          (scalar_type) this->points[jpoint].atom_weights(iatom);
    }
    }
    becke_w_gpu = becke_w_cpu;
  }

  // Energy calculation path
  if (compute_energy) {
    CudaMatrix<scalar_type> energy_gpu(this->number_of_points);

    // Run density and accumulation kernels for open-shell system
    if (compute_forces || compute_rmm) {
      // Compute density with forces or RMM updates
      gpu_compute_density_opened<scalar_type, true, true, false><<<threadGrid, threadBlock>>>(
             point_weights_gpu.data,this->number_of_points, function_values_transposed.data,
             gradient_values_transposed.data,hessian_values_transposed.data, group_m,
             partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
             partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data, rmm_input_gpu1.data, rmm_input_gpu2.data);
      
      // Accumulate energy and factors
      gpu_accumulate_point_open<scalar_type, true, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
             energy_gpu.data,
             factors_a_gpu.data, factors_b_gpu.data, point_weights_gpu.data,this->number_of_points,block_height,
             partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
             partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data, fortran_vars.fexc);
      
      // Handle CDFT if enabled
      if (my_cdft_vars.do_chrg || my_cdft_vars.do_spin ) {        
        cdft_factors_gpu.resize(my_cdft_vars.regions * this->number_of_points);
        cdft_factors_gpu.zero();

        CudaMatrix<uint> cdft_atoms(my_cdft_vars.atoms);
        CudaMatrix<uint> cdft_natom(my_cdft_vars.natom);    

        gpu_cdft_factors<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                            cdft_factors_gpu.data, cdft_natom.data, 
                                            cdft_atoms.data,  point_weights_gpu.data,
                                            becke_w_gpu.data, this->number_of_points,
                                            fortran_vars.atoms, my_cdft_vars.regions, my_cdft_vars.max_nat);
      }
    } else {
      // Compute density for energy only (no forces or RMM)
      gpu_compute_density_opened<scalar_type, true, false, false><<<threadGrid, threadBlock>>>(
             point_weights_gpu.data,this->number_of_points, function_values_transposed.data,
             gradient_values_transposed.data,hessian_values_transposed.data, group_m,
             partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
             partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data, rmm_input_gpu1.data, rmm_input_gpu2.data);
      
      // Accumulate energy only
      gpu_accumulate_point_open<scalar_type, true, false, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
             energy_gpu.data, factors_a_gpu.data, factors_b_gpu.data, point_weights_gpu.data,
             this->number_of_points, block_height,
             partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
             partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data, fortran_vars.fexc);
    }
    cudaAssertNoError("compute_density");

    // Transfer energy results back to host and accumulate
    HostMatrix<scalar_type> energy_cpu(energy_gpu);

    for (uint i = 0; i < this->number_of_points; i++) {
      energy    += energy_cpu(i);
    }

     // Handle Becke partitioning if enabled
     if (fortran_vars.becke) {   
      CudaMatrix<scalar_type> becke_dens_gpu(fortran_vars.atoms * this->number_of_points);
      CudaMatrix<scalar_type> becke_spin_gpu(fortran_vars.atoms * this->number_of_points);
      becke_dens_gpu.zero();
      becke_spin_gpu.zero();
      
      // Compute Becke partitioning for open-shell system
      gpu_compute_becke_os<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>
                          (becke_dens_gpu.data, becke_spin_gpu.data, partial_densities_a_gpu.data,
                           partial_densities_b_gpu.data, point_weights_gpu.data, becke_w_gpu.data,
                           this->number_of_points, fortran_vars.atoms, block_height);

      // Transfer results to host and accumulate
      HostMatrix<scalar_type> becke_dens_cpu(becke_dens_gpu);
      HostMatrix<scalar_type> becke_spin_cpu(becke_spin_gpu);
      for (unsigned int jpoint = 0; jpoint < this->number_of_points; jpoint++) {
        for (unsigned int iatom = 0; iatom < fortran_vars.atoms; iatom++) {
          becke_dens(iatom) += (double) becke_dens_cpu(jpoint * fortran_vars.atoms + iatom);
          becke_spin(iatom) += (double) becke_spin_cpu(jpoint * fortran_vars.atoms + iatom);
        }
      }
    }
  } else {
    // No energy calculation, but compute density for forces or RMM
    gpu_compute_density_opened<scalar_type, false, true, false><<<threadGrid, threadBlock>>>(
           point_weights_gpu.data, this->number_of_points, function_values_transposed.data,
           gradient_values_transposed.data,hessian_values_transposed.data, group_m,
           partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
           partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data, rmm_input_gpu1.data, rmm_input_gpu2.data);
    
    // Accumulate forces and RMM factors
    gpu_accumulate_point_open<scalar_type, false, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
           NULL,
           factors_a_gpu.data, factors_b_gpu.data, point_weights_gpu.data,this->number_of_points,block_height,
           partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
           partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data, fortran_vars.fexc);

    // Handle CDFT if enabled
    if (my_cdft_vars.do_chrg || my_cdft_vars.do_spin ) {           
      cdft_factors_gpu.resize(my_cdft_vars.regions * this->number_of_points);
      cdft_factors_gpu.zero();

      CudaMatrix<uint> cdft_atoms(my_cdft_vars.atoms);
      CudaMatrix<uint> cdft_natom(my_cdft_vars.natom);    

      gpu_cdft_factors<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                          cdft_factors_gpu.data, cdft_natom.data, 
                                          cdft_atoms.data,  point_weights_gpu.data,
                                          becke_w_gpu.data, this->number_of_points,
                                          fortran_vars.atoms, my_cdft_vars.regions, my_cdft_vars.max_nat);
    }
    cudaAssertNoError("compute_density");
  }

  timers.density.pause_and_sync();


  /* compute forces */
  if (compute_forces) {

    // Restore full density matrices - we need the full matrices for forces
    // (upper triangle elements were zeroed earlier)
    for (uint i=0; i<(group_m); i++) {
    for (uint j=0; j<(group_m); j++) {
      if((i>=group_m) || (j>=group_m) || (j > i)){
        rmm_input_a_cpu.data[COALESCED_DIMENSION(group_m)*i+j] =
                        rmm_input_a_cpu.data[COALESCED_DIMENSION(group_m)*j+i] ;
        rmm_input_b_cpu.data[COALESCED_DIMENSION(group_m)*i+j] =
                        rmm_input_b_cpu.data[COALESCED_DIMENSION(group_m)*j+i] ;
      }
    }
    }

 //   cudaMemcpyToArray(cuArray1, 0, 0,rmm_input_a_cpu.data,sizeof(scalar_type)*rmm_input_a_cpu.width*rmm_input_a_cpu.height, cudaMemcpyHostToDevice);
 //   cudaMemcpyToArray(cuArray2, 0, 0,rmm_input_b_cpu.data,sizeof(scalar_type)*rmm_input_b_cpu.width*rmm_input_b_cpu.height, cudaMemcpyHostToDevice);
  rmm_input_gpu1=rmm_input_a_cpu;
  rmm_input_gpu2=rmm_input_b_cpu;


    dim3 threads;
    timers.density_derivs.start_and_sync();
    
    // Configure threads for density derivatives
    threads = dim3(this->number_of_points);
    threadBlock = dim3(DENSITY_DERIV_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);

    // Allocate matrices for density derivatives
    CudaMatrix<vec_type4> dd_gpu_a(COALESCED_DIMENSION(this->number_of_points), this->total_nucleii());
    CudaMatrix<vec_type4> dd_gpu_b(COALESCED_DIMENSION(this->number_of_points), this->total_nucleii());
    dd_gpu_a.zero();
    dd_gpu_b.zero();
    CudaMatrixUInt nuc_gpu(this->func2local_nuc);

    // Compute density derivatives for both alpha and beta components
    gpu_compute_density_derivs_open<<<threadGrid, threadBlock>>>(function_values.data, gradient_values.data, nuc_gpu.data, dd_gpu_a.data, dd_gpu_b.data, this->number_of_points, group_m, this->total_nucleii(),rmm_input_gpu1.data,rmm_input_gpu2.data);

    cudaAssertNoError("density_derivs");
    timers.density_derivs.pause_and_sync();

    // Begin force calculation
    timers.forces.start_and_sync();
    CudaMatrix<vec_type4> forces_gpu_a(this->total_nucleii());
    CudaMatrix<vec_type4> forces_gpu_b(this->total_nucleii());
    forces_gpu_a.zero();
    forces_gpu_b.zero();

    // Configure threads for force computation
    threads = dim3(this->total_nucleii());
    threadBlock = dim3(FORCE_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);
    
    // Compute forces for both alpha and beta components
    gpu_compute_forces<<<threadGrid, threadBlock>>>(this->number_of_points, factors_a_gpu.data, dd_gpu_a.data, forces_gpu_a.data, this->total_nucleii());
    gpu_compute_forces<<<threadGrid, threadBlock>>>(this->number_of_points, factors_b_gpu.data, dd_gpu_b.data, forces_gpu_b.data, this->total_nucleii());

    cudaAssertNoError("forces");

    // Transfer results to host and accumulate forces
    HostMatrix<vec_type4> forces_cpu_a(forces_gpu_a);
    HostMatrix<vec_type4> forces_cpu_b(forces_gpu_b);

    for (uint i = 0; i < this->total_nucleii(); ++i) {
      vec_type4 atom_force_a = forces_cpu_a(i);
      vec_type4 atom_force_b = forces_cpu_b(i);
      uint global_nuc = this->local2global_nuc[i];

      // Sum contributions from alpha and beta electrons for each coordinate
      fort_forces_ms(global_nuc, 0) += atom_force_a.x + atom_force_b.x;
      fort_forces_ms(global_nuc, 1) += atom_force_a.y + atom_force_b.y;
      fort_forces_ms(global_nuc, 2) += atom_force_a.z + atom_force_b.z;
    }

    timers.forces.pause_and_sync();
  }

  /* compute RMM */
  timers.rmm.start_and_sync();
  if (compute_rmm) {
    // Configure thread organization for RMM updates
    threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);
    uint blocksPerRow = divUp(group_m, RMM_BLOCK_SIZE_XY);

    // Only use enough blocks for lower triangle (optimization)
    threadGrid = dim3(blocksPerRow*(blocksPerRow+1)/2);
    
    // Allocate output matrices for alpha and beta RMM updates
    CudaMatrix<scalar_type> rmm_output_a_gpu(COALESCED_DIMENSION(group_m), group_m);
    CudaMatrix<scalar_type> rmm_output_b_gpu(COALESCED_DIMENSION(group_m), group_m);

    rmm_output_a_gpu.zero();
    rmm_output_b_gpu.zero();

    // CDFT processing for charge and spin constraints
    CudaMatrix<scalar_type> cdft_Vc;
    CudaMatrix<scalar_type> cdft_Vs_a;
    CudaMatrix<scalar_type> cdft_Vs_b;
    
    // Handle charge constraints if enabled
    if (my_cdft_vars.do_chrg) {
      HostMatrix<scalar_type> cdft_Vc_cpu(my_cdft_vars.regions);
      cdft_Vc.resize(my_cdft_vars.regions);
      for (unsigned int i = 0; i < my_cdft_vars.regions; i++) {
        cdft_Vc_cpu(i) = (scalar_type) my_cdft_vars.Vc(i);
      }
      cdft_Vc = cdft_Vc_cpu;
      
      // Apply charge constraints to both alpha and beta factors
      gpu_cdft_factors_accum<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                                  cdft_factors_gpu.data, this->number_of_points,
                                                  my_cdft_vars.regions, cdft_Vc.data, factors_a_gpu.data);
      gpu_cdft_factors_accum<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                                  cdft_factors_gpu.data, this->number_of_points,
                                                  my_cdft_vars.regions, cdft_Vc.data, factors_b_gpu.data);
    }

    // Handle spin constraints if enabled
    if (my_cdft_vars.do_spin) {
      HostMatrix<scalar_type> cdft_Vs_cpu(my_cdft_vars.regions);
      cdft_Vs_a.resize(my_cdft_vars.regions);
      cdft_Vs_b.resize(my_cdft_vars.regions);
      
      // Beta electrons get positive potential
      for (unsigned int i = 0; i < my_cdft_vars.regions; i++) {
        cdft_Vs_cpu(i) = (scalar_type) my_cdft_vars.Vs(i);
      }
      cdft_Vs_b = cdft_Vs_cpu;

      // Alpha electrons get negative potential (opposite sign)
      for (unsigned int i = 0; i < my_cdft_vars.regions; i++) {
        cdft_Vs_cpu(i) = - (scalar_type) my_cdft_vars.Vs(i);
      }
      cdft_Vs_a = cdft_Vs_cpu;

      // Apply spin constraints to factors
      gpu_cdft_factors_accum<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                                  cdft_factors_gpu.data, this->number_of_points,
                                                  my_cdft_vars.regions, cdft_Vs_a.data, factors_a_gpu.data);
      gpu_cdft_factors_accum<scalar_type><<<threadGrid_accumulate, threadBlock_accumulate>>>(
                                                  cdft_factors_gpu.data, this->number_of_points,
                                                  my_cdft_vars.regions, cdft_Vs_b.data, factors_b_gpu.data);
    }

    // Update RMM with optimized kernel selection based on problem size
    // For calls with a single block (pretty common with cubes) don't bother doing the arithmetic to get block position in the matrix
    if (blocksPerRow > 1) {
      gpu_update_rmm<scalar_type,true><<<threadGrid, threadBlock>>>(factors_a_gpu.data, this->number_of_points,
                                                                    rmm_output_a_gpu.data, function_values.data,
                                                                    group_m);
      gpu_update_rmm<scalar_type,true><<<threadGrid, threadBlock>>>(factors_b_gpu.data, this->number_of_points,
                                                                    rmm_output_b_gpu.data, function_values.data,
                                                                    group_m);
    } else {
      gpu_update_rmm<scalar_type,false><<<threadGrid, threadBlock>>>(factors_a_gpu.data, this->number_of_points,
                                                                     rmm_output_a_gpu.data, function_values.data,
                                                                     group_m);
      gpu_update_rmm<scalar_type,false><<<threadGrid, threadBlock>>>(factors_b_gpu.data, this->number_of_points,
                                                                     rmm_output_b_gpu.data, function_values.data,
                                                                     group_m);
    }

    cudaAssertNoError("update_rmm");
    
    /*** Contribute this RMM to the total RMM ***/
    // Transfer results to host and add to global RMMs
    HostMatrix<scalar_type> rmm_output_a_cpu(rmm_output_a_gpu);
    HostMatrix<scalar_type> rmm_output_b_cpu(rmm_output_b_gpu);
    this->add_rmm_output(rmm_output_a_cpu, rmm_output_local_a);
    this->add_rmm_output(rmm_output_b_cpu, rmm_output_local_b);
  }
  timers.rmm.pause_and_sync();

  /* clear functions */
  // Free GPU memory if not needed for later use
  if(!(this->inGlobal)) {
    function_values.deallocate();
    gradient_values.deallocate();
    hessian_values_transposed.deallocate();
  }

  //Deshago el bind de textura de rmm
//  cudaUnbindTexture(rmm_input_gpu_tex); //Enroque el Unbind con el Free, asi parece mas logico. Nano
//  cudaUnbindTexture(rmm_input_gpu_tex2); //Enroque el Unbind con el Free, asi parece mas logico. Nano
//  cudaFreeArray(cuArray1);
//  cudaFreeArray(cuArray2);

  //uint free_memory, total_memory;
  //cudaGetMemoryInfo(free_memory, total_memory);
  //cout << "Maximum used memory: " << (double)max_used_memory / (1024 * 1024) << "MB (" << ((double)max_used_memory / total_memory) * 100.0 << "%)" << endl;
  //cudaPrintMemoryInfo();
}



/*******************************
 * Cube Functions
 *******************************/

/**
 * Computes basis functions, gradients, and Hessians at all points in the group.
 * 
 * This function evaluates all the quantum chemistry basis functions at each
 * spatial point, including their gradients if forces are needed and Hessians
 * if GGA functionals are used.
 * 
 * @param forces Whether to compute gradients for force calculation
 * @param gga Whether to compute Hessians for GGA functional evaluation
 */
template<class scalar_type>
void PointGroupGPU<scalar_type>::compute_functions(bool forces, bool gga)
{
  // Skip if functions are already in global memory
  if(this->inGlobal) //Ya las tengo en memoria? entonces salgo porque ya estan las 3 calculadas
    return;

  // Try to allocate space in global memory pool
  if(0 == GlobalMemoryPool::tryAlloc(this->size_in_gpu())) //1 si hubo error, 0 si pude reservar la memoria
    this->inGlobal=true;
    
  // Matrices for input data
  CudaMatrix<vec_type4> points_position_gpu;
  CudaMatrix<vec_type2> factor_ac_gpu;
  CudaMatrixUInt nuc_gpu;
  CudaMatrixUInt contractions_gpu;

  /** Load points from group **/
  {
    // Transfer point positions to GPU
    HostMatrix<vec_type4> points_position_cpu(this->number_of_points, 1);
    uint i = 0;
    for (vector<Point>::const_iterator p = this->points.begin(); p != this->points.end(); ++p, ++i) {
      points_position_cpu(i) = vec_type4(p->position.x, p->position.y, p->position.z, 0);
    }
    points_position_gpu = points_position_cpu;
  }
  
  /* Load group functions */
  // Calculate total number of functions with all components
  uint group_m = this->s_functions + this->p_functions * 3 + this->d_functions * 6;
  uint4 group_functions = make_uint4(this->s_functions, this->p_functions, this->d_functions, group_m);
  
  // Allocate and prepare contraction coefficient matrices
  HostMatrix<vec_type2> factor_ac_cpu(COALESCED_DIMENSION(group_m), MAX_CONTRACTIONS);
  HostMatrixUInt nuc_cpu(group_m, 1), contractions_cpu(group_m, 1);

  // TODO: hacer que functions.h itere por total_small_functions()... asi puedo hacer que
  // func2global_nuc sea de tamaño total_functions() y directamente copio esa matriz aca y en otros lados

  // Fill matrices with contraction coefficients and nuclear indices
  uint ii = 0;
  for (uint i = 0; i < this->total_functions_simple(); ++i) {
    uint inc = this->small_function_type(i);

    uint func = this->local2global_func[i];
    uint this_nuc = this->func2global_nuc(i);
    uint this_cont = fortran_vars.contractions(func);

    for (uint j = 0; j < inc; j++) {
      nuc_cpu(ii) = this_nuc;
      contractions_cpu(ii) = this_cont;
      for (unsigned int k = 0; k < this_cont; k++)
        factor_ac_cpu(ii, k) = vec_type2(fortran_vars.a_values(func, k), fortran_vars.c_values(func, k));
      ii++;
    }
  }
  
  // Transfer matrices to GPU
  factor_ac_gpu = factor_ac_cpu;
  nuc_gpu = nuc_cpu;
  contractions_gpu = contractions_cpu;

  // Matrix for Hessian values if needed
  CudaMatrix<vec_type<scalar_type,4> > hessian_values;
  
  /** Compute Functions **/
  // Allocate memory for function values and zero
  function_values.resize(COALESCED_DIMENSION(this->number_of_points), group_functions.w);
  function_values.zero();
  
  // Allocate memory for gradients if needed
  if (fortran_vars.do_forces || fortran_vars.gga) {
      gradient_values.resize(COALESCED_DIMENSION(this->number_of_points), group_functions.w);
      gradient_values.zero();
  }
  
  // Allocate memory for Hessians if needed
  if (fortran_vars.gga) {
      hessian_values.resize(COALESCED_DIMENSION(this->number_of_points), (group_functions.w) * 2);
      hessian_values.zero();
  }
  
  // Configure thread organization
  dim3 threads(this->number_of_points);
  dim3 threadBlock(FUNCTIONS_BLOCK_SIZE);
  dim3 threadGrid = divUp(threads, threadBlock);

  // Define parameters for compute_functions kernel
#define compute_functions_parameters \
  points_position_gpu.data,this->number_of_points,contractions_gpu.data,factor_ac_gpu.data,nuc_gpu.data,function_values.data,gradient_values.data,hessian_values.data,group_functions
  
  // Launch appropriate kernel variant based on parameters
  if (forces) {
    if (gga)
      gpu_compute_functions<scalar_type, true, true><<<threadGrid, threadBlock>>>(compute_functions_parameters);
    else
      gpu_compute_functions<scalar_type, true, false><<<threadGrid, threadBlock>>>(compute_functions_parameters);
  }
  else {
    if (gga)
      gpu_compute_functions<scalar_type, false, true><<<threadGrid, threadBlock>>>(compute_functions_parameters);
    else
      gpu_compute_functions<scalar_type, false, false><<<threadGrid, threadBlock>>>(compute_functions_parameters);
  }

  // Wait for kernel to complete
  cudaDeviceSynchronize();

  // Transpose Hessian values for better memory coalescing if needed
  if (fortran_vars.gga) {
    int transposed_width = COALESCED_DIMENSION(this->number_of_points);
    #define BLOCK_DIM 16
    dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);
    dim3 transpose_grid=dim3(transposed_width / BLOCK_DIM, divUp((group_m)*2, BLOCK_DIM), 1);
    hessian_values_transposed.resize((group_m) * 2, COALESCED_DIMENSION(this->number_of_points));
    transpose<<<transpose_grid, transpose_threads>>> (hessian_values_transposed.data,
        hessian_values.data, COALESCED_DIMENSION(this->number_of_points), (group_m)*2);
  }
  cudaAssertNoError("compute_functions");
}

/*******************************
 * Cube Weights
 *******************************/
 
/**
 * Computes integration weights for all points in the group.
 * 
 * This function calculates the appropriate integration weights for numerical
 * quadrature, based on the Becke partitioning scheme.
 */
template<class scalar_type>
void PointGroupGPU<scalar_type>::compute_weights(void)
{
  // Prepare matrices for point and atom data
  CudaMatrix<vec_type4> point_positions_gpu;
  CudaMatrix<vec_type4> atom_position_rm_gpu;
  {
    // Transfer point positions and atom to GPU
    HostMatrix<vec_type4> points_positions_cpu(this->number_of_points, 1);
		uint i = 0;
		for (vector<Point>::const_iterator p = this->points.begin(); p != this->points.end(); ++p, ++i) {
			points_positions_cpu(i) = vec_type4(p->position.x, p->position.y, p->position.z, p->atom);
		}
    point_positions_gpu = points_positions_cpu;

    // Transfer atom positions and radii to GPU
    HostMatrix<vec_type4> atom_position_rm_cpu(fortran_vars.atoms, 1);
    for (uint i = 0; i < fortran_vars.atoms; i++) {
      double3 atom_pos = fortran_vars.atom_positions(i);
      atom_position_rm_cpu(i) = vec_type4(atom_pos.x, atom_pos.y, atom_pos.z, fortran_vars.rm(i));
    }
    atom_position_rm_gpu = atom_position_rm_cpu;
  }

  // Transfer nuclear indices to GPU
  CudaMatrixUInt nucleii_gpu(this->local2global_nuc);

  // Allocate matrix for weights
  CudaMatrix<scalar_type> weights_gpu(this->number_of_points);
  
  // Configure thread organization
  dim3 threads(this->number_of_points);
  dim3 blockSize(WEIGHT_BLOCK_SIZE);
  dim3 gridSize = divUp(threads, blockSize);
  
  // Compute weights on GPU
  gpu_compute_weights<scalar_type><<<gridSize,blockSize>>>(
      this->number_of_points, point_positions_gpu.data, atom_position_rm_gpu.data, weights_gpu.data, nucleii_gpu.data, this->total_nucleii());
  cudaAssertNoError("compute_weights");

  // Transfer weights back to host and apply to points
  HostMatrix<scalar_type> weights_cpu(weights_gpu);
  uint i = 0;
  for (vector<Point>::iterator p =this->points.begin(); p != this->points.end(); ++p, ++i) {
    p->weight *= weights_cpu(i);
    }
}

// Instantiate template classes for single or double precision
#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupGPU<double>;
#else
template class PointGroup<float>;
template class PointGroupGPU<float>;
#endif

}
