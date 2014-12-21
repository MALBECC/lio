/* -*- mode: c -*- */
#include <fstream>
#include <iostream>
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

namespace G2G {
#if FULL_DOUBLE
texture<int2, 2, cudaReadModeElementType> rmm_input_gpu_tex;
texture<int2, 2, cudaReadModeElementType> rmm_input_gpu_tex2;
texture<int2, cudaTextureType2D, cudaReadModeElementType> qmmm_str_tex; // Texture for STR array (used in F(m,U))
#else
texture<float, 2, cudaReadModeElementType> rmm_input_gpu_tex;
texture<float, 2, cudaReadModeElementType> rmm_input_gpu_tex2;
texture<float, cudaTextureType2D, cudaReadModeElementType> qmmm_str_tex;
#endif
/** KERNELS **/
#include "gpu_variables.h"
#include "kernels/pot.h"
#include "kernels/pot_open.h"
#include "kernels/accumulate_point.h"
#include "kernels/energy.h"
#include "kernels/energy_open.h"
#include "kernels/energy_derivs.h"
#include "kernels/rmm.h"
#include "kernels/weight.h"
#include "kernels/functions.h"
#include "kernels/force.h"
#include "kernels/transpose.h"
#include "kernels/qmmm_forces.h"
#include "kernels/qmmm_energy.h"
#include "kernels/coulomb_forces.h"

using std::cout;
using std::vector;
using std::endl;

// Host function to set the constant
void gpu_set_variables(void) {
  cudaMemcpyToSymbol(gpu_normalization_factor, &fortran_vars.normalization_factor, sizeof(fortran_vars.normalization_factor), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(gpu_atoms, &fortran_vars.atoms, sizeof(fortran_vars.atoms), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(gpu_Iexch, &fortran_vars.iexch, sizeof(fortran_vars.iexch), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(gpu_m, &fortran_vars.m, sizeof(fortran_vars.m), 0, cudaMemcpyHostToDevice);

  // This is needed by d-d QM/MM forces calculations to know which orbital a thread maps to
  uint d_offset = fortran_vars.s_funcs + fortran_vars.p_funcs*3;
  cudaMemcpyToSymbol(gpu_d_offset, &d_offset, sizeof(d_offset), 0, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(gpu_dens_gauss, &fortran_vars.gaussians_dens, sizeof(fortran_vars.gaussians_dens), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(gpu_dens_s_gauss, &fortran_vars.s_gaussians_dens, sizeof(fortran_vars.s_gaussians_dens), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(gpu_dens_p_gauss, &fortran_vars.p_gaussians_dens, sizeof(fortran_vars.p_gaussians_dens), 0, cudaMemcpyHostToDevice);

  cudaAssertNoError("set_gpu_variables");
}

//
// Set up arrays needed for F(m,U) calculation in QM/MM kernels (STR and FAC) and send them to device
// FAC is small so it's put into constant memory
// STR is large and accessed (potentially) with a random access pattern in the first index
// TODO: putting STR into texture for now; need to see if there's a better way to access it in the kernel
//
template<class scalar_type>
void gpu_set_gamma_arrays()
{
  // Cast STR/FAC to appropriate type (float/double)
  HostMatrix<scalar_type> h_str(880,22), h_fac(17);
  for (uint i = 0; i < 880; i++) {
    for (uint j = 0; j < 22; j++) {
      h_str(i,j) = fortran_vars.str(i,j);
    }
  }
  for (uint i = 0; i < 17; i++) {
    h_fac(i) = fortran_vars.fac(i);
  }

  qmmm_str_tex.normalized = false;
  qmmm_str_tex.filterMode = cudaFilterModePoint;
  cudaMallocArray(&gammaArray,&qmmm_str_tex.channelDesc,880,22);//GAMMA_LENGTH,6);
  cudaMemcpyToArray(gammaArray,0,0,h_str.data,sizeof(scalar_type)*880*22,cudaMemcpyHostToDevice);
  //scalar_type* d_str_ptr;
  //cudaMalloc((void**)&d_str_ptr,880*22*sizeof(scalar_type));
  // STR data host->device
  //cudaMemcpy(d_str_ptr,h_str.data,h_str.bytes(),cudaMemcpyHostToDevice);
  // STR device pointer h->d
  //cudaMemcpyToSymbol(gpu_str,&d_str_ptr,sizeof(gpu_str),0,cudaMemcpyHostToDevice);

  // FAC data h->d
  cudaMemcpyToSymbol(gpu_fac,h_fac.data,h_fac.bytes(),0,cudaMemcpyHostToDevice);

  cudaAssertNoError("gpu_set_gamma_arrays");
}

template<class T> void gpu_set_atom_positions(const HostMatrix<T>& m) {
  cudaMemcpyToSymbol(gpu_atom_positions, m.data, m.bytes(), 0, cudaMemcpyHostToDevice);
}

//
// Tell the GPU how many MM atoms are being used this step
//
void gpu_set_clatoms(void)
{
  cudaMemcpyToSymbol(gpu_clatoms, &fortran_vars.clatoms, sizeof(fortran_vars.clatoms), 0, cudaMemcpyHostToDevice);

  cudaAssertNoError("gpu_set_clatoms");
}

#if FULL_DOUBLE
template void gpu_set_gamma_arrays<double>( void );
#else
template void gpu_set_gamma_arrays<float>( void );
#endif
template void gpu_set_atom_positions<double3>(const HostMatrix<double3>& m);
template void gpu_set_atom_positions<float3>(const HostMatrix<float3>& m);
//template<class scalar_type,true> __global__ void gpu_update_rmm(scalar_type* factors, uint points, scalar_type* rmm, scalar_type* function_values, uint m);
//template<class scalar_type,false> __global__ void gpu_update_rmm(scalar_type* factors, uint points, scalar_type* rmm, scalar_type* function_values, uint m);

template<class scalar_type>
void PointGroup<scalar_type>::solve(Timers& timers, bool compute_rmm, bool lda, bool compute_forces,
    bool compute_energy, double& energy,double& energy_i, double& energy_c, double& energy_c1,
    double& energy_c2, double* fort_forces_ptr, bool open){
  if(open) {
    solve_opened(timers, compute_rmm, lda, compute_forces, compute_energy, energy, energy_i, energy_c, energy_c1,
        energy_c2, fort_forces_ptr);
  }
  else {
    solve_closed(timers, compute_rmm, lda, compute_forces, compute_energy, energy, fort_forces_ptr);
  }

}

template<class scalar_type>
void PointGroup<scalar_type>::solve_closed(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, bool compute_energy, double& energy, double* fort_forces_ptr){
  //uint max_used_memory = 0;

  /*** Computo sobre cada cubo ****/
  CudaMatrix<scalar_type> point_weights_gpu;
  FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, fortran_vars.max_atoms);

  /** Compute this group's functions **/
  timers.functions.start_and_sync();
  compute_functions(compute_forces, !lda);
  timers.functions.pause_and_sync();

  uint group_m = total_functions();

  timers.density.start_and_sync();
  /** Load points from group **/
  HostMatrix<scalar_type> point_weights_cpu(number_of_points, 1);

  uint i = 0;
  for (vector<Point>::const_iterator p = points.begin(); p != points.end(); ++p, ++i) {
    point_weights_cpu(i) = p->weight;
  }
  point_weights_gpu = point_weights_cpu;

  dim3 threadBlock, threadGrid;
  /* compute density/factors */

  const int block_height= divUp(group_m, 2*DENSITY_BLOCK_SIZE);

  threadBlock = dim3(DENSITY_BLOCK_SIZE,1,1); // Hay que asegurarse que la cantidad de funciones este en rango
  threadGrid = dim3(number_of_points,block_height,1);

  CudaMatrix<scalar_type> factors_gpu;

  CudaMatrix<scalar_type> partial_densities_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dxyz_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd1_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd2_gpu;

  partial_densities_gpu.resize(COALESCED_DIMENSION(number_of_points), block_height);
  dxyz_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height);
  dd1_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );
  dd2_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );

  const dim3 threadGrid_accumulate(divUp(number_of_points,DENSITY_ACCUM_BLOCK_SIZE),1,1);
  const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE,1,1);

  if (compute_rmm || compute_forces)
    factors_gpu.resize(number_of_points);

  int transposed_width = COALESCED_DIMENSION(number_of_points);
  #define BLOCK_DIM 16
  dim3 transpose_grid(transposed_width / BLOCK_DIM, divUp((group_m),BLOCK_DIM));
  dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);

  CudaMatrix<scalar_type> function_values_transposed;
  CudaMatrix<vec_type<scalar_type,4> > gradient_values_transposed;

  function_values_transposed.resize(group_m, COALESCED_DIMENSION(number_of_points));

  if (fortran_vars.do_forces || fortran_vars.gga)
      gradient_values_transposed.resize( group_m,COALESCED_DIMENSION(number_of_points));

  transpose<<<transpose_grid, transpose_threads>>> (function_values_transposed.data,
      function_values.data, COALESCED_DIMENSION(number_of_points), group_m);

  if (fortran_vars.do_forces || fortran_vars.gga)
    transpose<<<transpose_grid, transpose_threads>>> (gradient_values_transposed.data,
        gradient_values.data, COALESCED_DIMENSION(number_of_points), group_m );


  HostMatrix<scalar_type> rmm_input_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
  get_rmm_input(rmm_input_cpu); //Achica la matriz densidad a la version reducida del grupo

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


  /*
   **********************************************************************
   * Pasando RDM (rmm) a texturas
   **********************************************************************
   */

  //Comentado porque ahora vamos a hacer esto a mano por la textura
  // TODO: pasarlo a un metodo dentro de matrix.cpp
  //rmm_input_gpu = rmm_input_cpu; //Aca copia de CPU a GPU

  cudaArray* cuArray;
  cudaMallocArray(&cuArray, &rmm_input_gpu_tex.channelDesc, rmm_input_cpu.width, rmm_input_cpu.height);
  cudaMemcpyToArray(cuArray, 0, 0, rmm_input_cpu.data, sizeof(scalar_type)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);
  cudaBindTextureToArray(rmm_input_gpu_tex, cuArray);

  rmm_input_gpu_tex.normalized = false;

  if (compute_energy) {
    CudaMatrix<scalar_type> energy_gpu(number_of_points);
#define compute_parameters \
    energy_gpu.data,factors_gpu.data,point_weights_gpu.data,number_of_points,function_values_transposed.data,gradient_values_transposed.data,hessian_values_transposed.data,group_m,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data
#define accumulate_parameters \
    energy_gpu.data,factors_gpu.data,point_weights_gpu.data,number_of_points,block_height,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data
    if (compute_forces || compute_rmm) {
      if (lda)
      {
          gpu_compute_density<scalar_type, true, true, true><<<threadGrid, threadBlock>>>(compute_parameters);
          gpu_accumulate_point<scalar_type, true, true, true><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
      }
      else
      {
          gpu_compute_density<scalar_type, true, true, false><<<threadGrid, threadBlock>>>(compute_parameters);
          gpu_accumulate_point<scalar_type, true, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
      }
    }
    else {
      if (lda)
      {
          gpu_compute_density<scalar_type, true, false, true><<<threadGrid, threadBlock>>>(compute_parameters);
          gpu_accumulate_point<scalar_type, true, false, true><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
      }
      else
      {
          gpu_compute_density<scalar_type, true, false, false><<<threadGrid, threadBlock>>>(compute_parameters);
          gpu_accumulate_point<scalar_type, true, false, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
      }
    }
    cudaAssertNoError("compute_density");

    HostMatrix<scalar_type> energy_cpu(energy_gpu);
    for (uint i = 0; i < number_of_points; i++) {
      energy += energy_cpu(i);
    } // TODO: hacer con un kernel?
  }
  else {
#undef compute_parameters
#undef accumulate_parameters

#define compute_parameters \
    NULL,factors_gpu.data,point_weights_gpu.data,number_of_points,function_values_transposed.data,gradient_values_transposed.data,hessian_values_transposed.data,group_m,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data
#define accumulate_parameters \
    NULL,factors_gpu.data,point_weights_gpu.data,number_of_points,block_height,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data
    if (lda)
    {
        gpu_compute_density<scalar_type, false, true, true><<<threadGrid, threadBlock>>>(compute_parameters);
        gpu_accumulate_point<scalar_type, false, true, true><<<threadGrid_accumulate, threadBlock_accumulate>>>(accumulate_parameters);
    }
    else
    {
        gpu_compute_density<scalar_type, false, true, false><<<threadGrid, threadBlock>>>(compute_parameters);
        gpu_accumulate_point<scalar_type, false, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>>(accumulate_parameters);
    }
    cudaAssertNoError("compute_density");
  }
#undef compute_parameters
#undef accumulate_parameters

  timers.density.pause_and_sync();

//************ Repongo los valores que puse a cero antes, para las fuerzas son necesarios (o por lo mens utiles)
  for (uint i=0; i<(group_m); i++) {
    for(uint j=0; j<(group_m); j++) {
      if((i>=group_m) || (j>=group_m) || (j > i))
      {
        rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*j+i] ;
      }
    }
  }

  cudaMemcpyToArray(cuArray, 0, 0,rmm_input_cpu.data,
      sizeof(scalar_type)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);

   dim3 threads;
  /* compute forces */
  if (compute_forces) {
    timers.density_derivs.start_and_sync();
    threads = dim3(number_of_points);
    threadBlock = dim3(DENSITY_DERIV_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);

    CudaMatrix<vec_type4> dd_gpu(COALESCED_DIMENSION(number_of_points), total_nucleii()); dd_gpu.zero();
    CudaMatrixUInt nuc_gpu(func2local_nuc);  // TODO: esto en realidad se podria guardar una sola vez durante su construccion

    gpu_compute_density_derivs<<<threadGrid, threadBlock>>>(
        function_values.data, gradient_values.data, nuc_gpu.data, dd_gpu.data, number_of_points, group_m, total_nucleii());
    cudaAssertNoError("density_derivs");
    timers.density_derivs.pause_and_sync();

    timers.forces.start_and_sync();
    CudaMatrix<vec_type4> forces_gpu(total_nucleii());

    threads = dim3(total_nucleii());
    threadBlock = dim3(FORCE_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);
    gpu_compute_forces<<<threadGrid, threadBlock>>>(
        number_of_points, factors_gpu.data, dd_gpu.data, forces_gpu.data, total_nucleii());
    cudaAssertNoError("forces");

    HostMatrix<vec_type4> forces_cpu(forces_gpu);

    for (uint i = 0; i < total_nucleii(); ++i) {
      vec_type4 atom_force = forces_cpu(i);
      uint global_nuc = local2global_nuc[i];
      fort_forces(global_nuc, 0) += atom_force.x;
      fort_forces(global_nuc, 1) += atom_force.y;
      fort_forces(global_nuc, 2) += atom_force.z;
    }
    timers.forces.pause_and_sync();
  }

  timers.rmm.start_and_sync();
  /* compute RMM */
  if (compute_rmm) {
    threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);
    uint blocksPerRow = divUp(group_m, RMM_BLOCK_SIZE_XY);
    // Only use enough blocks for lower triangle
    threadGrid = dim3(blocksPerRow*(blocksPerRow+1)/2);

    CudaMatrix<scalar_type> rmm_output_gpu(COALESCED_DIMENSION(group_m), group_m);
    // For calls with a single block (pretty common with cubes) don't bother doing the arithmetic to get block position in the matrix
    if (blocksPerRow > 1) {
        gpu_update_rmm<scalar_type,true><<<threadGrid, threadBlock>>>(factors_gpu.data, number_of_points, rmm_output_gpu.data, function_values.data, group_m);
    } else {
        gpu_update_rmm<scalar_type,false><<<threadGrid, threadBlock>>>(factors_gpu.data, number_of_points, rmm_output_gpu.data, function_values.data, group_m);
    }
    cudaAssertNoError("update_rmm");

    /*** Contribute this RMM to the total RMM ***/
    HostMatrix<scalar_type> rmm_output_cpu(rmm_output_gpu);
    add_rmm_output(rmm_output_cpu);
  }
  timers.rmm.pause_and_sync();

  /* clear functions */
  if(!(this->inGlobal)) {
    function_values.deallocate();
    gradient_values.deallocate();
    hessian_values_transposed.deallocate();
  }
  //Deshago el bind de textura de rmm
  cudaUnbindTexture(rmm_input_gpu_tex); //Enroque el Unbind con el Free, asi parece mas logico. Nano
  cudaFreeArray(cuArray);
}

//======================
// OPENSHELL
//======================
template<class scalar_type>
void PointGroup<scalar_type>::solve_opened(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, bool compute_energy, double& energy, double& energy_i, double& energy_c, double& energy_c1, double& energy_c2, double* fort_forces_ptr){

  /*** Computo sobre cada cubo ****/
  CudaMatrix<scalar_type> point_weights_gpu;
  FortranMatrix<double> fort_forces(fort_forces_ptr, fortran_vars.atoms, 3, fortran_vars.max_atoms);

  /** Compute this group's functions **/
  timers.functions.start_and_sync();
  compute_functions(compute_forces, !lda); //<<<<==========
  timers.functions.pause_and_sync();

  uint group_m = total_functions();

  timers.density.start_and_sync();
  /** Load points from group **/
  HostMatrix<scalar_type> point_weights_cpu(number_of_points, 1);

  uint i = 0;
  for (vector<Point>::const_iterator p = points.begin(); p != points.end(); ++p, ++i) {
    point_weights_cpu(i) = p->weight;
  }
  point_weights_gpu = point_weights_cpu;

//<<===========================>>//
  dim3 threadBlock, threadGrid;
  /* compute density/factors */

  const int block_height= divUp(group_m,2*DENSITY_BLOCK_SIZE);

  threadBlock = dim3(DENSITY_BLOCK_SIZE,1,1); // Hay que asegurarse que la cantidad de funciones este en rango
  threadGrid = dim3(number_of_points,block_height,1);

  CudaMatrix<scalar_type> factors_a_gpu;
  CudaMatrix<scalar_type> factors_b_gpu;

  /*
  dxyz_gpu; gradiente
  dd1_gpu;  hessiano ii
  dd2_gpu;  hessiano ij
  */

  CudaMatrix<scalar_type> partial_densities_a_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dxyz_a_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd1_a_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd2_a_gpu;

  CudaMatrix<scalar_type> partial_densities_b_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dxyz_b_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd1_b_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd2_b_gpu;

  /*
   **********************************************************************
   * Transposiciones de matrices para la coalescencia mejorada en density
   **********************************************************************
   */

  CudaMatrix<scalar_type> function_values_transposed;
  CudaMatrix<vec_type<scalar_type,4> > gradient_values_transposed;

  int transposed_width = COALESCED_DIMENSION(number_of_points);

  function_values_transposed.resize(group_m, COALESCED_DIMENSION(number_of_points));
  if (fortran_vars.do_forces || fortran_vars.gga)
      gradient_values_transposed.resize( group_m,COALESCED_DIMENSION(number_of_points));

  #define BLOCK_DIM 16
  dim3 transpose_grid(transposed_width / BLOCK_DIM, divUp((group_m),BLOCK_DIM));
  dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);

  transpose<<<transpose_grid, transpose_threads>>> (function_values_transposed.data, function_values.data,  COALESCED_DIMENSION(number_of_points),group_m   );
  if (fortran_vars.do_forces || fortran_vars.gga)
      transpose<<<transpose_grid, transpose_threads>>> (gradient_values_transposed.data, gradient_values.data, COALESCED_DIMENSION(number_of_points), group_m );

  partial_densities_a_gpu.resize(COALESCED_DIMENSION(number_of_points), block_height);
  dxyz_a_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height);
  dd1_a_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );
  dd2_a_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );

  partial_densities_b_gpu.resize(COALESCED_DIMENSION(number_of_points), block_height);
  dxyz_b_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height);
  dd1_b_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );
  dd2_b_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );

  const dim3 threadGrid_accumulate(divUp(number_of_points,DENSITY_ACCUM_BLOCK_SIZE),1,1);
  const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE,1,1);

  if (compute_rmm || compute_forces) {
  	factors_a_gpu.resize(number_of_points);
  	factors_b_gpu.resize(number_of_points);
  }

  HostMatrix<scalar_type> rmm_input_a_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
  HostMatrix<scalar_type> rmm_input_b_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
  get_rmm_input(rmm_input_a_cpu,rmm_input_b_cpu); //Achica las matrices densidad (Up,Down) a la version reducida del grupo

  for (uint i=0; i<(group_m+DENSITY_BLOCK_SIZE); i++) {
    for(uint j=0; j<COALESCED_DIMENSION(group_m); j++) {
      if((i>=group_m) || (j>=group_m) || (j > i)) {
        rmm_input_a_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
        rmm_input_b_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
      }
    }
  }

  /*
  **********************************************************************
  * Pasando RDM (rmm) a texturas/
  **********************************************************************
  */


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

  if (compute_energy) {
    CudaMatrix<scalar_type> energy_gpu(number_of_points);
    CudaMatrix<scalar_type> energy_i_gpu(number_of_points);
    CudaMatrix<scalar_type> energy_c_gpu(number_of_points);
    CudaMatrix<scalar_type> energy_c1_gpu(number_of_points);
    CudaMatrix<scalar_type> energy_c2_gpu(number_of_points);

    if (compute_forces || compute_rmm) {
      gpu_compute_density_opened<scalar_type, true, true, false><<<threadGrid, threadBlock>>>(
             point_weights_gpu.data,number_of_points, function_values_transposed.data,
             gradient_values_transposed.data,hessian_values_transposed.data, group_m,
             partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
             partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data);

      gpu_accumulate_point_open<scalar_type, true, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
             energy_gpu.data,energy_i_gpu.data,energy_c_gpu.data,energy_c1_gpu.data,energy_c2_gpu.data,
             factors_a_gpu.data, factors_b_gpu.data, point_weights_gpu.data,number_of_points,block_height,
             partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
             partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data);
    }
    else {
      gpu_compute_density_opened<scalar_type, true, false, false><<<threadGrid, threadBlock>>>(
             point_weights_gpu.data,number_of_points, function_values_transposed.data,
             gradient_values_transposed.data,hessian_values_transposed.data, group_m,
             partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
             partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data);
      gpu_accumulate_point_open<scalar_type, true, false, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
             energy_gpu.data, energy_i_gpu.data,energy_c_gpu.data,energy_c1_gpu.data,energy_c2_gpu.data,
             factors_a_gpu.data, factors_b_gpu.data, point_weights_gpu.data,number_of_points,block_height,
             partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
             partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data);
    }
    cudaAssertNoError("compute_density");

    HostMatrix<scalar_type> energy_cpu(energy_gpu);
    HostMatrix<scalar_type> energy_i_cpu(energy_i_gpu);
    HostMatrix<scalar_type> energy_c_cpu(energy_c_gpu);
    HostMatrix<scalar_type> energy_c1_cpu(energy_c1_gpu);
    HostMatrix<scalar_type> energy_c2_cpu(energy_c2_gpu);
    for (uint i = 0; i < number_of_points; i++) {
      energy    += energy_cpu(i);
      energy_i  += energy_i_cpu(i);
      energy_c  += energy_c_cpu(i);
      energy_c1 += energy_c1_cpu(i);
      energy_c2 += energy_c2_cpu(i);
    } // TODO: hacer con un kernel?
  }
  else {
    gpu_compute_density_opened<scalar_type, false, true, false><<<threadGrid, threadBlock>>>(
           point_weights_gpu.data, number_of_points, function_values_transposed.data,
           gradient_values_transposed.data,hessian_values_transposed.data, group_m,
           partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
           partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data);
    gpu_accumulate_point_open<scalar_type, false, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (
           NULL,NULL,NULL,NULL,NULL,
           factors_a_gpu.data, factors_b_gpu.data, point_weights_gpu.data,number_of_points,block_height,
           partial_densities_a_gpu.data, dxyz_a_gpu.data, dd1_a_gpu.data, dd2_a_gpu.data,
           partial_densities_b_gpu.data, dxyz_b_gpu.data, dd1_b_gpu.data, dd2_b_gpu.data);
    cudaAssertNoError("compute_density");
  }

  timers.density.pause_and_sync();

//************ Repongo los valores que puse a cero antes, para las fuerzas son necesarios (o por lo menos utiles)
  for (uint i=0; i<(group_m); i++) {
    for(uint j=0; j<(group_m); j++) {
      if((i>=group_m) || (j>=group_m) || (j > i)){
        rmm_input_a_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=rmm_input_a_cpu.data[COALESCED_DIMENSION(group_m)*j+i] ;
        rmm_input_b_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=rmm_input_b_cpu.data[COALESCED_DIMENSION(group_m)*j+i] ;
      }
    }
  }

  cudaMemcpyToArray(cuArray1, 0, 0,rmm_input_a_cpu.data,sizeof(scalar_type)*rmm_input_a_cpu.width*rmm_input_a_cpu.height, cudaMemcpyHostToDevice);
  cudaMemcpyToArray(cuArray2, 0, 0,rmm_input_b_cpu.data,sizeof(scalar_type)*rmm_input_b_cpu.width*rmm_input_b_cpu.height, cudaMemcpyHostToDevice);

//**********************************************

   dim3 threads;
  /* compute forces */

  if (compute_forces) {
    timers.density_derivs.start_and_sync();
    threads = dim3(number_of_points);
    threadBlock = dim3(DENSITY_DERIV_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);

    CudaMatrix<vec_type4> dd_gpu_a(COALESCED_DIMENSION(number_of_points), total_nucleii());
    CudaMatrix<vec_type4> dd_gpu_b(COALESCED_DIMENSION(number_of_points), total_nucleii());
    dd_gpu_a.zero();
    dd_gpu_b.zero();
    CudaMatrixUInt nuc_gpu(func2local_nuc);  // TODO: esto en realidad se podria guardar una sola vez durante su construccion

    // Kernel
    gpu_compute_density_derivs_open<<<threadGrid, threadBlock>>>(function_values.data, gradient_values.data, nuc_gpu.data, dd_gpu_a.data, dd_gpu_b.data, number_of_points, group_m, total_nucleii());

    cudaAssertNoError("density_derivs");
    timers.density_derivs.pause_and_sync();

    timers.forces.start_and_sync();
    CudaMatrix<vec_type4> forces_gpu_a(total_nucleii());
    CudaMatrix<vec_type4> forces_gpu_b(total_nucleii());

    threads = dim3(total_nucleii());
    threadBlock = dim3(FORCE_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);
    // Kernel
    gpu_compute_forces<<<threadGrid, threadBlock>>>(number_of_points, factors_a_gpu.data, dd_gpu_a.data, forces_gpu_a.data, total_nucleii());
    gpu_compute_forces<<<threadGrid, threadBlock>>>(number_of_points, factors_b_gpu.data, dd_gpu_b.data, forces_gpu_b.data, total_nucleii());

    cudaAssertNoError("forces");

    HostMatrix<vec_type4> forces_cpu_a(forces_gpu_a);
    HostMatrix<vec_type4> forces_cpu_b(forces_gpu_b);

    for (uint i = 0; i < total_nucleii(); ++i) {
      vec_type4 atom_force_a = forces_cpu_a(i);
      vec_type4 atom_force_b = forces_cpu_b(i);
      uint global_nuc = local2global_nuc[i];

      fort_forces(global_nuc, 0)=fort_forces(global_nuc, 0) + atom_force_a.x + atom_force_b.x;
      fort_forces(global_nuc, 1)=fort_forces(global_nuc, 1) + atom_force_a.y + atom_force_b.y;
      fort_forces(global_nuc, 2)=fort_forces(global_nuc, 2) + atom_force_a.z + atom_force_b.z;
    }

    timers.forces.pause_and_sync();
  }

  /* compute RMM */
  timers.rmm.start_and_sync();
  if (compute_rmm) {
	//cout<<"EMPEZANDO GPU_UPDATE_RMM"<<endl;
	//threads = dim3(group_m, group_m);
    	threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);
    	//threadGrid = divUp(threads, threadBlock);
        uint blocksPerRow = divUp(group_m, RMM_BLOCK_SIZE_XY);
        // Only use enough blocks for lower triangle
        threadGrid = dim3(blocksPerRow*(blocksPerRow+1)/2);

    	CudaMatrix<scalar_type> rmm_output_a_gpu(COALESCED_DIMENSION(group_m), group_m);
    	CudaMatrix<scalar_type> rmm_output_b_gpu(COALESCED_DIMENSION(group_m), group_m);
    	// Kernel
//	cout<<"alpha"<<endl;
        // For calls with a single block (pretty common with cubes) don't bother doing the arithmetic to get block position in the matrix
        if (blocksPerRow > 1) {
	    gpu_update_rmm<scalar_type,true><<<threadGrid, threadBlock>>>(factors_a_gpu.data, number_of_points, rmm_output_a_gpu.data, function_values.data, group_m);
//	cout<<endl;
//        cout<<"beta"<<endl;
            gpu_update_rmm<scalar_type,true><<<threadGrid, threadBlock>>>(factors_b_gpu.data, number_of_points, rmm_output_b_gpu.data, function_values.data, group_m);
        } else {
	    gpu_update_rmm<scalar_type,false><<<threadGrid, threadBlock>>>(factors_a_gpu.data, number_of_points, rmm_output_a_gpu.data, function_values.data, group_m);
            gpu_update_rmm<scalar_type,false><<<threadGrid, threadBlock>>>(factors_b_gpu.data, number_of_points, rmm_output_b_gpu.data, function_values.data, group_m);
        }
    	//cout<<endl;
        cudaAssertNoError("update_rmm");

    	/*** Contribute this RMM to the total RMM ***/
    	HostMatrix<scalar_type> rmm_output_a_cpu(rmm_output_a_gpu);
    	HostMatrix<scalar_type> rmm_output_b_cpu(rmm_output_b_gpu);
    	//add_rmm_open_output(rmm_output_a_cpu,rmm_output_b_cpu);
    	add_rmm_output_a(rmm_output_a_cpu);
    	add_rmm_output_b(rmm_output_b_cpu);

    //CudaMatrix<scalar_type> rmm_output_a_gpu(COALESCED_DIMENSION(group_m), group_m);
//    CudaMatrix<scalar_type> rmm_output_b_gpu(COALESCED_DIMENSION(group_m), group_m);
//    gpu_update_rmm<scalar_type,true><<<threadGrid, threadBlock>>>(factors_a_gpu.data, number_of_points, rmm_output_a_gpu.data, function_values.data, group_m);
//    gpu_update_rmm<scalar_type,true><<<threadGrid, threadBlock>>>(factors_b_gpu.data, number_of_points, rmm_output_b_gpu.data, function_values.data, group_m);
//    cudaAssertNoError("update_rmm");

    /*** Contribute this RMM to the total RMM ***/
//    HostMatrix<scalar_type> rmm_output_a_cpu(rmm_output_a_gpu);
//    HostMatrix<scalar_type> rmm_output_b_cpu(rmm_output_b_gpu);
    //add_rmm_open_output(rmm_output_a_cpu,rmm_output_b_cpu);
//    add_rmm_output_a(rmm_output_a_cpu);
//    add_rmm_output_b(rmm_output_b_cpu);
  }
  timers.rmm.pause_and_sync();

  /* clear functions */
  if(!(this->inGlobal)) {
    function_values.deallocate();
    gradient_values.deallocate();
    hessian_values_transposed.deallocate();
  }

  //Deshago el bind de textura de rmm
  cudaUnbindTexture(rmm_input_gpu_tex); //Enroque el Unbind con el Free, asi parece mas logico. Nano
  cudaUnbindTexture(rmm_input_gpu_tex2); //Enroque el Unbind con el Free, asi parece mas logico. Nano
  cudaFreeArray(cuArray1);
  cudaFreeArray(cuArray2);
}

/*******************************
 * Cube Functions
 *******************************/

template<class scalar_type>
void PointGroup<scalar_type>::compute_functions(bool forces, bool gga)
{
  if(this->inGlobal) //Ya las tengo en memoria? entonces salgo porque ya estan las 3 calculadas
    return;

  if(0 == globalMemoryPool::tryAlloc(this->size_in_gpu())) //1 si hubo error, 0 si pude reservar la memoria
    this->inGlobal=true;
  CudaMatrix<vec_type4> points_position_gpu;
  CudaMatrix<vec_type2> factor_ac_gpu;
  CudaMatrixUInt nuc_gpu;
  CudaMatrixUInt contractions_gpu;

  /** Load points from group **/
  {
    HostMatrix<vec_type4> points_position_cpu(number_of_points, 1);
    uint i = 0;
    for (vector<Point>::const_iterator p = points.begin(); p != points.end(); ++p, ++i) {
      points_position_cpu(i) = vec_type4(p->position.x, p->position.y, p->position.z, 0);
    }
    points_position_gpu = points_position_cpu;
  }
  /* Load group functions */
  uint group_m = s_functions + p_functions * 3 + d_functions * 6;
  uint4 group_functions = make_uint4(s_functions, p_functions, d_functions, group_m);
  HostMatrix<vec_type2> factor_ac_cpu(COALESCED_DIMENSION(group_m), MAX_CONTRACTIONS);
  HostMatrixUInt nuc_cpu(group_m, 1), contractions_cpu(group_m, 1);

  // TODO: hacer que functions.h itere por total_small_functions()... asi puedo hacer que
  // func2global_nuc sea de tamaño total_functions() y directamente copio esa matriz aca y en otros lados

  uint ii = 0;
  for (uint i = 0; i < total_functions_simple(); ++i) {
    uint inc = small_function_type(i);

    uint func = local2global_func[i];
    uint this_nuc = func2global_nuc(i);
    uint this_cont = fortran_vars.contractions(func);

    for (uint j = 0; j < inc; j++) {
      nuc_cpu(ii) = this_nuc;
      contractions_cpu(ii) = this_cont;
      for (unsigned int k = 0; k < this_cont; k++)
        factor_ac_cpu(ii, k) = vec_type2(fortran_vars.a_values(func, k), fortran_vars.c_values(func, k));
      ii++;
    }
  }
  factor_ac_gpu = factor_ac_cpu;
  nuc_gpu = nuc_cpu;
  contractions_gpu = contractions_cpu;

  CudaMatrix<vec_type<scalar_type,4> > hessian_values;
  /** Compute Functions **/
  function_values.resize(COALESCED_DIMENSION(number_of_points), group_functions.w);
  if (fortran_vars.do_forces || fortran_vars.gga)
      gradient_values.resize(COALESCED_DIMENSION(number_of_points), group_functions.w);
  if (fortran_vars.gga)
      hessian_values.resize(COALESCED_DIMENSION(number_of_points), (group_functions.w) * 2);

  dim3 threads(number_of_points);
  dim3 threadBlock(FUNCTIONS_BLOCK_SIZE);
  dim3 threadGrid = divUp(threads, threadBlock);

 // cout << "points: " << threads.x << " " << threadGrid.x << " " << threadBlock.x << endl;
#define compute_functions_parameters \
  points_position_gpu.data,number_of_points,contractions_gpu.data,factor_ac_gpu.data,nuc_gpu.data,function_values.data,gradient_values.data,hessian_values.data,group_functions
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

  if (fortran_vars.gga) {
    int transposed_width = COALESCED_DIMENSION(number_of_points);
    #define BLOCK_DIM 16
    dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);
    dim3 transpose_grid=dim3(transposed_width / BLOCK_DIM, divUp((group_m)*2, BLOCK_DIM), 1);

    hessian_values_transposed.resize((group_m) * 2, COALESCED_DIMENSION(number_of_points));
    transpose<<<transpose_grid, transpose_threads>>> (hessian_values_transposed.data,
        hessian_values.data, COALESCED_DIMENSION(number_of_points), (group_m)*2);
  }
  cudaAssertNoError("compute_functions");
}

/*******************************
 * Cube Weights
 *******************************/
template<class scalar_type>
void PointGroup<scalar_type>::compute_weights(void)
{
  CudaMatrix<vec_type4> point_positions_gpu;
  CudaMatrix<vec_type4> atom_position_rm_gpu;
  {
    HostMatrix<vec_type4> points_positions_cpu(number_of_points, 1);
		uint i = 0;
		for (vector<Point>::const_iterator p = points.begin(); p != points.end(); ++p, ++i) {
			points_positions_cpu(i) = vec_type4(p->position.x, p->position.y, p->position.z, p->atom);
		}
    point_positions_gpu = points_positions_cpu;

    HostMatrix<vec_type4> atom_position_rm_cpu(fortran_vars.atoms, 1);
    for (uint i = 0; i < fortran_vars.atoms; i++) {
      double3 atom_pos = fortran_vars.atom_positions(i);
      atom_position_rm_cpu(i) = vec_type4(atom_pos.x, atom_pos.y, atom_pos.z, fortran_vars.rm(i));
    }
    atom_position_rm_gpu = atom_position_rm_cpu;
  }

  CudaMatrixUInt nucleii_gpu(local2global_nuc);

  CudaMatrix<scalar_type> weights_gpu(number_of_points);
  dim3 threads(number_of_points);
  dim3 blockSize(WEIGHT_BLOCK_SIZE);
  dim3 gridSize = divUp(threads, blockSize);
  gpu_compute_weights<scalar_type><<<gridSize,blockSize>>>(
      number_of_points, point_positions_gpu.data, atom_position_rm_gpu.data, weights_gpu.data, nucleii_gpu.data, total_nucleii());
  cudaAssertNoError("compute_weights");

  #if REMOVE_ZEROS
  std::vector<Point> nonzero_points;
  uint nonzero_number_of_points = 0;
  #endif

  uint ceros = 0;

  HostMatrix<scalar_type> weights_cpu(weights_gpu);
  uint i = 0;
  for (vector<Point>::iterator p = points.begin(); p != points.end(); ++p, ++i) {
    p->weight *= weights_cpu(i);

    if (p->weight == 0.0) {
      ceros++;
    }
    #if REMOVE_ZEROS
    else {
      nonzero_points.push_back(*p);
      nonzero_number_of_points++;
    }
    #endif
  }

  //cout << "ceros: " << ceros << "/" << group.number_of_points << " (" << (ceros / (double)group.number_of_points) * 100 << "%)" << endl;

  #if REMOVE_ZEROS
  points = nonzero_points;
  number_of_points = nonzero_number_of_points;
  #endif
}

template class PointGroup<double>;
template class PointGroup<float>;

#define NUM_TERM_TYPES 6 // 6 types when using s,p,and d functions: s-s,p-s,p-p,d-s,d-p,d-d
#define MAX_TERM_TYPE 6

//
// Main QM/MM routine
// If forces = true, calculate gradients of QM/MM operator and return in qm_forces and mm_forces
// If forces = false, calculate Fock matrix elements of QM/MM operator (returned in RMM(M11) back in Fortran)
//                    and QM/MM energy of the current density (nuc-nuc in Ens, e-nuc in Es)
//
template <class scalar_type,bool forces> void g2g_qmmm(double* qm_forces, double* mm_forces, double& Ens, double& Es)
{
  uint i,j,ni,nj;
  uint i_orbitals, j_orbitals;
  uint nuc_i,nuc_j;
  vec_type<double,3> A,B,AmB;
  double ai,aj;
  double dsq,ksi,zeta;
  uint num_terms=0, total_num_terms = 0;
  std::vector<uint> func_code, local_dens;

  std::vector<uint> local2globaldens; // Maps reduced density/Fock -> full density/Fock
  std::vector<scalar_type> dens_values;

  uint term_type_counts[NUM_TERM_TYPES]; // Number of threads for a particular type of term (0 = s-s, 1 = p-s, etc)
  uint term_type_offsets[NUM_TERM_TYPES]; // Offsets into the input arrays for each term type
  term_type_offsets[0] = 0; // s-s starts at 0
  uint i_begin, i_end, j_begin, j_end;
  uint tmp_ind = 0;
  uint s_start = 0, p_start = fortran_vars.s_funcs, d_start = fortran_vars.s_funcs + fortran_vars.p_funcs*3, m = fortran_vars.m;
  uint dd_orb;
  if (forces) { dd_orb = 1; } // 1 for d-d means 6 threads use one of dxx, dyx, etc for func i
  else        { dd_orb = 6; }

  //                                       s-s        p-s        p-p        d-s        d-p        d-d
  uint i_begin_vals[MAX_TERM_TYPE]   = { s_start,   p_start,   p_start,   d_start,   d_start,   d_start};
  uint i_end_vals[MAX_TERM_TYPE]     = { p_start,   d_start,   d_start,   m,         m,         m      };
  uint j_begin_vals[MAX_TERM_TYPE]   = { s_start,   s_start,   p_start,   s_start,   p_start,   d_start};
  uint j_end_vals[MAX_TERM_TYPE]     = { p_start-1, p_start-1, d_start-1, p_start-1, d_start-1, m-1    };
  uint i_orbital_vals[MAX_TERM_TYPE] = { 1,         3,         3,         6,         6,         dd_orb }; 
  uint j_orbital_vals[MAX_TERM_TYPE] = { 1,         1,         3,         1,         3,         6      };

  uint local_dens_ind, num_dens_terms = 0, total_dens_terms = 0;
  uint dens_counts[NUM_TERM_TYPES], dens_offsets[NUM_TERM_TYPES];
  dens_offsets[0] = 0;
  uint tmp_dens_ind = 0;

  Timer nuc,check,prep,kernel,down,reduce;

  //
  // First, the energy/gradient of the nuclear-nuclear interaction between QM and MM centers
  //
  nuc.start_and_sync();
  if (forces) {
    for (i = 0; i < fortran_vars.atoms; i++) {
      double3 qm_pos = fortran_vars.atom_positions(i);
      for (j = 0; j < fortran_vars.clatoms; j++) {
        double3 mm_pos = fortran_vars.clatom_positions(j);
        double3 diff = qm_pos - mm_pos;
        double dist = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
        dist = sqrt(dist);

        double prefactor = -fortran_vars.clatom_charges(j) * (fortran_vars.atom_types(i)+1) / pow(dist,3.0);
        qm_forces[i + 0 * fortran_vars.atoms] += prefactor * diff.x;
        qm_forces[i + 1 * fortran_vars.atoms] += prefactor * diff.y;
        qm_forces[i + 2 * fortran_vars.atoms] += prefactor * diff.z;
        mm_forces[j + 0 * fortran_vars.clatoms] -= prefactor * diff.x;
        mm_forces[j + 1 * fortran_vars.clatoms] -= prefactor * diff.y;
        mm_forces[j + 2 * fortran_vars.clatoms] -= prefactor * diff.z;
      }
    }
  } else {
    for (i = 0; i < fortran_vars.atoms; i++) {
      double3 qm_pos = fortran_vars.atom_positions(i);
      for (j = 0; j < fortran_vars.clatoms; j++) {
        double3 mm_pos = fortran_vars.clatom_positions(j);
        double3 diff = qm_pos - mm_pos;
        double dist = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
        dist = sqrt(dist);

        double E = fortran_vars.clatom_charges(j) * (fortran_vars.atom_types(i)+1) / dist;
        Ens += E;
      }
    }
  }
  nuc.pause();

  //
  // Do check between all basis primitives to find those with significant overlap
  // Check the resulting Gaussian argument from two primitives to the rmax parameter; only use primitives within that cut-off
  //
  // A single thread gets mapped to a pair of significant primitives 
  // We set up here arrays that tell which two functions/two primitives a thread is calculating
  // We also pick out the density matrix elements for significant functions here
  //
  check.start();
  for (uint current_term_type = 0; current_term_type < NUM_TERM_TYPES; current_term_type++) {

    term_type_counts[current_term_type] = 0;
    i_begin = i_begin_vals[current_term_type]; i_end = i_end_vals[current_term_type];
    j_begin = j_begin_vals[current_term_type]; j_end = j_end_vals[current_term_type];
    i_orbitals = i_orbital_vals[current_term_type];
    j_orbitals = j_orbital_vals[current_term_type];

    dens_counts[current_term_type] = 0;
    local_dens_ind = 0;

    // We pad the input arrays between term types, so the offsets for each term type need to be tracked
    if (current_term_type > 0) {
      tmp_ind += COALESCED_DIMENSION(term_type_counts[current_term_type-1]);
      term_type_offsets[current_term_type] = tmp_ind;
      tmp_dens_ind += COALESCED_DIMENSION(dens_counts[current_term_type-1]);
      dens_offsets[current_term_type] = tmp_dens_ind;
    }

    // function i, center A
    for (i = i_begin; i < i_end; i += i_orbitals) {
      nuc_i = fortran_vars.nucleii(i) - 1;
      A = fortran_vars.atom_positions(nuc_i);
      // function j, center B
      for (j = j_begin; j <= ((i > j_end)? j_end : i); j += j_orbitals) {
        nuc_j = fortran_vars.nucleii(j) - 1;
        B = fortran_vars.atom_positions(nuc_j);
        AmB = A - B;
        dsq = length2(AmB);
        bool use_funcs = false; // Do these two functions have any significant primitive pairs?
        // primitive ni, function i
        for (ni = 0; ni < fortran_vars.contractions(i); ni++) {
          // primitive nj, function j
          for (nj = 0; nj < fortran_vars.contractions(j); nj++) {
            ai = fortran_vars.a_values(i,ni);
            aj = fortran_vars.a_values(j,nj);
            zeta = ai + aj;
            ksi = ai * aj / zeta;
            total_num_terms += (i==j)? i_orbitals*(i_orbitals+1)/2 : i_orbitals * j_orbitals;

            if (dsq*ksi < fortran_vars.rmax) {
              use_funcs = true;
              num_terms++;
              term_type_counts[current_term_type]++;

              // Encoding which two primitives this thread will calculate a force term for in one number
              // NOTE: with a long integer, we can use this scheme up to an m of about 9000
              //       if m needs to ever be larger than that, we need to break this up into multiple arrays
              uint this_func_code = nj;                                                         // First, primitive # nj in the lowest part
              this_func_code     += ni * MAX_CONTRACTIONS;                                      // Primitive # ni after the space for nj
              this_func_code     += j  * MAX_CONTRACTIONS * MAX_CONTRACTIONS;                   // Function # j after space for primitives
              this_func_code     += i  * MAX_CONTRACTIONS * MAX_CONTRACTIONS * fortran_vars.m;  // Finally, function # i in the highest part

              func_code.push_back(this_func_code); // Which primitives the thread represents
              local_dens.push_back(local_dens_ind); // Which part of the (reduced) density matrix the thread needs

            }
          }
        }

        total_dens_terms += (i==j)? i_orbitals*(i_orbitals+1)/2 : i_orbitals * j_orbitals;
        // dens_values is a reduced density matrix that only keeps the elements of functions with significant primitive pairs
        // local2globaldens maps from a reduced density/Fock index back to the full matrix
        if (use_funcs) {
          for (uint i_orbital = 0; i_orbital < i_orbitals; i_orbital++) {
            uint j_orbital_finish = (i==j)? i_orbital+1 : j_orbitals;
            for (uint j_orbital = 0; j_orbital < j_orbital_finish; j_orbital++) {
              num_dens_terms++;
              dens_counts[current_term_type]++;

              uint dens_ind = (i+i_orbital) + (2*fortran_vars.m-((j+j_orbital)+1))*(j+j_orbital)/2;
              dens_values.push_back(fortran_vars.rmm_input_ndens1.data[dens_ind]);
              if (!forces) {
                local2globaldens.push_back(dens_ind);
              }
              local_dens_ind++;
            }
          }
        }
      }
    }
    // Pad the input arrays so the next term type has an aligned offset
    for (j = term_type_counts[current_term_type]; j < COALESCED_DIMENSION(term_type_counts[current_term_type]); j++) {
      func_code.push_back(func_code[term_type_offsets[current_term_type]]); // Use the first code from this term type
      local_dens.push_back(local_dens[term_type_offsets[current_term_type]]);
    }
    for (j = dens_counts[current_term_type]; j < COALESCED_DIMENSION(dens_counts[current_term_type]); j++) {
      dens_values.push_back(dens_values[dens_offsets[current_term_type]]);
      if (!forces) {
        local2globaldens.push_back(local2globaldens[dens_offsets[current_term_type]]);
      }
    }
  }
  check.pause();

  //std::cout << "[G2G_QMMM] Number of threads: " << num_terms << std::endl;
  //std::cout << "[G2G_QMMM] Total Gaussian pairs: " << total_num_terms << std::endl;
  //std::cout << "[G2G_QMMM] Number of significant density elements: " << num_dens_terms << std::endl;
  //std::cout << "[G2G_QMMM] Total density elements: " << total_dens_terms << std::endl;

  prep.start();

  // Pad the input so that out-of-range threads do a dummy calculation (same as the first thread), rather than branching and idling
  for (i = 0; i < QMMM_BLOCK_SIZE - (COALESCED_DIMENSION(term_type_counts[NUM_TERM_TYPES-1]) % QMMM_BLOCK_SIZE); i++) {
    func_code.push_back(func_code[term_type_offsets[NUM_TERM_TYPES-1]]);
    local_dens.push_back(local_dens[term_type_offsets[NUM_TERM_TYPES-1]]);
  }
  for (i = 0; i < QMMM_BLOCK_SIZE - (COALESCED_DIMENSION(dens_counts[NUM_TERM_TYPES-1]) % QMMM_BLOCK_SIZE); i++) {
    dens_values.push_back(dens_values[dens_offsets[NUM_TERM_TYPES-1]]);
  }

  HostMatrix<vec_type<scalar_type, 2> > factor_ac_cpu(COALESCED_DIMENSION(fortran_vars.m), MAX_CONTRACTIONS);
  HostMatrixUInt nuc_cpu(fortran_vars.m, 1);
  CudaMatrix<vec_type<scalar_type, 2> > factor_ac_gpu;
  CudaMatrixUInt nuc_gpu;

  //
  // Set up device arrays for function values and mapping function -> nuclei
  //
  // TODO: tried contracting the nuc / function value arrays to a size of total_funcs rather than m, seems to slow down the kernel
  // Doing this requires some more indexing math in the kernel, but need to test on bigger test case
  uint localfunc = 0, func = 0;
  while (func < fortran_vars.m) {
    nuc_cpu(localfunc) = fortran_vars.nucleii(func) - 1;
    for (uint k = 0; k < fortran_vars.contractions(func); k++) {
      factor_ac_cpu(localfunc, k) = vec_type<scalar_type, 2>(fortran_vars.a_values(func, k), fortran_vars.c_values(func, k));
    }
    func++;
    localfunc++;
  }
  factor_ac_gpu = factor_ac_cpu;
  nuc_gpu = nuc_cpu;

  //
  // Send the MM atom positions and charges to the device
  //
  CudaMatrix<vec_type<scalar_type, 3> > clatom_pos_gpu;
  CudaMatrix<scalar_type> clatom_chg_gpu;
  {
    HostMatrix<vec_type<scalar_type,3> > clatom_pos_cpu(fortran_vars.clatoms,1);
    HostMatrix<scalar_type> clatom_chg_cpu(fortran_vars.clatoms,1);
    for (i = 0; i < fortran_vars.clatoms; i++) {
      clatom_pos_cpu(i) = vec_type<scalar_type,3>(fortran_vars.clatom_positions(i).x,fortran_vars.clatom_positions(i).y,fortran_vars.clatom_positions(i).z);
      clatom_chg_cpu(i) = fortran_vars.clatom_charges(i);
    }
    clatom_pos_gpu = clatom_pos_cpu;
    clatom_chg_gpu = clatom_chg_cpu;
  }

  //
  // Send input arrays (thread->primitive map and thread->density map) to the device
  //
  CudaMatrixUInt dev_func_code(func_code), dev_local_dens(local_dens);
  //
  // Send reduced density matrix to the device
  //
  CudaMatrix<scalar_type> dev_dens_values(dens_values);

  //
  // Allocate output arrays on device
  // Currently, each block in the kernel reduces its output, so the output arrays have length (# block)
  // 
  uint partial_out_size = 0, max_partial_size = 0;
  uint out_offsets[NUM_TERM_TYPES];
  out_offsets[0] = 0;
  // Output arrays (probably) don't need to be padded for alignment as only one (or three) threads per block write to them
  for (i = 0; i < NUM_TERM_TYPES; i++) {
    uint this_count = divUp(term_type_counts[i],QMMM_BLOCK_SIZE);
    if (this_count > max_partial_size) max_partial_size = this_count;
    partial_out_size += this_count;
    if (i+1<NUM_TERM_TYPES) { out_offsets[i+1] = partial_out_size; }
  }
  CudaMatrix<vec_type<scalar_type,3> > gpu_partial_mm_forces, gpu_partial_qm_forces;
  CudaMatrix<scalar_type> gpu_partial_fock;
  //
  // Forces: output is partial QM and MM forces
  //
  if (forces)
  {
    gpu_partial_mm_forces.resize(COALESCED_DIMENSION(partial_out_size), fortran_vars.clatoms);
    gpu_partial_qm_forces.resize(COALESCED_DIMENSION(partial_out_size), fortran_vars.atoms);
  //
  // Fock: ouptut is partial Fock elements
  //
  }
  else {
    // The partial Fock matrix is partitioned by term type, so the second (partial) dimension needs to be as big as largest count of a single term type
    gpu_partial_fock.resize(dens_values.size(),max_partial_size);
    dim3 threads(dens_values.size(),max_partial_size);
    dim3 blockSize(32,4);
    dim3 gridSize = divUp(threads,blockSize);
    //
    // Zero the partial Fock matrix on the GPU
    //
    zero_fock<scalar_type><<<gridSize,blockSize>>>(gpu_partial_fock.data,dens_values.size(),max_partial_size);
  }

  //
  // When calculating energies, the energy gets reduced per-block; we figure out the offets/counts of different term types into the partial output energy array here
  //
  uint energies_offsets[NUM_TERM_TYPES];
  uint energies_size = 0;
  CudaMatrix<scalar_type> gpu_qmmm_partial_energies;
  if (!forces) {
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      energies_offsets[i] = energies_size;
      energies_size += divUp(dens_counts[i],QMMM_REDUCE_BLOCK_SIZE);
    }
    gpu_qmmm_partial_energies.resize(energies_size,1);
  }

  //
  // The STR table for F(m,U) calculation is being accessed via texture fetches
  //
  cudaBindTextureToArray(qmmm_str_tex,gammaArray);
  prep.pause_and_sync();

  //
  // Forces kernel
  //
  if (forces) {

    kernel.start();
#define qmmm_forces_parameters \
  term_type_counts[i], factor_ac_gpu.data, nuc_gpu.data, dev_dens_values.data+dens_offset, dev_func_code.data+offset,dev_local_dens.data+offset, \
  gpu_partial_mm_forces.data+force_offset, gpu_partial_qm_forces.data+force_offset, COALESCED_DIMENSION(partial_out_size),clatom_pos_gpu.data,clatom_chg_gpu.data
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_TERM_TYPES];
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (i = 0; i < NUM_TERM_TYPES; i++)
    {
      uint offset = term_type_offsets[i];
      uint dens_offset = dens_offsets[i];
      uint force_offset = out_offsets[i];
      dim3 threads = term_type_counts[i];
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads, blockSize);
      switch (i) {
        case 0: gpu_qmmm_forces<scalar_type,0><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 1: gpu_qmmm_forces<scalar_type,1><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 2: gpu_qmmm_forces<scalar_type,2><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 3: gpu_qmmm_forces<scalar_type,3><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 4: gpu_qmmm_forces<scalar_type,4><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 5: gpu_qmmm_forces<scalar_type,5><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
      }
    }
    cudaDeviceSynchronize();
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }
    kernel.pause();
  }
  //
  // Energy/Fock kernel
  //
  else {

    kernel.start();

#define qmmm_fock_parameters \
  term_type_counts[i], factor_ac_gpu.data, nuc_gpu.data, dev_func_code.data+offset,dev_local_dens.data+offset, /*num_dens_terms,*/ \
  gpu_partial_fock.data+fock_offset, dens_values.size()/*COALESCED_DIMENSION(num_dens_terms)*/,clatom_pos_gpu.data,clatom_chg_gpu.data//,fock_out_offset
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_TERM_TYPES];
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (i = 0; i < NUM_TERM_TYPES; i++)
    {
      uint offset = term_type_offsets[i];
      uint fock_offset = dens_offsets[i];
      dim3 threads = term_type_counts[i];
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads, blockSize);
      switch (i) {
        case 0: gpu_qmmm_fock<scalar_type,0><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 1: gpu_qmmm_fock<scalar_type,1><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 2: gpu_qmmm_fock<scalar_type,2><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 3: gpu_qmmm_fock<scalar_type,3><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 4: gpu_qmmm_fock<scalar_type,4><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 5: gpu_qmmm_fock<scalar_type,5><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
      }
      
      //
      // Reduce the partial Fock terms for a particular term type as soon as that kernel is done; also calculate partial energies for that type
      //
      dim3 reduceThreads = dens_counts[i];
      dim3 reduceBlockSize(QMMM_REDUCE_BLOCK_SIZE);
      dim3 reduceGridSize = divUp(reduceThreads,reduceBlockSize);
      gpu_qmmm_fock_reduce<scalar_type><<<reduceGridSize,reduceBlockSize,0,stream[i]>>>( gpu_partial_fock.data+fock_offset, dev_dens_values.data+fock_offset, gpu_qmmm_partial_energies.data+energies_offsets[i],
                                                                                         dens_values.size(), max_partial_size, dens_counts[i] );
    }
    cudaDeviceSynchronize();
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }
    kernel.pause();
  }

  if (forces) {
    //
    // Download the partial forces from the device 
    // TODO: this could maybe be done asynchronously with the kernels; as one term type finishes we can download its forces, etc
    //
    down.start();
    HostMatrix<vec_type<scalar_type,3> > cpu_partial_mm_forces(gpu_partial_mm_forces), cpu_partial_qm_forces(gpu_partial_qm_forces);
    down.pause_and_sync();

    //
    // Accumulate partial output
    //
    // TODO: need to think about how to accumulate individual force terms
    // Currently, we reduce on a per-block basis in the kernel, then accumulate the block results here on the host
    //
    // The energy partial results are being reduced on-device and that works very well, could probably do that for forces too
    //
    reduce.start();
    for (i = 0; i < fortran_vars.atoms; i++) {
      for (j = 0; j < partial_out_size; j++) {
        qm_forces[i + 0 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).x;
        qm_forces[i + 1 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).y;
        qm_forces[i + 2 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).z;
      }
    }
    for (i = 0; i < fortran_vars.clatoms; i++) {
      for (j = 0; j < partial_out_size; j++) {
        mm_forces[i + 0 * fortran_vars.clatoms] += cpu_partial_mm_forces(j,i).x;
        mm_forces[i + 1 * fortran_vars.clatoms] += cpu_partial_mm_forces(j,i).y;
        mm_forces[i + 2 * fortran_vars.clatoms] += cpu_partial_mm_forces(j,i).z;
      }
    }
    reduce.pause();
  } else {
    //
    // Download reduced Fock matrix and partially reduced energies
    // The Fock matrix has been reduced to the first row of the output, so we only want that much of the device array
    //
    down.start();
    HostMatrix<scalar_type> cpu_fock(dens_values.size());
    cudaMemcpy(cpu_fock.data,gpu_partial_fock.data,cpu_fock.bytes(),cudaMemcpyDeviceToHost);
    HostMatrix<scalar_type> cpu_partial_energies(gpu_qmmm_partial_energies);
    down.pause_and_sync();

    //
    // Send Fock elements back to RMM(M11) and do final reduction of e-nuc energies into Es
    //
    reduce.start();
    for (uint t = 0; t < NUM_TERM_TYPES; t++) {
      for (i = dens_offsets[t]; i < dens_offsets[t] + dens_counts[t]; i++) {
        uint dens_ind = local2globaldens[i];
        fortran_vars.rmm_1e_output(dens_ind) += cpu_fock(i);
      }
    }
    for (i = 0; i < energies_size; i++) {
      Es += cpu_partial_energies(i);
    }
    reduce.pause();
  }

  //cout << "[G2G_QMMM] nuc-nuc: " << nuc << " overlap check: " << check << " kernel prep: " << prep << endl;
  //cout << "[G2G_QMMM] kernel: " << kernel << " download: " << down << " host reduction: " << reduce << endl;

  cudaUnbindTexture(qmmm_str_tex);

  cudaAssertNoError("qmmm");
}

template<class scalar_type>
void clean_gamma( void ) {
  //scalar_type* d_str_ptr;
  //cudaMemcpyFromSymbol(&d_str_ptr,gpu_str,sizeof(d_str_ptr));
  //cudaFree(d_str_ptr);
  cudaFreeArray(gammaArray);

  cudaAssertNoError("clean_gamma");
}
#if FULL_DOUBLE
template void g2g_qmmm<double,true>(double* qm_forces, double* mm_forces, double& Ens, double& Es);
template void g2g_qmmm<double,false>(double* qm_forces, double* mm_forces, double& Ens, double& Es);
template void clean_gamma<double>( void );
#else
template void g2g_qmmm<float,true>(double* qm_forces, double* mm_forces, double& Ens, double& Es);
template void g2g_qmmm<float,false>(double* qm_forces, double* mm_forces, double& Ens, double& Es);
template void clean_gamma<float>( void );
#endif

#define NUM_COULOMB_TERM_TYPES 6 // 6 types when using s,p,and d functions: s-s,p-s,p-p,d-s,d-p,d-d

//
// Main QM/MM routine
// If forces = true, calculate gradients of QM/MM operator and return in qm_forces and mm_forces
// If forces = false, calculate Fock matrix elements of QM/MM operator (returned in RMM(M11) back in Fortran)
//                    and QM/MM energy of the current density (nuc-nuc in Ens, e-nuc in Es)
//
template <class scalar_type,bool forces> void g2g_coulomb(double* qm_forces, double& Es)
{
  uint i,j,ni,nj;
  uint i_orbitals, j_orbitals;
  uint nuc_i,nuc_j;
  vec_type<double,3> A,B,AmB;
  double ai,aj;
  double dsq,ksi,zeta;
  uint num_terms=0, total_num_terms = 0;
  std::vector<uint> func_code, local_dens;

  //std::vector<uint> local2globaldens; // Maps reduced density/Fock -> full density/Fock
  std::vector<scalar_type> dens_values;

  uint term_type_counts[NUM_COULOMB_TERM_TYPES]; // Number of threads for a particular type of term (0 = s-s, 1 = p-s, etc)
  uint term_type_offsets[NUM_COULOMB_TERM_TYPES]; // Offsets into the input arrays for each term type
  term_type_offsets[0] = 0; // s-s starts at 0
  uint i_begin, i_end, j_begin, j_end;
  uint tmp_ind = 0;
  uint s_start = 0, p_start = fortran_vars.s_funcs, d_start = fortran_vars.s_funcs + fortran_vars.p_funcs*3, m = fortran_vars.m;
  uint dd_orb;
  if (forces) { dd_orb = 6; } // 1 for d-d means 6 threads use one of dxx, dyx, etc for func i
  else        { dd_orb = 6; }

  //                                       s-s        p-s        p-p        d-s        d-p        d-d
  uint i_begin_vals[MAX_TERM_TYPE]   = { s_start,   p_start,   p_start,   d_start,   d_start,   d_start};
  uint i_end_vals[MAX_TERM_TYPE]     = { p_start,   d_start,   d_start,   m,         m,         m      };
  uint j_begin_vals[MAX_TERM_TYPE]   = { s_start,   s_start,   p_start,   s_start,   p_start,   d_start};
  uint j_end_vals[MAX_TERM_TYPE]     = { p_start-1, p_start-1, d_start-1, p_start-1, d_start-1, m-1    };
  uint i_orbital_vals[MAX_TERM_TYPE] = { 1,         3,         3,         6,         6,         dd_orb }; 
  uint j_orbital_vals[MAX_TERM_TYPE] = { 1,         1,         3,         1,         3,         6      };

  uint local_dens_ind, num_dens_terms = 0, total_dens_terms = 0;
  uint dens_counts[NUM_COULOMB_TERM_TYPES], dens_offsets[NUM_COULOMB_TERM_TYPES];
  dens_offsets[0] = 0;
  uint tmp_dens_ind = 0;

  Timer check,prep,kernel,down,reduce;

  //
  // Do check between all basis primitives to find those with significant overlap
  // Check the resulting Gaussian argument from two primitives to the rmax parameter; only use primitives within that cut-off
  //
  // A single thread gets mapped to a pair of significant primitives 
  // We set up here arrays that tell which two functions/two primitives a thread is calculating
  // We also pick out the density matrix elements for significant functions here
  //
  check.start();
  for (uint current_term_type = 0; current_term_type < NUM_COULOMB_TERM_TYPES; current_term_type++) {

    term_type_counts[current_term_type] = 0;
    i_begin = i_begin_vals[current_term_type]; i_end = i_end_vals[current_term_type];
    j_begin = j_begin_vals[current_term_type]; j_end = j_end_vals[current_term_type];
    i_orbitals = i_orbital_vals[current_term_type];
    j_orbitals = j_orbital_vals[current_term_type];

    dens_counts[current_term_type] = 0;
    local_dens_ind = 0;

    // We pad the input arrays between term types, so the offsets for each term type need to be tracked
    if (current_term_type > 0) {
      tmp_ind += COALESCED_DIMENSION(term_type_counts[current_term_type-1]);
      term_type_offsets[current_term_type] = tmp_ind;
      tmp_dens_ind += COALESCED_DIMENSION(dens_counts[current_term_type-1]);
      dens_offsets[current_term_type] = tmp_dens_ind;
    }

    // function i, center A
    for (i = i_begin; i < i_end; i += i_orbitals) {
      nuc_i = fortran_vars.nucleii(i) - 1;
      A = fortran_vars.atom_positions(nuc_i);
      // function j, center B
      for (j = j_begin; j <= ((i > j_end)? j_end : i); j += j_orbitals) {
        nuc_j = fortran_vars.nucleii(j) - 1;
        B = fortran_vars.atom_positions(nuc_j);
        AmB = A - B;
        dsq = length2(AmB);
        bool use_funcs = false; // Do these two functions have any significant primitive pairs?
        // primitive ni, function i
        for (ni = 0; ni < fortran_vars.contractions(i); ni++) {
          // primitive nj, function j
          for (nj = 0; nj < fortran_vars.contractions(j); nj++) {
            ai = fortran_vars.a_values(i,ni);
            aj = fortran_vars.a_values(j,nj);
            zeta = ai + aj;
            ksi = ai * aj / zeta;
            total_num_terms += (i==j)? i_orbitals*(i_orbitals+1)/2 : i_orbitals * j_orbitals;

            if (dsq*ksi < fortran_vars.rmax) {
              use_funcs = true;
              num_terms++;
              term_type_counts[current_term_type]++;

              // Encoding which two primitives this thread will calculate a force term for in one number
              // NOTE: with a long integer, we can use this scheme up to an m of about 9000
              //       if m needs to ever be larger than that, we need to break this up into multiple arrays
              uint this_func_code = nj;                                                         // First, primitive # nj in the lowest part
              this_func_code     += ni * MAX_CONTRACTIONS;                                      // Primitive # ni after the space for nj
              this_func_code     += j  * MAX_CONTRACTIONS * MAX_CONTRACTIONS;                   // Function # j after space for primitives
              this_func_code     += i  * MAX_CONTRACTIONS * MAX_CONTRACTIONS * fortran_vars.m;  // Finally, function # i in the highest part

              func_code.push_back(this_func_code); // Which primitives the thread represents
              local_dens.push_back(local_dens_ind); // Which part of the (reduced) density matrix the thread needs

            }
          }
        }

        total_dens_terms += (i==j)? i_orbitals*(i_orbitals+1)/2 : i_orbitals * j_orbitals;
        // dens_values is a reduced density matrix that only keeps the elements of functions with significant primitive pairs
        // local2globaldens maps from a reduced density/Fock index back to the full matrix
        if (use_funcs) {
          for (uint i_orbital = 0; i_orbital < i_orbitals; i_orbital++) {
            uint j_orbital_finish = (i==j)? i_orbital+1 : j_orbitals;
            for (uint j_orbital = 0; j_orbital < j_orbital_finish; j_orbital++) {
              num_dens_terms++;
              dens_counts[current_term_type]++;

              uint dens_ind = (i+i_orbital) + (2*fortran_vars.m-((j+j_orbital)+1))*(j+j_orbital)/2;
              dens_values.push_back(fortran_vars.rmm_input_ndens1.data[dens_ind]);
              if (!forces) {
              //  local2globaldens.push_back(dens_ind);
              }
              local_dens_ind++;
            }
          }
        }
      }
    }
    // Pad the input arrays so the next term type has an aligned offset
    for (j = term_type_counts[current_term_type]; j < COALESCED_DIMENSION(term_type_counts[current_term_type]); j++) {
      func_code.push_back(func_code[term_type_offsets[current_term_type]]); // Use the first code from this term type
      local_dens.push_back(local_dens[term_type_offsets[current_term_type]]);
    }
    for (j = dens_counts[current_term_type]; j < COALESCED_DIMENSION(dens_counts[current_term_type]); j++) {
      dens_values.push_back(dens_values[dens_offsets[current_term_type]]);
      if (!forces) {
      //  local2globaldens.push_back(local2globaldens[dens_offsets[current_term_type]]);
      }
    }
  }
  check.pause();

  std::cout << "[G2G_COULOMB] Number of threads: " << num_terms << std::endl;
  std::cout << "[G2G_COULOMB] Total Gaussian pairs: " << total_num_terms << std::endl;
  std::cout << "[G2G_COULOMB] Number of significant density elements: " << num_dens_terms << std::endl;
  std::cout << "[G2G_COULOMB] Total density elements: " << total_dens_terms << std::endl;

  prep.start();

  // Pad the input so that out-of-range threads do a dummy calculation (same as the first thread), rather than branching and idling
  for (i = 0; i < QMMM_BLOCK_SIZE - (COALESCED_DIMENSION(term_type_counts[NUM_COULOMB_TERM_TYPES-1]) % QMMM_BLOCK_SIZE); i++) {
    func_code.push_back(func_code[term_type_offsets[NUM_COULOMB_TERM_TYPES-1]]);
    local_dens.push_back(local_dens[term_type_offsets[NUM_COULOMB_TERM_TYPES-1]]);
  }
  for (i = 0; i < QMMM_BLOCK_SIZE - (COALESCED_DIMENSION(dens_counts[NUM_COULOMB_TERM_TYPES-1]) % QMMM_BLOCK_SIZE); i++) {
    dens_values.push_back(dens_values[dens_offsets[NUM_COULOMB_TERM_TYPES-1]]);
  }

  HostMatrix<vec_type<scalar_type, 2> > factor_ac_cpu(COALESCED_DIMENSION(fortran_vars.m), MAX_CONTRACTIONS);
  HostMatrixUInt nuc_cpu(fortran_vars.m, 1);
  CudaMatrix<vec_type<scalar_type, 2> > factor_ac_gpu;
  CudaMatrixUInt nuc_gpu;

  //
  // Set up device arrays for function values and mapping function -> nuclei
  //
  // TODO: tried contracting the nuc / function value arrays to a size of total_funcs rather than m, seems to slow down the kernel
  // Doing this requires some more indexing math in the kernel, but need to test on bigger test case
  uint localfunc = 0, func = 0;
  while (func < fortran_vars.m) {
    nuc_cpu(localfunc) = fortran_vars.nucleii(func) - 1;
    for (uint k = 0; k < fortran_vars.contractions(func); k++) {
      factor_ac_cpu(localfunc, k) = vec_type<scalar_type, 2>(fortran_vars.a_values(func, k), fortran_vars.c_values(func, k));
    }
    func++;
    localfunc++;
  }
  factor_ac_gpu = factor_ac_cpu;
  nuc_gpu = nuc_cpu;

  std::vector<vec_type<scalar_type,2> > factor_ac_dens_cpu;
  std::vector<scalar_type> fit_dens_cpu;
  std::vector<vec_type<scalar_type,3> > nuc_dens_cpu;
  std::vector<uint> nuc_ind_dens_cpu;
  //HostMatrix<scalar_type> fit_dens_factor_c_cpu(fortran_vars.);
  for (func = 0; func < fortran_vars.m_dens; func++) {
    uint nuc_ind = fortran_vars.nucleii_dens(func) - 1;
    double3 nuc_pos = fortran_vars.atom_positions(nuc_ind);
    for (uint k = 0; k < fortran_vars.contractions_dens(func); k++) {
      factor_ac_dens_cpu.push_back(vec_type<scalar_type,2>(fortran_vars.a_values_dens(func,k),fortran_vars.c_values_dens(func,k)));
      nuc_dens_cpu.push_back(vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
      nuc_ind_dens_cpu.push_back(nuc_ind);
      fit_dens_cpu.push_back(fortran_vars.af_input_ndens1(func));
    }
  }
  CudaMatrix<vec_type<scalar_type,2> > factor_ac_dens_gpu(factor_ac_dens_cpu);
  CudaMatrix<vec_type<scalar_type, 3> > nuc_dens_gpu(nuc_dens_cpu);
  CudaMatrix<scalar_type> fit_dens_gpu(fit_dens_cpu);
  CudaMatrix<uint> nuc_ind_dens_gpu(nuc_ind_dens_cpu);

  //
  // Send input arrays (thread->primitive map and thread->density map) to the device
  //
  CudaMatrixUInt dev_func_code(func_code), dev_local_dens(local_dens);
  //
  // Send reduced density matrix to the device
  //
  CudaMatrix<scalar_type> dev_dens_values(dens_values);

  //
  // Allocate output arrays on device
  // Currently, each block in the kernel reduces its output, so the output arrays have length (# block)
  // 
  uint partial_out_size = 0;//, max_partial_size = 0;
  uint out_offsets[NUM_COULOMB_TERM_TYPES];
  // Output arrays (probably) don't need to be padded for alignment as only one (or three) threads per block write to them
  for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
    uint this_count = divUp(term_type_counts[i],QMMM_BLOCK_SIZE);
    //if (this_count > max_partial_size) max_partial_size = this_count;
    out_offsets[i] = partial_out_size;
    partial_out_size += this_count;
  }
  CudaMatrix<vec_type<scalar_type,3> > gpu_partial_qm_forces;
  //CudaMatrix<scalar_type> gpu_partial_fock;
  //
  // Forces: output is partial QM and MM forces
  //
  if (forces)
  {
    gpu_partial_qm_forces.resize(COALESCED_DIMENSION(partial_out_size), fortran_vars.atoms);
    dim3 threads(COALESCED_DIMENSION(partial_out_size),fortran_vars.atoms);
    dim3 blockSize(32,4);
    dim3 gridSize = divUp(threads,blockSize);
    zero_forces<scalar_type><<<gridSize,blockSize>>>(gpu_partial_qm_forces.data,COALESCED_DIMENSION(partial_out_size),fortran_vars.atoms);
    //cudaMemset(gpu_partial_qm_forces.data, 0, COALESCED_DIMENSION(partial_out_size) * fortran_vars.atoms * sizeof(vec_type<scalar_type,3>));

  //
  // Fock: ouptut is partial Fock elements
  //
  }
  /*else {
    // The partial Fock matrix is partitioned by term type, so the second (partial) dimension needs to be as big as largest count of a single term type
    gpu_partial_fock.resize(dens_values.size(),max_partial_size);
    dim3 threads(dens_values.size(),max_partial_size);
    dim3 blockSize(32,4);
    dim3 gridSize = divUp(threads,blockSize);
    //
    // Zero the partial Fock matrix on the GPU
    //
    zero_fock<scalar_type><<<gridSize,blockSize>>>(gpu_partial_fock.data,dens_values.size(),max_partial_size);
  }*/

  //
  // When calculating energies, the energy gets reduced per-block; we figure out the offets/counts of different term types into the partial output energy array here
  //
  /*uint energies_offsets[NUM_COULOMB_TERM_TYPES];
  uint energies_size = 0;
  CudaMatrix<scalar_type> gpu_qmmm_partial_energies;
  if (!forces) {
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      energies_offsets[i] = energies_size;
      energies_size += divUp(dens_counts[i],QMMM_REDUCE_BLOCK_SIZE);
    }
    gpu_qmmm_partial_energies.resize(energies_size,1);
  }*/

  //
  // The STR table for F(m,U) calculation is being accessed via texture fetches
  //
  cudaBindTextureToArray(qmmm_str_tex,gammaArray);
  prep.pause_and_sync();

  //
  // Forces kernel
  //
  if (forces) {

    kernel.start();
#define coulomb_forces_parameters \
  term_type_counts[i], factor_ac_gpu.data, nuc_gpu.data, dev_dens_values.data+dens_offset, dev_func_code.data+offset,dev_local_dens.data+offset, \
  gpu_partial_qm_forces.data+force_offset, COALESCED_DIMENSION(partial_out_size),factor_ac_dens_gpu.data,nuc_dens_gpu.data,nuc_ind_dens_gpu.data,fit_dens_gpu.data //,clatom_pos_gpu.data,clatom_chg_gpu.data
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_COULOMB_TERM_TYPES];
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++)
    {
      uint offset = term_type_offsets[i];
      uint dens_offset = dens_offsets[i];
      uint force_offset = out_offsets[i];
      dim3 threads = term_type_counts[i];
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads, blockSize);
      switch (i) {
        case 0: gpu_coulomb_forces<scalar_type,0><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 1: gpu_coulomb_forces<scalar_type,1><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 2: gpu_coulomb_forces<scalar_type,2><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 3: gpu_coulomb_forces<scalar_type,3><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 4: gpu_coulomb_forces<scalar_type,4><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 5: gpu_coulomb_forces<scalar_type,5><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
      }
    }
    cudaDeviceSynchronize();
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }
    kernel.pause();
  }
  //
  // Energy/Fock kernel
  //
  /*else {

    kernel.start();

#define qmmm_fock_parameters \
  term_type_counts[i], factor_ac_gpu.data, nuc_gpu.data, dev_func_code.data+offset,dev_local_dens.data+offset, \
  gpu_partial_fock.data+fock_offset, dens_values.size(),clatom_pos_gpu.data,clatom_chg_gpu.data//,fock_out_offset
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_COULOMB_TERM_TYPES];
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++)
    {
      uint offset = term_type_offsets[i];
      uint fock_offset = dens_offsets[i];
      dim3 threads = term_type_counts[i];
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads, blockSize);
      switch (i) {
        case 0: gpu_qmmm_fock<scalar_type,0><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 1: gpu_qmmm_fock<scalar_type,1><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 2: gpu_qmmm_fock<scalar_type,2><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 3: gpu_qmmm_fock<scalar_type,3><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 4: gpu_qmmm_fock<scalar_type,4><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 5: gpu_qmmm_fock<scalar_type,5><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
      }
      
      //
      // Reduce the partial Fock terms for a particular term type as soon as that kernel is done; also calculate partial energies for that type
      //
      dim3 reduceThreads = dens_counts[i];
      dim3 reduceBlockSize(QMMM_REDUCE_BLOCK_SIZE);
      dim3 reduceGridSize = divUp(reduceThreads,reduceBlockSize);
      gpu_qmmm_fock_reduce<scalar_type><<<reduceGridSize,reduceBlockSize,0,stream[i]>>>( gpu_partial_fock.data+fock_offset, dev_dens_values.data+fock_offset, gpu_qmmm_partial_energies.data+energies_offsets[i],
                                                                                         dens_values.size(), max_partial_size, dens_counts[i] );
    }
    cudaDeviceSynchronize();
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }
    kernel.pause();
  }*/

  if (forces) {
    //
    // Download the partial forces from the device 
    // TODO: this could maybe be done asynchronously with the kernels; as one term type finishes we can download its forces, etc
    //
    down.start();
    HostMatrix<vec_type<scalar_type,3> > cpu_partial_qm_forces(gpu_partial_qm_forces);
    down.pause_and_sync();

    //
    // Accumulate partial output
    //
    // TODO: need to think about how to accumulate individual force terms
    // Currently, we reduce on a per-block basis in the kernel, then accumulate the block results here on the host
    //
    // The energy partial results are being reduced on-device and that works very well, could probably do that for forces too
    //
    reduce.start();
    for (i = 0; i < fortran_vars.atoms; i++) {
      for (j = 0; j < partial_out_size; j++) {
        qm_forces[i + 0 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).x;
        qm_forces[i + 1 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).y;
        qm_forces[i + 2 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).z;
      }
    }
    reduce.pause();
  } /*else {
    //
    // Download reduced Fock matrix and partially reduced energies
    // The Fock matrix has been reduced to the first row of the output, so we only want that much of the device array
    //
    down.start();
    HostMatrix<scalar_type> cpu_fock(dens_values.size());
    cudaMemcpy(cpu_fock.data,gpu_partial_fock.data,cpu_fock.bytes(),cudaMemcpyDeviceToHost);
    HostMatrix<scalar_type> cpu_partial_energies(gpu_qmmm_partial_energies);
    down.pause_and_sync();

    //
    // Send Fock elements back to RMM(M11) and do final reduction of e-nuc energies into Es
    //
    reduce.start();
    for (uint t = 0; t < NUM_COULOMB_TERM_TYPES; t++) {
      for (i = dens_offsets[t]; i < dens_offsets[t] + dens_counts[t]; i++) {
        uint dens_ind = local2globaldens[i];
        fortran_vars.rmm_1e_output(dens_ind) += cpu_fock(i);
      }
    }
    for (i = 0; i < energies_size; i++) {
      Es += cpu_partial_energies(i);
    }
    reduce.pause();
  }*/

  cout << "[G2G_COULOMB] overlap check: " << check << " kernel prep: " << prep << endl;
  cout << "[G2G_COULOMB] kernel: " << kernel << " download: " << down << " host reduction: " << reduce << endl;

  cudaUnbindTexture(qmmm_str_tex);

  cudaAssertNoError("qmmm");
}

#if FULL_DOUBLE
template void g2g_coulomb<double,true>(double* qm_forces, double& Es);
template void g2g_coulomb<double,false>(double* qm_forces, double& Es);
#else
template void g2g_coulomb<float,true>(double* qm_forces, double& Es);
template void g2g_coulomb<float,false>(double* qm_forces, double& Es);
#endif

}
