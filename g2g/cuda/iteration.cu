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
#else
texture<float, 2, cudaReadModeElementType> rmm_input_gpu_tex;
texture<float, 2, cudaReadModeElementType> rmm_input_gpu_tex2;
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

using std::cout;
using std::vector;
using std::endl;

// Host function to set the constant
void gpu_set_variables(void) {
  cudaMemcpyToSymbol(gpu_normalization_factor, &fortran_vars.normalization_factor, sizeof(fortran_vars.normalization_factor), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(gpu_atoms, &fortran_vars.atoms, sizeof(fortran_vars.atoms), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(gpu_Iexch, &fortran_vars.iexch, sizeof(fortran_vars.iexch), 0, cudaMemcpyHostToDevice);

  cudaAssertNoError("set_gpu_variables");
}

template<class T> void gpu_set_atom_positions(const HostMatrix<T>& m) {
  cudaMemcpyToSymbol(gpu_atom_positions, m.data, m.bytes(), 0, cudaMemcpyHostToDevice);
}

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

}
