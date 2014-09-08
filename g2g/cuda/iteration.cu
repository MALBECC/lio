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
#else
texture<float, 2, cudaReadModeElementType> rmm_input_gpu_tex;
#endif
/** KERNELS **/
#include "gpu_variables.h"
#include "kernels/pot.h"
#include "kernels/accumulate_point.h"
#include "kernels/energy.h"
#include "kernels/energy_derivs.h"
#include "kernels/rmm.h"
#include "kernels/weight.h"
#include "kernels/functions.h"
#include "kernels/force.h"
#include "kernels/transpose.h"

using std::cout;
using std::vector;
using std::endl;

//Definicion global para poder tener acceso



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

template void gpu_set_atom_positions<float3>(const HostMatrix<float3>& m);
template void gpu_set_atom_positions<double3>(const HostMatrix<double3>& m);

template<class scalar_type>
void PointGroup<scalar_type>::solve(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, bool compute_energy,
                                    double& energy, double* fort_forces_ptr)
{
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

  const int block_height= divUp(group_m,2*DENSITY_BLOCK_SIZE);

  threadBlock = dim3(DENSITY_BLOCK_SIZE,1,1); // Hay que asegurarse que la cantidad de funciones este en rango
  threadGrid = dim3(number_of_points,block_height,1);

  CudaMatrix<scalar_type> factors_gpu;

  CudaMatrix<scalar_type> partial_densities_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dxyz_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd1_gpu;
  CudaMatrix<vec_type<scalar_type,4> > dd2_gpu;

  /*
   **********************************************************************
   * Transposiciones de matrices para la coalescencia mejorada en density
   **********************************************************************
   */


  CudaMatrix<scalar_type>   function_values_transposed_gpu;
  CudaMatrix<vec_type<scalar_type,4> > gradient_values_transposed_gpu;
  CudaMatrix<vec_type<scalar_type,4> > hessian_values_transposed_gpu;

  int transposed_width = COALESCED_DIMENSION(number_of_points);

  function_values_transposed_gpu.resize(group_m, COALESCED_DIMENSION(number_of_points));
  if (fortran_vars.do_forces || fortran_vars.gga)
      gradient_values_transposed_gpu.resize( group_m,COALESCED_DIMENSION(number_of_points));
  if (fortran_vars.gga)
      hessian_values_transposed_gpu.resize((group_m) * 2, COALESCED_DIMENSION(number_of_points));

  #define BLOCK_DIM 16
  dim3 transpose_grid(transposed_width / BLOCK_DIM, divUp((group_m),BLOCK_DIM));
  dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);

  transpose<<<transpose_grid, transpose_threads>>> (function_values_transposed_gpu.data, function_values.data,  COALESCED_DIMENSION(number_of_points),group_m   );

  if (fortran_vars.do_forces || fortran_vars.gga)
    transpose_vec<<<transpose_grid, transpose_threads>>> (gradient_values_transposed_gpu.data, gradient_values.data, COALESCED_DIMENSION(number_of_points), group_m );

  transpose_grid=dim3(transposed_width / BLOCK_DIM, divUp((group_m)*2, BLOCK_DIM), 1);

  if (fortran_vars.gga)
    transpose_vec<<<transpose_grid, transpose_threads>>> (hessian_values_transposed_gpu.data, hessian_values.data, COALESCED_DIMENSION(number_of_points), (group_m)*2);


  partial_densities_gpu.resize(COALESCED_DIMENSION(number_of_points), block_height);
  dxyz_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height);
  dd1_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );
  dd2_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );

  const dim3 threadGrid_accumulate(divUp(number_of_points,DENSITY_ACCUM_BLOCK_SIZE),1,1);
  const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE,1,1);

  if (compute_rmm || compute_forces) factors_gpu.resize(number_of_points);

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
  cudaMallocArray(&cuArray, &rmm_input_gpu_tex.channelDesc, rmm_input_cpu.width,rmm_input_cpu.height);
#if FULL_DOUBLE
  cudaMemcpyToArray(cuArray, 0, 0,rmm_input_cpu.data,sizeof(int2)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);
#else
  cudaMemcpyToArray(cuArray, 0, 0,rmm_input_cpu.data,sizeof(float)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);
#endif
  cudaBindTextureToArray(rmm_input_gpu_tex, cuArray);

  rmm_input_gpu_tex.normalized = false;

  if (compute_energy) {
    CudaMatrix<scalar_type> energy_gpu(number_of_points);
#define compute_parameters \
    energy_gpu.data,factors_gpu.data,point_weights_gpu.data,number_of_points,function_values_transposed_gpu.data,gradient_values_transposed_gpu.data,hessian_values_transposed_gpu.data,group_m,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data
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
    NULL,factors_gpu.data,point_weights_gpu.data,number_of_points,function_values_transposed_gpu.data,gradient_values_transposed_gpu.data,hessian_values_transposed_gpu.data,group_m,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data
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

  function_values_transposed_gpu.deallocate();
  gradient_values_transposed_gpu.deallocate();
  hessian_values_transposed_gpu.deallocate();

  timers.density.pause_and_sync();

//************ Repongo los valores que puse a cero antes, para las fuerzas son necesarios (o por lo mens utiles)
  for (uint i=0; i<(group_m); i++)
  {
    for(uint j=0; j<(group_m); j++)
    {
      if((i>=group_m) || (j>=group_m) || (j > i))
      {
        rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*j+i] ;
      }
    }
  }
#if FULL_DOUBLE
  cudaMemcpyToArray(cuArray, 0, 0,rmm_input_cpu.data,sizeof(int2)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);
#else
  cudaMemcpyToArray(cuArray, 0, 0,rmm_input_cpu.data,sizeof(float)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);
#endif

//**********************************************

   dim3 threads;
  /* compute forces */
  if (compute_forces) {
    timers.density_derivs.start_and_sync();
    threads = dim3(number_of_points);
    threadBlock = dim3(DENSITY_DERIV_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);

    CudaMatrix<vec_type4> dd_gpu(COALESCED_DIMENSION(number_of_points), total_nucleii()); dd_gpu.zero();
    CudaMatrixUInt nuc_gpu(func2local_nuc);  // TODO: esto en realidad se podria guardar una sola vez durante su construccion

    gpu_compute_density_derivs<<<threadGrid, threadBlock>>>(function_values.data, gradient_values.data, nuc_gpu.data, dd_gpu.data, number_of_points, group_m, total_nucleii());
    cudaAssertNoError("density_derivs");
    timers.density_derivs.pause_and_sync();

    timers.forces.start_and_sync();
    CudaMatrix<vec_type4> forces_gpu(total_nucleii());

    threads = dim3(total_nucleii());
    threadBlock = dim3(FORCE_BLOCK_SIZE);
    threadGrid = divUp(threads, threadBlock);
    gpu_compute_forces<<<threadGrid, threadBlock>>>(number_of_points, factors_gpu.data, dd_gpu.data, forces_gpu.data, total_nucleii());
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

  /* compute RMM */
  timers.rmm.start_and_sync();
  if (compute_rmm) {
    threads = dim3(group_m, group_m);
    threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);
    threadGrid = divUp(threads, threadBlock);

    CudaMatrix<scalar_type> rmm_output_gpu(COALESCED_DIMENSION(group_m), group_m);
    gpu_update_rmm<<<threadGrid, threadBlock>>>(factors_gpu.data, number_of_points, rmm_output_gpu.data, function_values.data, group_m);
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
    hessian_values.deallocate();
  }
  //Deshago el bind de textura de rmm
  cudaUnbindTexture(rmm_input_gpu_tex); //Enroque el Unbind con el Free, asi parece mas logico. Nano
  cudaFreeArray(cuArray);
}

