#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math_constants.h>
#include <string>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>

#include "../../../g2g/common.h"
#include "../../../g2g/init.h"
#include "../../../g2g/cuda/cuda_extra.h"
#include "../../../g2g/matrix.h"
#include "../../../g2g/timer.h"
#include "../../../g2g/partition.h"
#include "../../../g2g/scalar_vector_types.h"
#include "../../../g2g/global_memory_pool.h"

#include "../../../g2g/libxc/libxcproxy.h"

void solve_closed( bool compute_rmm, bool lda, bool compute_forces, bool compute_energy,
    double& energy,    G2G::HostMatrix<double>& fort_forces_ms,
    int inner_threads, G2G::HostMatrix<double>& rmm_output_local) 
{

//  int device;
//  cudaGetDevice(&device);
//  current_device = device;

  /*** Computo sobre cada cubo ****/
//  CudaMatrix<float> point_weights_gpu;

  /** Compute this group's functions **/
//  timers.functions.start_and_sync();
//  compute_functions(compute_forces, !lda);
//  timers.functions.pause_and_sync();

//  uint group_m = this->total_functions();
//  unit group_m = 10;
//  int number_of_points = 10;

//  timers.density.start_and_sync();
  /** Load points from group **/
//  HostMatrix<float> point_weights_cpu(number_of_points, 1);
//  HostMatrix<float> point_weights_cpu(number_of_points, 1);


//  uint i = 0;
//  for (vector<Point>::const_iterator p = this->points.begin(); p != this->points.end(); ++p, ++i) {
//    point_weights_cpu(i) = p->weight;
//  }
//  point_weights_gpu = point_weights_cpu;

//  dim3 threadBlock, threadGrid;
  /* compute density/factors */

//  const int block_height= divUp(group_m, 2*DENSITY_BLOCK_SIZE);

//  threadBlock = dim3(DENSITY_BLOCK_SIZE,1,1); // Hay que asegurarse que la cantidad de funciones este en rango
//  threadGrid = dim3(number_of_points,block_height,1);

//  CudaMatrix<float> partial_densities_gpu;
//  CudaMatrix< vec_type<float,4> > dxyz_gpu;
//  CudaMatrix< vec_type<float,4> > dd1_gpu;
//  CudaMatrix< vec_type<float,4> > dd2_gpu;

//  partial_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
//  dxyz_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height);
//  dd1_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );
//  dd2_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );

//  const dim3 threadGrid_accumulate(divUp(number_of_points,DENSITY_ACCUM_BLOCK_SIZE),1,1);
//  const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE,1,1);

//  CudaMatrix<float> factors_gpu;
//  if (compute_rmm || compute_forces)
//    factors_gpu.resize(number_of_points);

//  int transposed_width = COALESCED_DIMENSION(number_of_points);
  #define BLOCK_DIM 16
//  dim3 transpose_grid(transposed_width / BLOCK_DIM, divUp((group_m),BLOCK_DIM), 1);
//  dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);

//  CudaMatrix<float> function_values_transposed;
//  CudaMatrix<vec_type<float,4> > gradient_values_transposed;

  // Probar si esta intercalado al pedo.
//  function_values_transposed.resize(group_m, COALESCED_DIMENSION(this->number_of_points));

//  if (fortran_vars.do_forces || fortran_vars.gga)
//      gradient_values_transposed.resize( group_m,COALESCED_DIMENSION(this->number_of_points));

//  transpose<<<transpose_grid, transpose_threads>>> (function_values_transposed.data,
//      function_values.data, COALESCED_DIMENSION(this->number_of_points), group_m);

//  if (fortran_vars.do_forces || fortran_vars.gga)
//    transpose<<<transpose_grid, transpose_threads>>> (gradient_values_transposed.data,
//        gradient_values.data, COALESCED_DIMENSION(this->number_of_points), group_m );
  // fin intercalado al pedo

//  HostMatrix<float> rmm_input_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
//  get_rmm_input(rmm_input_cpu); //Achica la matriz densidad a la version reducida del grupo

//  for (uint i=0; i<(group_m+DENSITY_BLOCK_SIZE); i++)
//  {
//    for(uint j=0; j<COALESCED_DIMENSION(group_m); j++)
//    {
//      if((i>=group_m) || (j>=group_m) || (j > i))
//      {
//        rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
//      }
//    }
//  }


  /*
   **********************************************************************
   * Pasando RDM (rmm) a texturas
   **********************************************************************
   */
//  cudaArray* cuArray;
//  cudaMallocArray(&cuArray, &rmm_input_gpu_tex.channelDesc, rmm_input_cpu.width, rmm_input_cpu.height);
//  cudaMemcpyToArray(cuArray, 0, 0, rmm_input_cpu.data, sizeof(float)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);
//  cudaBindTextureToArray(rmm_input_gpu_tex, cuArray);

//  rmm_input_gpu_tex.normalized = false;

//#if USE_LIBXC
//  const int nspin = XC_UNPOLARIZED;
//  const int functionalExchange = fortran_vars.ex_functional_id;
//  const int functionalCorrelation = fortran_vars.ec_functional_id;
//  LibxcProxy<float,3> libxcProxy(functionalExchange, functionalCorrelation, nspin);
//#endif


//  if (compute_energy) {
//    CudaMatrix<float> energy_gpu(number_of_points);

//#define compute_parameters \
//        energy_gpu.data, factors_gpu.data, point_weights_gpu.data, number_of_points, function_values_transposed.data, \
//        gradient_values_transposed.data, hessian_values_transposed.data, group_m, partial_densities_gpu.data, dxyz_gpu.data, \
//        dd1_gpu.data,dd2_gpu.data

//#define accumulate_parameters \
//        energy_gpu.data, factors_gpu.data, point_weights_gpu.data, number_of_points, block_height, \
//        partial_densities_gpu.data, dxyz_gpu.data, dd1_gpu.data, dd2_gpu.data

// VER QUE PASA SI SACAMOS COMPUTE_FACTOR Y COMPUTE ENERGY DE gpu_compute_density
//    if (compute_forces || compute_rmm) {
//      if (lda)
//      {
//          gpu_compute_density<float, true, true, true><<<threadGrid, threadBlock>>>(compute_parameters);
//          gpu_accumulate_point<float, true, true, true><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
//      }
//      else
//      {
//          gpu_compute_density<float, true, true, false><<<threadGrid, threadBlock>>>(compute_parameters);
//	  // TODO: aca tiene que ir el proxy a libxc.
//          gpu_accumulate_point<float, true, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
//      }
//    }
//    else {
//      if (lda)
//      {
//          gpu_compute_density<float, true, false, true><<<threadGrid, threadBlock>>>(compute_parameters);
//          gpu_accumulate_point<float, true, false, true><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
//      }
//      else
//      {
//          gpu_compute_density<float, true, false, false><<<threadGrid, threadBlock>>>(compute_parameters);
//	// TODO: aca tiene q ir el proxy a libxc.
//          gpu_accumulate_point<float, true, false, false><<<threadGrid_accumulate, threadBlock_accumulate>>> (accumulate_parameters);
//      }
//    }
//    cudaAssertNoError("compute_density");

//    HostMatrix<float> energy_cpu(energy_gpu);
//    for (uint i = 0; i < number_of_points; i++) {
//      energy += energy_cpu(i);
//    }
//  }
//  else {
//#undef compute_parameters
//#undef accumulate_parameters

//#define compute_parameters \
//    NULL,factors_gpu.data,point_weights_gpu.data,number_of_points,function_values_transposed.data,gradient_values_transposed.data,hessian_values_transposed.data,group_m,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data
//#define accumulate_parameters \
//    NULL,factors_gpu.data,point_weights_gpu.data,number_of_points,block_height,partial_densities_gpu.data,dxyz_gpu.data,dd1_gpu.data,dd2_gpu.data
//    if (lda)
//    {
//        gpu_compute_density<float, false, true, true><<<threadGrid, threadBlock>>>(compute_parameters);
//        gpu_accumulate_point<float, false, true, true><<<threadGrid_accumulate, threadBlock_accumulate>>>(accumulate_parameters);
//    }
//    else
//    {
//        gpu_compute_density<float, false, true, false><<<threadGrid, threadBlock>>>(compute_parameters);
//        gpu_accumulate_point<float, false, true, false><<<threadGrid_accumulate, threadBlock_accumulate>>>(accumulate_parameters);
//    }
//    cudaAssertNoError("compute_density");
//  }
//#undef compute_parameters
//#undef accumulate_parameters

//  timers.density.pause_and_sync();
  /* compute forces */
//  if (compute_forces) {
    //************ Repongo los valores que puse a cero antes, para las fuerzas son necesarios (o por lo mens utiles)
//    for (uint i=0; i<(group_m); i++) {
//      for(uint j=0; j<(group_m); j++) {
//        if((i>=group_m) || (j>=group_m) || (j > i))
//        {
//          rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=rmm_input_cpu.data[COALESCED_DIMENSION(group_m)*j+i] ;
//        }
//      }
//    }

//    timers.density_derivs.start_and_sync();
//    cudaMemcpyToArray(cuArray, 0, 0,rmm_input_cpu.data,
//      sizeof(float)*rmm_input_cpu.width*rmm_input_cpu.height, cudaMemcpyHostToDevice);

//    timers.density_derivs.start_and_sync();
//    dim3 threads = dim3(number_of_points);
//    threadBlock = dim3(DENSITY_DERIV_BLOCK_SIZE);
//    threadGrid = divUp(threads, threadBlock);

//    CudaMatrix<vec_type4> dd_gpu(COALESCED_DIMENSION(number_of_points), total_nucleii); dd_gpu.zero();
//    CudaMatrixUInt nuc_gpu(func2local_nuc);  // TODO: esto en realidad se podria guardar una sola vez durante su construccion

//    gpu_compute_density_derivs<<<threadGrid, threadBlock>>>(
//        function_values.data, gradient_values.data, nuc_gpu.data, dd_gpu.data, number_of_points, group_m, total_nucleii());
//    cudaAssertNoError("density_derivs");
//    timers.density_derivs.pause_and_sync();

//    timers.forces.start_and_sync();
//    CudaMatrix<vec_type4> forces_gpu(total_nucleii());

//    threads = dim3(total_nucleii());
//    threadBlock = dim3(FORCE_BLOCK_SIZE);
//    threadGrid = divUp(threads, threadBlock);
//    gpu_compute_forces<<<threadGrid, threadBlock>>>(
//        number_of_points, factors_gpu.data, dd_gpu.data, forces_gpu.data, total_nucleii());
//    cudaAssertNoError("forces");

//    HostMatrix<vec_type4> forces_cpu(forces_gpu);

//    for (uint i = 0; i < total_nucleii(); ++i) {
//      vec_type4 atom_force = forces_cpu(i);
//      uint global_nuc = local2global_nuc[i];
//      fort_forces_ms(global_nuc, 0) += atom_force.x;
//      fort_forces_ms(global_nuc, 1) += atom_force.y;
//      fort_forces_ms(global_nuc, 2) += atom_force.z;

//    }
//    timers.forces.pause_and_sync();
//  }

//  timers.rmm.start_and_sync();
  /* compute RMM */
//  if (compute_rmm) {
//    threadBlock = dim3(RMM_BLOCK_SIZE_XY, RMM_BLOCK_SIZE_XY);
//    uint blocksPerRow = divUp(group_m, RMM_BLOCK_SIZE_XY);
    // Only use enough blocks for lower triangle
//    threadGrid = dim3(blocksPerRow*(blocksPerRow+1)/2);

//    CudaMatrix<float> rmm_output_gpu(COALESCED_DIMENSION(group_m), group_m);
    // For calls with a single block (pretty common with cubes) don't bother doing the arithmetic to get block position in the matrix
//    if (blocksPerRow > 1) {
//        gpu_update_rmm<float,true><<<threadGrid, threadBlock>>>(factors_gpu.data, number_of_points, rmm_output_gpu.data, function_values.data, group_m);
//    } else {
//        gpu_update_rmm<float,false><<<threadGrid, threadBlock>>>(factors_gpu.data, number_of_points, rmm_output_gpu.data, function_values.data, group_m);
//    }
//    cudaAssertNoError("update_rmm");

    /*** Contribute this  RMM to the total RMM ***/
//    HostMatrix<float> rmm_output_cpu(rmm_output_gpu);
//    add_rmm_output(rmm_output_cpu, rmm_output_local);

//  }
//  timers.rmm.pause_and_sync();

  /* clear functions */
//  if(!(inGlobal)) {
//    function_values.deallocate();
//    gradient_values.deallocate();
//    hessian_values_transposed.deallocate();
//  }
  //Deshago el bind de textura de rmm
//  cudaUnbindTexture(rmm_input_gpu_tex); //Enroque el Unbind con el Free, asi parece mas logico. Nano
//  cudaFreeArray(cuArray);
}

//////////////////////////////////////
//// KERNELS
__global__ void kernel_matrix_test0002 (float* const theMatrix, int n, int m)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n+m) {
	// Estoy en rango.
	//printf("%f,",theMatrix[idx]);
    }
}

__global__ void kernel_matrix_test0003 (G2G::vec_type<float,4>* dxyz,
    G2G::vec_type<float,4>* dd1, G2G::vec_type<float,4>* dd2, int n, int m)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n+m) {
	// Estoy en rango.
	//printf("%f,",theMatrix[idx]);
    }
}

__global__ void kernel_matrix_test0004 (G2G::vec_type<double,4>* dxyz,
    G2G::vec_type<double,4>* dd1, G2G::vec_type<double,4>* dd2, int n, int m)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n+m) {
	// Estoy en rango.
	//printf("%f,",theMatrix[idx]);
    }
}


//////////////////////////////////////
//// TESTS

void matrix_test0001 () 
{
    printf("**  matrix_test0001  **\n");

    // Creamos una cuda matrix que ya
    // reserva memoria memoria en el device
    G2G::CudaMatrix<float> aHostMatrix();

}

void matrix_test0002 () 
{
    printf("**  matrix_test0002  **\n");
    int n = 5;
    int m = 5;
    G2G::CudaMatrix<float> partial_densities_gpu;
//  CudaMatrix< vec_type<float,4> > dxyz_gpu;
//  CudaMatrix< vec_type<float,4> > dd1_gpu;
//  CudaMatrix< vec_type<float,4> > dd2_gpu;

    //partial_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
    partial_densities_gpu.resize(n,m);
//  dxyz_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height);
//  dd1_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );
//  dd2_gpu.resize(COALESCED_DIMENSION(number_of_points),block_height );


    // Launch the Vector Add CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    //vectorAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, numElements);
    // Llamamos al kernel
    kernel_matrix_test0002<<<blocksPerGrid,threadsPerBlock>>> (partial_densities_gpu.data, n, m);

}


void matrix_test0003 ()
{
    printf("**  matrix_test0003  **\n");

    cudaError_t err = cudaSuccess;

    uint n = 5;
    uint m = 5;

    // Allocate the host input vector A
    uint size = (n+m) * sizeof(float); 
    float *h_A = (float *)malloc(size);

    G2G::CudaMatrix<float> partial_densities_gpu;
    partial_densities_gpu.resize(n,m);
    partial_densities_gpu.zero();

    // Initialize the host input vectors
    for (int i = 0; i < n+m; ++i)
    {
        h_A[i] = -1;
    }

    // Print the original data
    printf("Array original data:\n");
    for (int i = 0; i < n+m; ++i)
    {
        printf ("%f,", h_A[i]);
    }
    printf("\n");

    partial_densities_gpu.resize(n,m);

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

    // Llamamos al kernel
    kernel_matrix_test0002<<<blocksPerGrid,threadsPerBlock>>> (partial_densities_gpu.data, n, m);

    // Copy the device result vector in device memory to the host result vector
    // in host memory.
    printf("Copy output data from the CUDA device to the host memory\n");
    err = cudaMemcpy(h_A, partial_densities_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector A from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Print the results
    printf("Array with data from the matrix in CUDA\n");
    for (int i = 0; i < n+m; ++i)
    {
        printf ("%f,", h_A[i]);
    }
    printf("\n");

    // We modify the array data.
    for (int i = 0; i < n+m; ++i)
    {
        h_A[i] = 0.1 * i;
    }

    // La volvemos a enviar al kernel modificada.
    err = cudaMemcpy(partial_densities_gpu.data, h_A, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector A from host to device!\n");
        exit(EXIT_FAILURE);
    }

    // Reset the array.
    for (int i = 0; i < n+m; ++i)
    {
        h_A[i] = -1;
    }

    // Call the kernel again.
    kernel_matrix_test0002<<<blocksPerGrid,threadsPerBlock>>> (partial_densities_gpu.data, n, m);

    printf("Copy output data from the CUDA device to the host memory 2\n");
    err = cudaMemcpy(h_A, partial_densities_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector A from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Print the results
    printf("Array with data from the matrix in CUDA 2\n");
    for (int i = 0; i < n+m; ++i)
    {
        printf ("%f,", h_A[i]);
    }
    printf("\n");


    // Free memory.
    free(h_A);
}


void libxc_cpu_accumulate_point_local(LibxcProxy<float, 4>* libxcProxy, 
    float* const energy, float* const factor, const float* const point_weights,
    uint points, int block_height, float* partial_density, 
    G2G::vec_type<float,4>* dxyz, G2G::vec_type<float,4>* dd1, G2G::vec_type<float,4>* dd2) 
{
    //printf("libxc_cpu_accumulate_point_local()\n");
    float point_weight = 0.0f;
    float y2a, exc_corr, exc_c, exc_x;

    float _partial_density(0.0f);
    G2G::vec_type<float,4> _dxyz, _dd1, _dd2;

    _dxyz = _dd1 = _dd2 = G2G::vec_type<float,4>(0.0f,0.0f,0.0f,0.0f);

    for(int i=0; i<points; i++) {
	point_weight = point_weights[i];
	for(int j=0; j<block_height; j++) {
	    const int this_row = j*points+i;
	    _partial_density += partial_density[this_row];
	    _dxyz += dxyz[this_row];
	    _dd1 += dd1[this_row];
	    _dd2 += dd2[this_row];
	}
	
	// TODO: aca va la llamada al proxy.
	//calc_ggaCS_in<scalar_type, 4>(_partial_density, _dxyz, _dd1, _dd2, exc_x, exc_c, y2a, 9);
	if (libxcProxy != NULL) {
    	    libxcProxy->doGGA(_partial_density, _dxyz, _dd1, _dd2, exc_x, exc_c, y2a);
	}

	exc_corr = exc_x + exc_c;

	//if(compute_energy) {
	//    energy[i] = (_partial_density * point_weight) * exc_corr;
	//}
	energy[i] = (_partial_density * point_weight) * exc_corr;

	//if(compute_factor) {
	//    factor[i] = point_weight * y2a;
	//}
	factor[i] = point_weight * y2a;
    }

}


void matrix_test0004 ()
{
    printf("**  matrix_test0004  **\n");
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    G2G::CudaMatrix< G2G::vec_type<float,4> > dxyz_gpu;
    G2G::CudaMatrix< G2G::vec_type<float,4> > dd1_gpu;
    G2G::CudaMatrix< G2G::vec_type<float,4> > dd2_gpu;

    // Partial_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
    dxyz_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu.resize(COALESCED_DIMENSION(number_of_points),m);

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    // Call the kernel
    kernel_matrix_test0003<<<blocksPerGrid,threadsPerBlock>>> (dxyz_gpu.data, dd1_gpu.data, dd2_gpu.data, n, m);

    // Copy back to host the matrix data in gpu memory for the same datatype vec_type<float,4> and
    // call the new accumulate_point_cpu function to see what happens.

    G2G::vec_type<float,4>* dxyz_cpu;
    G2G::vec_type<float,4>* dd1_cpu;
    G2G::vec_type<float,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = (n+m) * sizeof(G2G::vec_type<float,4>); 
    dxyz_cpu = (G2G::vec_type<float,4> *)malloc(size);
    dd1_cpu = (G2G::vec_type<float,4> *)malloc(size);
    dd2_cpu = (G2G::vec_type<float,4> *)malloc(size);    

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Now the arrays for energy, factors, point_weight and partial_density
    float *energy_gpu = NULL;
    float *factor_gpu = NULL;
    float *point_weights_gpu = NULL;
    float *partial_density_gpu = NULL;

    // Create the arrays in CUDA memory.
    err = cudaMalloc((void**)&energy_gpu, size);
    if (err != cudaSuccess) 
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu, size);
    if (err != cudaSuccess) 
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu, size);
    if (err != cudaSuccess) 
    {
	printf("Failed to allocate vector point_weights_gpu!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu, size);
    if (err != cudaSuccess) 
    {
	printf("Failed to allocate vector partial_density_gpu!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu, -1, size);
    cudaMemset(factor_gpu, -2, size);
    cudaMemset(point_weights_gpu, -3, size);
    cudaMemset(partial_density_gpu, -4, size);

    // Allocate the host input vectors
    float *energy_cpu = (float *)malloc(size);
    float *factor_cpu = (float *)malloc(size);
    float *point_weights_cpu = (float *)malloc(size);
    float *partial_density_cpu = (float *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu, energy_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu, factor_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(point_weights_cpu, point_weights_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_cpu, partial_density_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Call the CUDA KERNEL
//    libxc_cpu_accumulate_point_local(NULL, energy_cpu, factor_cpu, point_weights_cpu, 
//	number_of_points, 1, partial_density_cpu,
//	dxyz_cpu, dd1_cpu, dd2_cpu);

    // TODO: now copy back the results to the gpu.
    
    err = cudaMemcpy(energy_gpu, energy_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_gpu, factor_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_cpu from host to device!\n");
        exit(EXIT_FAILURE);
    }
    
    err = cudaMemcpy(point_weights_gpu, point_weights_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_gpu, partial_density_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }
    
    // Copy the matrix data.
    err = cudaMemcpy(dxyz_gpu.data, dxyz_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dxyz_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dd1_gpu.data, dd1_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd1_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dd2_gpu.data, dd2_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd2_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }

    // Free jack :)
    free(energy_cpu);
    free(factor_cpu);
    free(point_weights_cpu);
    free(partial_density_cpu);

    free(dd1_cpu);
    free(dd2_cpu);
    free(dxyz_cpu);

    cudaFree(energy_gpu);
    cudaFree(factor_gpu);
    cudaFree(point_weights_gpu);
    cudaFree(partial_density_gpu);

}


void matrix_test0005 ()
{
    printf("**  matrix_test0005 - con libxcProxy  **\n");
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    G2G::CudaMatrix< G2G::vec_type<float,4> > dxyz_gpu;
    G2G::CudaMatrix< G2G::vec_type<float,4> > dd1_gpu;
    G2G::CudaMatrix< G2G::vec_type<float,4> > dd2_gpu;

    // Partial_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
    dxyz_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu.resize(COALESCED_DIMENSION(number_of_points),m);

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    // Call the kernel
    kernel_matrix_test0003<<<blocksPerGrid,threadsPerBlock>>> (dxyz_gpu.data, dd1_gpu.data, dd2_gpu.data, n, m);

    // Copy back to host the matrix data in gpu memory for the same datatype vec_type<float,4> and
    // call the new accumulate_point_cpu function to see what happens.

    G2G::vec_type<float,4>* dxyz_cpu;
    G2G::vec_type<float,4>* dd1_cpu;
    G2G::vec_type<float,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = (n+m) * sizeof(G2G::vec_type<float,4>); 
    dxyz_cpu = (G2G::vec_type<float,4> *)malloc(size);
    dd1_cpu = (G2G::vec_type<float,4> *)malloc(size);
    dd2_cpu = (G2G::vec_type<float,4> *)malloc(size);    

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Now the arrays for energy, factors, point_weight and partial_density
    float *energy_gpu = NULL;
    float *factor_gpu = NULL;
    float *point_weights_gpu = NULL;
    float *partial_density_gpu = NULL;

    // Create the arrays in CUDA memory.
    err = cudaMalloc((void**)&energy_gpu, size);
    if (err != cudaSuccess) 
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu, size);
    if (err != cudaSuccess) 
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu, size);
    if (err != cudaSuccess) 
    {
	printf("Failed to allocate vector point_weights_gpu!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu, size);
    if (err != cudaSuccess) 
    {
	printf("Failed to allocate vector partial_density_gpu!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu, -1, size);
    cudaMemset(factor_gpu, -2, size);
    cudaMemset(point_weights_gpu, -3, size);
    cudaMemset(partial_density_gpu, -4, size);

    // Allocate the host input vectors
    float *energy_cpu = (float *)malloc(size);
    float *factor_cpu = (float *)malloc(size);
    float *point_weights_cpu = (float *)malloc(size);
    float *partial_density_cpu = (float *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu, energy_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu, factor_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(point_weights_cpu, point_weights_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_cpu, partial_density_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    const int nspin = 1;
    const int functionalExchange = 101;
    const int functionalCorrelation = 130;
    LibxcProxy<float,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    // Call the CUDA KERNEL
//    libxc_cpu_accumulate_point_local(&libxcProxy, energy_cpu, factor_cpu, point_weights_cpu, 
//	number_of_points, 1, partial_density_cpu,
//	dxyz_cpu, dd1_cpu, dd2_cpu);

    // TODO: now copy back the results to the gpu.
    
    err = cudaMemcpy(energy_gpu, energy_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_gpu, factor_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_cpu from host to device!\n");
        exit(EXIT_FAILURE);
    }
    
    err = cudaMemcpy(point_weights_gpu, point_weights_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_gpu, partial_density_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }
    
    // Copy the matrix data.
    err = cudaMemcpy(dxyz_gpu.data, dxyz_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dxyz_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dd1_gpu.data, dd1_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd1_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(dd2_gpu.data, dd2_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd2_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }

    // Free jack :)
    free(energy_cpu);
    free(factor_cpu);
    free(point_weights_cpu);
    free(partial_density_cpu);

    free(dd1_cpu);
    free(dd2_cpu);
    free(dxyz_cpu);

    cudaFree(energy_gpu);
    cudaFree(factor_gpu);
    cudaFree(point_weights_gpu);
    cudaFree(partial_density_gpu);

}


void matrix_test0006 ()
{
    printf("**  matrix_test0006 - con libxcProxy 2  **\n");
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    G2G::CudaMatrix< G2G::vec_type<float,4> > dxyz_gpu;
    G2G::CudaMatrix< G2G::vec_type<float,4> > dd1_gpu;
    G2G::CudaMatrix< G2G::vec_type<float,4> > dd2_gpu;

    // Partial_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
    dxyz_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu.resize(COALESCED_DIMENSION(number_of_points),m);

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    // Call the kernel
    kernel_matrix_test0003<<<blocksPerGrid,threadsPerBlock>>> (dxyz_gpu.data, dd1_gpu.data, dd2_gpu.data, n, m);

    // Copy back to host the matrix data in gpu memory for the same datatype vec_type<float,4> and
    // call the new accumulate_point_cpu function to see what happens.

    G2G::vec_type<float,4>* dxyz_cpu;
    G2G::vec_type<float,4>* dd1_cpu;
    G2G::vec_type<float,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = (n+m) * sizeof(G2G::vec_type<float,4>);
    dxyz_cpu = (G2G::vec_type<float,4> *)malloc(size);
    dd1_cpu = (G2G::vec_type<float,4> *)malloc(size);
    dd2_cpu = (G2G::vec_type<float,4> *)malloc(size);

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Now the arrays for energy, factors, point_weight and partial_density
    float *energy_gpu = NULL;
    float *factor_gpu = NULL;
    float *point_weights_gpu = NULL;
    float *partial_density_gpu = NULL;

    // Create the arrays in CUDA memory.
    err = cudaMalloc((void**)&energy_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu, -1, size);
    cudaMemset(factor_gpu, -2, size);
    cudaMemset(point_weights_gpu, -3, size);
    cudaMemset(partial_density_gpu, -4, size);

    // Allocate the host input vectors
    float *energy_cpu = (float *)malloc(size);
    float *factor_cpu = (float *)malloc(size);
    float *point_weights_cpu = (float *)malloc(size);
    float *partial_density_cpu = (float *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu, energy_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu, factor_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(point_weights_cpu, point_weights_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_cpu, partial_density_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    const int nspin = 1;
    const int functionalExchange = 101;
    const int functionalCorrelation = 130;
    LibxcProxy<float,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    // Call the CUDA KERNEL
//    libxc_cpu_accumulate_point(&libxcProxy, energy_cpu, factor_cpu, point_weights_cpu, 
//	number_of_points, 1, partial_density_cpu,
//	dxyz_cpu, dd1_cpu, dd2_cpu);

//    <float, true, true, false> 
//    libxc_cpu_accumulate_point<float, true, true, false>(&libxcProxy, energy_cpu, 
//	factor_cpu, point_weights_cpu,
//        number_of_points, 1, partial_density_cpu, 
//	dxyz_cpu, dd1_cpu, dd2_cpu);

    // TODO: now copy back the results to the gpu.

    err = cudaMemcpy(energy_gpu, energy_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_gpu, factor_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_cpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(point_weights_gpu, point_weights_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_gpu, partial_density_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    // Copy the matrix data.
    err = cudaMemcpy(dxyz_gpu.data, dxyz_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dxyz_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_gpu.data, dd1_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd1_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_gpu.data, dd2_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd2_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }

    // Free jack :)
    free(energy_cpu);
    free(factor_cpu);
    free(point_weights_cpu);
    free(partial_density_cpu);

    free(dd1_cpu);
    free(dd2_cpu);
    free(dxyz_cpu);

    cudaFree(energy_gpu);
    cudaFree(factor_gpu);
    cudaFree(point_weights_gpu);
    cudaFree(partial_density_gpu);

}


void matrix_test0007 ()
{
    printf("**  matrix_test0007 - con libxcProxy 3  **\n");
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu;

    // Partial_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
    dxyz_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu.resize(COALESCED_DIMENSION(number_of_points),m);

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    // Call the kernel
    kernel_matrix_test0004<<<blocksPerGrid,threadsPerBlock>>> (dxyz_gpu.data, dd1_gpu.data, dd2_gpu.data, n, m);

    // Copy back to host the matrix data in gpu memory for the same datatype vec_type<double,4> and
    // call the new accumulate_point_cpu function to see what happens.

    G2G::vec_type<double,4>* dxyz_cpu;
    G2G::vec_type<double,4>* dd1_cpu;
    G2G::vec_type<double,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = (n+m) * sizeof(G2G::vec_type<double,4>);
    dxyz_cpu = (G2G::vec_type<double,4> *)malloc(size);
    dd1_cpu = (G2G::vec_type<double,4> *)malloc(size);
    dd2_cpu = (G2G::vec_type<double,4> *)malloc(size);

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Now the arrays for energy, factors, point_weight and partial_density
    double *energy_gpu = NULL;
    double *factor_gpu = NULL;
    double *point_weights_gpu = NULL;
    double *partial_density_gpu = NULL;

    // Create the arrays in CUDA memory.
    err = cudaMalloc((void**)&energy_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu, -1, size);
    cudaMemset(factor_gpu, -2, size);
    cudaMemset(point_weights_gpu, -3, size);
    cudaMemset(partial_density_gpu, -4, size);

    // Allocate the host input vectors
    double *energy_cpu = (double *)malloc(size);
    double *factor_cpu = (double *)malloc(size);
    double *point_weights_cpu = (double *)malloc(size);
    double *partial_density_cpu = (double *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu, energy_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu, factor_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(point_weights_cpu, point_weights_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_cpu, partial_density_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    const int nspin = 1;
    const int functionalExchange = 101;
    const int functionalCorrelation = 130;
    LibxcProxy<double,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    // Call the CUDA KERNEL
//    libxc_cpu_accumulate_point(&libxcProxy, energy_cpu, factor_cpu, point_weights_cpu, 
//	number_of_points, 1, partial_density_cpu,
//	dxyz_cpu, dd1_cpu, dd2_cpu);

//    <double, true, true, false> 
//    libxc_cpu_accumulate_point<double, true, true, false>(&libxcProxy, energy_cpu, 
//	factor_cpu, point_weights_cpu,
//        number_of_points, 1, partial_density_cpu, 
//	dxyz_cpu, dd1_cpu, dd2_cpu);

    // TODO: now copy back the results to the gpu.

    err = cudaMemcpy(energy_gpu, energy_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_gpu, factor_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_cpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(point_weights_gpu, point_weights_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_gpu, partial_density_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    // Copy the matrix data.
    err = cudaMemcpy(dxyz_gpu.data, dxyz_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dxyz_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_gpu.data, dd1_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd1_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_gpu.data, dd2_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd2_cpu A from host to device!\n");
        exit(EXIT_FAILURE);
    }

    // Free jack :)
    free(energy_cpu);
    free(factor_cpu);
    free(point_weights_cpu);
    free(partial_density_cpu);

    free(dd1_cpu);
    free(dd2_cpu);
    free(dxyz_cpu);

    cudaFree(energy_gpu);
    cudaFree(factor_gpu);
    cudaFree(point_weights_gpu);
    cudaFree(partial_density_gpu);

}

void matrix_test0008()
{
  printf("matrix_test0008()\n");
  G2G::CudaMatrix<double> point_weights_gpu;
  G2G::HostMatrix<double> point_weights_cpu(5, 1);

  for (int i=0; i<5; i++){
    point_weights_cpu(i) = 0.001*i;
  }

  point_weights_gpu = point_weights_cpu;
}

void matrix_test0009()
{
  printf("matrix_test0009()\n");
  //typedef vec_type<scalar_type,2> vec_type2;
  //typedef vec_type<scalar_type,3> vec_type3;
  typedef G2G::vec_type<float,4> vec_type4;
  //G2G::CudaMatrix<scalar_type> function_values;
  //G2G::CudaMatrix<vec_type4> gradient_values;
  //G2G::CudaMatrix<vec_type4> hessian_values_transposed;

  G2G::CudaMatrix<vec_type4> point_weights_gpu;
  G2G::HostMatrix<vec_type4> point_weights_cpu(5, 1);

  G2G::vec_type<float,4> one(1,1,1,1);

  for (int i=0; i<5; i++){
    point_weights_cpu(i).x = one.x;
  }

  point_weights_gpu = point_weights_cpu;
}

void matrix_test0010()
{
    printf("matrix_test0010()\n");

    G2G::CudaMatrix< G2G::vec_type<float,4> > point_weights_gpu;
    G2G::HostMatrix< G2G::vec_type<float,4> > point_weights_cpu(5, 1);

    G2G::vec_type<float,4> one(1,1,1,1);

    for (int i=0; i<5; i++){
	point_weights_cpu(i).x = one.x;
    }

    point_weights_gpu = point_weights_cpu;

    //int vec_size = point_weights_cpu.elements() * sizeof (G2G::vec_type<float,4>);
    int vec_size = 10;
    //unsigned int vec_size = point_weights_cpu.bytes();
    G2G::vec_type<float,4>* vectors = (G2G::vec_type<float,4>*)malloc(vec_size);

    for (int i=0; i<5; i++) {
	vectors[i].x = 2;
    }
}

/////////////////////////////////////
//// MAIN
int main(int argc, char **argv)
{
    printf("*************************\n");
    printf("**  Matrix Unit Tests  **\n");
    printf("*************************\n");

    //solve_closed(false, false, false, true, 0, NULL, 1,NULL);
    matrix_test0001();
    matrix_test0002();
    matrix_test0003();
    matrix_test0004();
    matrix_test0005();
    matrix_test0006();
    matrix_test0007();
    matrix_test0008();

    printf("*************************\n");
    printf("**      Test End       **\n");
    printf("*************************\n");

    return 0;
}