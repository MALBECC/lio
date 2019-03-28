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

////////////////////////
// Dummy solve_closed
//
void solve_closed( bool compute_rmm, bool lda, bool compute_forces, bool compute_energy,
    double& energy,    G2G::HostMatrix<double>& fort_forces_ms,
    int inner_threads, G2G::HostMatrix<double>& rmm_output_local) 
{
    // This is used to simulate a call inside the tests
}

//////////////////////////////////////
//// KERNEL UTILS
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

/////////////////////////////////////////////
// Test: matrix_test0001
//
// We test the CudaMatrix constructor
//
void matrix_test0001 () 
{
    printf("**  matrix_test0001  **\n");

    // Creamos una cuda matrix que ya
    // reserva memoria memoria en el device
    G2G::CudaMatrix<float> aHostMatrix();

}

/////////////////////////////////////////////
// Test: matrix_test0002
//
// We test the CudaMatrix 'resize' function
//
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

/////////////////////////////////////////////
// Test: matrix_test0003
//
// We test the copy of a single CudaMatrix
// into a kernel.
//
//

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

/////////////////////////////////////////////
// Test: libxc_cpu_accumulate_point_local
//
// Test to simulate the accumulate point
// function before we call the CPU version
// of the LibxcProxy component.
//

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

/////////////////////////////////////////////
// Test: matrix_test0004
//
// We test the use of the CudaMatrix class with
// with the LibxcProxy component. We only test
// the alloc and free of the parameters
// needed by the LibxcProxy component.
//
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

/////////////////////////////////////////////
// Test: matrix_test0005
//
// We test the use of the CudaMatrix class with
// with the LibxcProxy component.
// We set all the parameters needed by the.
// LibxcProxy component.
//
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
    const int functionalExchange = 1101;
    const int functionalCorrelation = 1130;
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

/////////////////////////////////////////////
// Test: matrix_test0006
//
// We test the use of the CudaMatrix class with
// with the LibxcProxy component.
//
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
    const int functionalExchange = 1101;
    const int functionalCorrelation = 1130;
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

/////////////////////////////////////////////
// Test: matrix_test0007
//
// We test the use of the CudaMatrix class with
// with the LibxcProxy component.
//
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
    const int functionalExchange = 1101;
    const int functionalCorrelation = 1130;
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

/////////////////////////////////////////////
// Test: matrix_test00008
//
// Simple matrix test to check how the matrix
// is copied to cuda device.
//
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

/////////////////////////////////////////////
// Test: matrix_test00009
//
// Simple matrix test to check how the matrix
// is copied to cuda device.
//
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

/////////////////////////////////////////////
// Test: matrix_test00010
//
// Simple matrix test to check how the matrix
// is copied to cuda device.
//
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

    try {
	matrix_test0001();
	matrix_test0002();
	matrix_test0003();
	matrix_test0004();
	matrix_test0005();
	matrix_test0006();
        matrix_test0007();
	matrix_test0008();
    } catch (int e) {
	printf("An exception ocurred: %u \n", e);
	exit (EXIT_FAILURE);
    }
    printf("*************************\n");
    printf("**      Test End       **\n");
    printf("*************************\n");

    return 0;
}