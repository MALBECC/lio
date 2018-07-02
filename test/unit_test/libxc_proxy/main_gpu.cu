#include <iostream>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math_constants.h>
#include <float.h>
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
#include "../../../g2g/libxc/libxc_accumulate_point.h"
#include "../commons/test_input.h"
#include <typeinfo>

using namespace std;

using std::cout;
using std::endl;

__constant__ double TWO = 2;

void accumulate_data_for_libxc_test0001()
{
#if FULL_DOUBLE
    printf("** accumulate_data_for_libxc_test0001 **\n");

    cudaError_t err = cudaSuccess;
    uint n = 1000;
    uint m = 1000;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in;

    // Accumulate
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_accum;

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),1);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),1);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),1);

    dxyz_gpu_in.zero();
    dd1_gpu_in.zero();
    dd2_gpu_in.zero();

    dxyz_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),1);
    dd1_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),1);
    dd2_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),1);

    dxyz_gpu_accum.zero();
    dd1_gpu_accum.zero();
    dd2_gpu_accum.zero();

    // Now the arrays for energy, factors, point_weight and partial_density
    double *point_weights_gpu_in = NULL;
    double *partial_density_gpu_in = NULL;
    // Accum
    double *accumulated_density_gpu = NULL;

    // Now the arrays for energy, factors, point_weight and partial_density
    double *energy_gpu_in = NULL;
    double *factor_gpu_in = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);

    err = cudaMalloc((void**)&energy_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_in!\n");
    }

    err = cudaMalloc((void**)&accumulated_density_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_accum!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu_in, 0, size);
    cudaMemset(factor_gpu_in, 0, size);
    cudaMemset(point_weights_gpu_in, 0, size);
    cudaMemset(partial_density_gpu_in, 0, size);
    cudaMemset(accumulated_density_gpu, 0, size);

    // Launch the CUDA Kernel
    //int numElements = number_of_points;
    //int threadsPerBlock = 32;
    //int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    //uint block_height = 1;

    // Create the libxcproxy
    const int nspin = 1;
    const int functionalExchange = 1101;
    const int functionalCorrelation = 1130;
    LibxcProxy<double,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    // Calculate exc_corr and y2a
    libxc_exchange_correlation_gpu<double, true, true, false> (&libxcProxy,
	energy_gpu_in,
	factor_gpu_in,
	number_of_points,
	accumulated_density_gpu,
	dxyz_gpu_accum.data,
        dd1_gpu_accum.data,
	dd2_gpu_accum.data);

    // Check and print the results.

    // Free memory
    cudaFree(energy_gpu_in);
    cudaFree(factor_gpu_in);
    cudaFree(point_weights_gpu_in);
    cudaFree(partial_density_gpu_in);
    cudaFree(accumulated_density_gpu);
#endif
}

template <class T, int width>
__global__ void funcionDeMierda(
		    T* ex, double* exchange,
		    T* ec, double* correlation,
		    T* vrho, double* vrhoC,
		    T* vsigma, double* vsigmaC,
		    T* v2rho, double* v2rhoC,
		    T* v2rhosigma, double* v2rhosigmaC,
		    T* v2sigma, double* v2sigmaC,
		    T* y2a,
		    T* sigma,
		    G2G::vec_type<T, width>* grad,
		    G2G::vec_type<T, width>* hess1,
		    G2G::vec_type<T, width>* hess2,
		    int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
        printf("%i %lf %lf %lf\n",i, grad[i].x, grad[i].y, grad[i].z);

	ex[i] = exchange[i];
	ec[i] = correlation[i];

	// Merge the results for the derivatives.
	vrho[i] += vrhoC[i];
        vsigma[i] += vsigmaC[i];
        v2rho[i] += v2rhoC[i];
        v2rhosigma[i] += v2rhosigmaC[i];
        v2sigma[i] += v2sigmaC[i];
        // Now, compute y2a value.

	y2a[i] = vrho[i] - (2 * sigma[i] * v2rhosigma[i]
            + 2 * (hess1[i].x + hess1[i].y + hess1[i].z) * vsigma[i]
            + 4 * v2sigma[i] * (grad[i].x * grad[i].x * hess1[i].x + 
				grad[i].y * grad[i].y * hess1[i].y + 
				grad[i].z * grad[i].z * hess1[i].z + 
				2 * grad[i].x * grad[i].y * hess2[i].x + 
				2 * grad[i].x * grad[i].z * hess2[i].y + 
				2 * grad[i].y * grad[i].z * hess2[i].z));

    }
}

void joinResultsTest0001() {
    printf("joinResultsTest0001()\n");
    // Gather the results.
    // Variables for the Kernels
    int number_of_points = 10;
    int threadsPerBlock = 256;
    int blocksPerGrid = (number_of_points + threadsPerBlock - 1) / threadsPerBlock;

    cudaError_t err = cudaSuccess;
    int array_size = sizeof(double) * number_of_points;
    
    double* rho = NULL;
    err = cudaMalloc((void **)&rho, array_size);
    if (err != cudaSuccess) {
	fprintf(stderr, "Failed to allocate device rho! \n");
    }

    double* sigma = NULL;
    err = cudaMalloc((void**)&sigma, array_size);
    if (err != cudaSuccess) {
	fprintf(stderr, "Failed to allocate device sigma! \n");
    }

    double* exchange = NULL;
    err = cudaMalloc((void **)&exchange, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device exchange!\n");
        exit(EXIT_FAILURE);
    }

    double* correlation = NULL;
    err = cudaMalloc((void **)&correlation, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device correlation!\n");
        exit(EXIT_FAILURE);
    }
    // Clean arrays
    cudaMemset(exchange, 0, array_size);
    cudaMemset(correlation, 0, array_size);

    // The outputs for exchange
    double* vrho = NULL;
    double* vsigma = NULL;
    double* v2rho = NULL;
    double* v2rhosigma = NULL;
    double* v2sigma = NULL;

    err = cudaMalloc((void **)&vrho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vrho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&vsigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vsigma!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhosigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhosigma!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2sigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2sigma!\n");
        exit(EXIT_FAILURE);
    }
    // Clear arrays
    cudaMemset(vrho, 0, array_size);
    cudaMemset(vsigma, 0, array_size);
    cudaMemset(v2rho, 0, array_size);
    cudaMemset(v2rhosigma, 0, array_size);
    cudaMemset(v2sigma, 0, array_size);

    // The outputs for correlation
    double* vrhoC = NULL;
    double* vsigmaC = NULL;
    double* v2rhoC = NULL;
    double* v2rhosigmaC = NULL;
    double* v2sigmaC = NULL;

    err = cudaMalloc((void **)&vrhoC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vrhoC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&vsigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vsigmaC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhoC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhoC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhosigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhosigmaC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2sigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2sigmaC!\n");
        exit(EXIT_FAILURE);
    }

    double* ex = NULL;
    err = cudaMalloc((void **)&ex, array_size);
    if (err != cudaSuccess) {
	fprintf(stderr, "Failed to allocate device ex! \n");
    }

    double* ec = NULL;
    err = cudaMalloc((void**)&ec, array_size);
    if (err != cudaSuccess) {
	fprintf(stderr, "Failed to allocate device ec! \n");
    }

    double* y2a = NULL;
    err = cudaMalloc((void **)&y2a, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device y2a!\n");
        exit(EXIT_FAILURE);
    }

    // More data
    G2G::CudaMatrix< G2G::vec_type<double,4> > grad;
    G2G::CudaMatrix< G2G::vec_type<double,4> > hess1;
    G2G::CudaMatrix< G2G::vec_type<double,4> > hess2;

    grad.resize(COALESCED_DIMENSION(number_of_points),1);
    hess1.resize(COALESCED_DIMENSION(number_of_points),1);
    hess2.resize(COALESCED_DIMENSION(number_of_points),1);

    grad.zero();
    hess1.zero();
    hess2.zero();

    // Clear arrays
    cudaMemset(vrhoC, 0, array_size);
    cudaMemset(vsigmaC, 0, array_size);
    cudaMemset(v2rhoC, 0, array_size);
    cudaMemset(v2rhosigmaC, 0, array_size);
    cudaMemset(v2sigmaC, 0, array_size);
    cudaMemset(ex, 0, array_size);
    cudaMemset(ec, 0, array_size);
    cudaMemset(y2a, 0, array_size);

    // Gather the results.
    funcionDeMierda<double, 4><<<blocksPerGrid, threadsPerBlock>>>(
	ex, exchange,
	ec, correlation,
	vrho, vrhoC,
	vsigma, vsigmaC,
	v2rho, v2rhoC,
	v2rhosigma, v2rhosigmaC,
	v2sigma, v2sigmaC,
	y2a,
	sigma,
	grad.data,
	hess1.data,
	hess2.data,
	number_of_points);


    // Print the fucking results
    double* ex_cpu = (double*)malloc(array_size);
    double* ec_cpu = (double*)malloc(array_size);
    double* y2a_cpu = (double*)malloc(array_size);

    cudaMemcpy(ex_cpu, ex, array_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(ec_cpu, ec, array_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(y2a_cpu, y2a, array_size, cudaMemcpyDeviceToHost);

    printf("Results \n");
    for (int j=0; j<number_of_points; j++) {
	printf("%i %lf %lf %lf \n",j, ex_cpu[j], ec_cpu[j], y2a_cpu[j]);
    }

    // Free device memory.
    if (rho != NULL) {
	cudaFree(rho);
    }
    if (sigma != NULL) {
	cudaFree(sigma);
    }
    if (exchange != NULL) {
	cudaFree(exchange);
    }
    if (correlation != NULL) {
	cudaFree(correlation);
    }
    if (vrho != NULL) {
        cudaFree(vrho);
    }
    if (vsigma != NULL) {
	cudaFree(vsigma);
    }
    if (v2rho != NULL) {
	cudaFree(v2rho);
    }
    if (v2rhosigma != NULL) {
	cudaFree(v2rhosigma);
    }
    if (v2sigma != NULL) {
	cudaFree(v2sigma);
    }
    if (vrhoC != NULL) {
        cudaFree(vrhoC);
    }
    if (vsigmaC != NULL) {
	cudaFree(vsigmaC);
    }
    if (v2rhoC != NULL) {
	cudaFree(v2rhoC);
    }
    if (v2rhosigmaC != NULL) {
	cudaFree(v2rhosigmaC);
    }
    if (v2sigmaC != NULL) {
	cudaFree(v2sigmaC);
    }

    free(ex_cpu);
    free(ec_cpu);
    free(y2a_cpu);
}

#if FULL_DOUBLE
void doGGA_gpu (int number_of_points,
    double *dens_cpu,
    double *contracted_gradient_cpu,
    G2G::vec_type<double,4>* grad,
    G2G::vec_type<double,4>* hess1,
    G2G::vec_type<double,4>* hess2) {

    printf("doGAA_gpu(%i)\n", number_of_points);

    /////////////////////////////
    // CUDA ARRAYS
    double* accumulated_density_gpu = NULL;
    cudaMalloc ((void**)&accumulated_density_gpu, sizeof(double)*number_of_points);

    double* contracted_gradient_gpu = NULL;
    cudaMalloc ((void**)&contracted_gradient_gpu, sizeof(double)*number_of_points);

    G2G::vec_type<double,WIDTH>* dxyz_gpu = NULL;
    cudaMalloc((void**)&dxyz_gpu, sizeof(G2G::vec_type<double,4>)*number_of_points);

    G2G::vec_type<double,WIDTH>* dd1_gpu = NULL;
    cudaMalloc((void**)&dd1_gpu, sizeof(G2G::vec_type<double,4>)*number_of_points);

    G2G::vec_type<double,WIDTH>* dd2_gpu = NULL;
    cudaMalloc((void**)&dd2_gpu, sizeof(G2G::vec_type<double,4>)*number_of_points);

    double* ex_gpu = NULL;
    cudaMalloc ((void**)&ex_gpu, sizeof(double)*number_of_points);

    double* ec_gpu = NULL;
    cudaMalloc ((void**)&ec_gpu, sizeof(double)*number_of_points);

    double* y2a_gpu = NULL;
    cudaMalloc ((void**)&y2a_gpu, sizeof(double)*number_of_points);

    //////////////////////////////
    // SET CUDA ARRAYS VALUES
    cudaMemset(accumulated_density_gpu, 0, sizeof(double)*number_of_points);
    cudaMemset(contracted_gradient_gpu, 0, sizeof(double)*number_of_points);
    cudaMemset(ex_gpu, 0, sizeof(double)*number_of_points);
    cudaMemset(ex_gpu, 0, sizeof(double)*number_of_points);
    cudaMemset(y2a_gpu, 0, sizeof(double)*number_of_points);

    cudaMemcpy(accumulated_density_gpu, dens_cpu, sizeof(double)*number_of_points, cudaMemcpyHostToDevice);
    cudaMemcpy(contracted_gradient_gpu, contracted_gradient_cpu, sizeof(double)*number_of_points, cudaMemcpyHostToDevice);

    cudaMemcpy(dxyz_gpu, grad, sizeof(G2G::vec_type<double,4>)*number_of_points, cudaMemcpyHostToDevice);
    cudaMemcpy(dd1_gpu, hess1, sizeof(G2G::vec_type<double,4>)*number_of_points, cudaMemcpyHostToDevice);
    cudaMemcpy(dd2_gpu, hess2, sizeof(G2G::vec_type<double,4>)*number_of_points, cudaMemcpyHostToDevice);

    //////////////////////////////
    // CREATE THE PROXY
    const int nspin = 1;
    const int functionalExchange = 1101;
    const int functionalCorrelation = 1130;
    LibxcProxy<double,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    //////////////////////////////
    // MAKE THE CALLS
    libxcProxy.doGGA(
	accumulated_density_gpu,
	number_of_points,
	contracted_gradient_gpu,
	dxyz_gpu,
	dd1_gpu,
	dd2_gpu,
	ex_gpu,
	ec_gpu,
	y2a_gpu);

    /////////////////////////////
    // PRINT THE RESULTS
    double* energy_cpu = (double*)malloc(sizeof(double)*number_of_points);
    double* y2a_cpu = (double*)malloc(sizeof(double)*number_of_points);
    double* ex_cpu = (double*)malloc(sizeof(double)*number_of_points);
    double* ec_cpu = (double*)malloc(sizeof(double)*number_of_points);

    memset(ec_cpu,0,sizeof(double)*number_of_points);
    memset(ex_cpu,0,sizeof(double)*number_of_points);
    memset(y2a_cpu,0,sizeof(double)*number_of_points);

    cudaMemcpy(ec_cpu, ec_gpu, sizeof(double)*number_of_points, cudaMemcpyDeviceToHost);
    cudaMemcpy(ex_cpu, ex_gpu, sizeof(double)*number_of_points, cudaMemcpyDeviceToHost);
    cudaMemcpy(y2a_cpu, y2a_gpu, sizeof(double)*number_of_points, cudaMemcpyDeviceToHost);

    print_accumulate_point_data (NULL, NULL, NULL, ex_cpu, y2a_cpu, NULL, NULL, number_of_points);

    /////////////////////////////
    // FREE MEMORY
    cudaFree(ec_gpu);
    cudaFree(ex_gpu);
    cudaFree(y2a_gpu);
    cudaFree(accumulated_density_gpu);
    cudaFree(contracted_gradient_gpu);
    cudaFree(dxyz_gpu);
    cudaFree(dd1_gpu);
    cudaFree(dd2_gpu);

    free(ec_cpu);
    free(ex_cpu);
    free(y2a_cpu);
}
#endif

void doGGA_gpu_float (const int number_of_points,
    float *dens_cpu,
    float *contracted_gradient_cpu,
    G2G::vec_type<float,4>* grad,
    G2G::vec_type<float,4>* hess1,
    G2G::vec_type<float,4>* hess2) {

    printf("doGAA_gpu(%i)\n", number_of_points);

    /////////////////////////////
    // CUDA ARRAYS
    float* accumulated_density_gpu = NULL;
    cudaMalloc ((void**)&accumulated_density_gpu, sizeof(float)*number_of_points);

    float* contracted_gradient_gpu = NULL;
    cudaMalloc ((void**)&contracted_gradient_gpu, sizeof(float)*number_of_points);

    G2G::vec_type<float,WIDTH>* dxyz_gpu = NULL;
    cudaMalloc((void**)&dxyz_gpu, sizeof(G2G::vec_type<float,4>)*number_of_points);

    G2G::vec_type<float,WIDTH>* dd1_gpu = NULL;
    cudaMalloc((void**)&dd1_gpu, sizeof(G2G::vec_type<float,4>)*number_of_points);

    G2G::vec_type<float,WIDTH>* dd2_gpu = NULL;
    cudaMalloc((void**)&dd2_gpu, sizeof(G2G::vec_type<float,4>)*number_of_points);

    float* ex_gpu = NULL;
    cudaMalloc ((void**)&ex_gpu, sizeof(float)*number_of_points);

    float* ec_gpu = NULL;
    cudaMalloc ((void**)&ec_gpu, sizeof(float)*number_of_points);

    float* y2a_gpu = NULL;
    cudaMalloc ((void**)&y2a_gpu, sizeof(float)*number_of_points);

    //////////////////////////////
    // SET CUDA ARRAYS VALUES
    cudaMemset(accumulated_density_gpu, 0, sizeof(float)*number_of_points);
    cudaMemset(contracted_gradient_gpu, 0, sizeof(float)*number_of_points);
    cudaMemset(ex_gpu, 0, sizeof(float)*number_of_points);
    cudaMemset(ex_gpu, 0, sizeof(float)*number_of_points);
    cudaMemset(y2a_gpu, 0, sizeof(float)*number_of_points);

    cudaMemcpy(accumulated_density_gpu, dens_cpu, sizeof(float)*number_of_points, cudaMemcpyHostToDevice);
    cudaMemcpy(contracted_gradient_gpu, contracted_gradient_cpu, sizeof(float)*number_of_points, cudaMemcpyHostToDevice);

    cudaMemcpy(dxyz_gpu, grad, sizeof(G2G::vec_type<float,4>)*number_of_points, cudaMemcpyHostToDevice);
    cudaMemcpy(dd1_gpu, hess1, sizeof(G2G::vec_type<float,4>)*number_of_points, cudaMemcpyHostToDevice);
    cudaMemcpy(dd2_gpu, hess2, sizeof(G2G::vec_type<float,4>)*number_of_points, cudaMemcpyHostToDevice);

    //////////////////////////////
    // CREATE THE PROXY
    const int nspin = 1;
    const int functionalExchange = 1101;
    const int functionalCorrelation = 1130;
    LibxcProxy<float,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    //////////////////////////////
    // MAKE THE CALLS
    libxcProxy.doGGA(
	accumulated_density_gpu,
	number_of_points,
	contracted_gradient_gpu,
	dxyz_gpu,
	dd1_gpu,
	dd2_gpu,
	ex_gpu,
	ec_gpu,
	y2a_gpu);

    /////////////////////////////
    // PRINT THE RESULTS
    float* energy_cpu = (float*)malloc(sizeof(float)*number_of_points);
    float* y2a_cpu = (float*)malloc(sizeof(float)*number_of_points);
    float* ex_cpu = (float*)malloc(sizeof(float)*number_of_points);
    float* ec_cpu = (float*)malloc(sizeof(float)*number_of_points);

    memset(ec_cpu,0,sizeof(float)*number_of_points);
    memset(ex_cpu,0,sizeof(float)*number_of_points);
    memset(y2a_cpu,0,sizeof(float)*number_of_points);

    cudaMemcpy(ec_cpu, ec_gpu, sizeof(float)*number_of_points, cudaMemcpyDeviceToHost);
    cudaMemcpy(ex_cpu, ex_gpu, sizeof(float)*number_of_points, cudaMemcpyDeviceToHost);
    cudaMemcpy(y2a_cpu, y2a_gpu, sizeof(float)*number_of_points, cudaMemcpyDeviceToHost);

    //print_accumulate_point_data (NULL, NULL, NULL, ex_cpu, y2a_cpu, NULL, NULL, number_of_points);

    /////////////////////////////
    // FREE MEMORY
    cudaFree(ec_gpu);
    cudaFree(ex_gpu);
    cudaFree(y2a_gpu);
    cudaFree(accumulated_density_gpu);
    cudaFree(contracted_gradient_gpu);
    cudaFree(dxyz_gpu);
    cudaFree(dd1_gpu);
    cudaFree(dd2_gpu);

    free(ec_cpu);
    free(ex_cpu);
    free(y2a_cpu);
}


/////////////////////////////
// Proxy test
template <class T>
void proxyTest0001d()
{
#if FULL_DOUBLE
    printf("proxyTest0001() \n");

    ////////////////////////////////
    // PARAMS SETUP
    
    int number_of_points[9] = {221,227,256,537,1796,4007,2910,2910,3492};
    T* dens_cpu[9] = {dens_221,dens_227,dens_256,dens_537,dens_1796,dens_4007,dens_2910_1,dens_2910_2,dens_3492};
    T* contracted_gradients_cpu[9] = {contracted_grad_221,contracted_grad_227,contracted_grad_256,contracted_grad_537,contracted_grad_1796,contracted_grad_4007,contracted_grad_2910_1,contracted_grad_2910_2,contracted_grad_3492};
    G2G::vec_type<T,4>* grads[9] = {grad_221,grad_227,grad_256,grad_537,grad_1796,grad_4007,grad_2910_1,grad_2910_2,grad_3492};
    G2G::vec_type<T,4>* hess1s[9] = {hess1_221,hess1_227,hess1_256,hess1_537,hess1_1796,hess1_4007,hess1_2910_1,hess1_2910_2,hess1_3492};
    G2G::vec_type<T,4>* hess2s[9] = {hess2_221,hess2_227,hess2_256,hess2_537,hess2_1796,hess2_4007,hess2_2910_1,hess2_2910_2,hess2_3492};
    
    /*
    int number_of_points[1] = {10};
    T* dens_cpu[1] = {dens_10};
    T* contracted_gradients_cpu[1] = {contracted_grad_10};
    G2G::vec_type<T,4>* grads[1] = {grad_10};,grad_10};
    G2G::vec_type<T,4>* hess1s[1] = {hess1_10};
    G2G::vec_type<T,4>* hess2s[1] = {hess2_10};
    */

    for (int i=0; i<9; i++) {
        doGGA_gpu (number_of_points[i], 
	    dens_cpu[i], 
	    contracted_gradients_cpu[i],
	    grads[i],
	    hess1s[i],
	    hess2s[i]);
    }
#endif
}

template <class T>
void proxyTest0001f()
{
    printf("proxyTest0001() \n");

    ////////////////////////////////
    // PARAMS SETUP

    int number_of_points[9] = {221,227,256,537,1796,4007,2910,2910,3492};
    T* dens_cpu[9] = {dens_221_f,dens_227_f,dens_256_f,dens_537_f,dens_1796_f,dens_4007_f,dens_2910_1_f,dens_2910_2_f,dens_3492_f};
    T* contracted_gradients_cpu[9] = {contracted_grad_221_f,contracted_grad_227_f,contracted_grad_256_f,contracted_grad_537_f,contracted_grad_1796_f,contracted_grad_4007_f,contracted_grad_2910_1_f,contracted_grad_2910_2_f,contracted_grad_3492_f};
    G2G::vec_type<T,4>* grads[9] = {grad_221_f,grad_227_f,grad_256_f,grad_537_f,grad_1796_f,grad_4007_f,grad_2910_1_f,grad_2910_2_f,grad_3492_f};
    G2G::vec_type<T,4>* hess1s[9] = {hess1_221_f,hess1_227_f,hess1_256_f,hess1_537_f,hess1_1796_f,hess1_4007_f,hess1_2910_1_f,hess1_2910_2_f,hess1_3492_f};
    G2G::vec_type<T,4>* hess2s[9] = {hess2_221_f,hess2_227_f,hess2_256_f,hess2_537_f,hess2_1796_f,hess2_4007_f,hess2_2910_1_f,hess2_2910_2_f,hess2_3492_f};

    for (int i=0; i<9; i++) {
        doGGA_gpu_float (number_of_points[i], 
	    dens_cpu[i], 
	    contracted_gradients_cpu[i],
	    grads[i],
	    hess1s[i],
	    hess2s[i]);
    }
}


//////////////////////////////
// Test for FLOATS values
void proxyTest0002()
{
    printf("proxyTest0002() - FLOATS VERSION\n");

    ////////////////////////////////
    // PARAMS SETUP
    int number_of_points[9] = {221,227,256,537,1796,4007,2910,2910,3492};
    float* dens_cpu[9] = {dens_221_f,dens_227_f,dens_256_f,dens_537_f,dens_1796_f,dens_4007_f,dens_2910_1_f,dens_2910_2_f,dens_3492_f};
    float* contracted_gradients_cpu[9] = {contracted_grad_221_f,contracted_grad_227_f,contracted_grad_256_f,contracted_grad_537_f,contracted_grad_1796_f,contracted_grad_4007_f,contracted_grad_2910_1_f,contracted_grad_2910_2_f,contracted_grad_3492_f};
    G2G::vec_type<float,4>* grads[9]  = {grad_221_f,grad_227_f,grad_256_f,grad_537_f,grad_1796_f,grad_4007_f,grad_2910_1_f,grad_2910_2_f,grad_3492_f};
    G2G::vec_type<float,4>* hess1s[9] = {hess1_221_f,hess1_227_f,hess1_256_f,hess1_537_f,hess1_1796_f,hess1_4007_f,hess1_2910_1_f,hess1_2910_2_f,hess1_3492_f};
    G2G::vec_type<float,4>* hess2s[9] = {hess2_221_f,hess2_227_f,hess2_256_f,hess2_537_f,hess2_1796_f,hess2_4007_f,hess2_2910_1_f,hess2_2910_2_f,hess2_3492_f};

    for (int i=0; i<9; i++) {
        doGGA_gpu_float (number_of_points[i], 
	    dens_cpu[i], 
	    contracted_gradients_cpu[i],
	    grads[i],
	    hess1s[i],
	    hess2s[i]);
    }
}

//////////////////////////////////
// Conversion KERNELS
__global__ void floatToDouble(float* input, double* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	output[i] = (double)input[i];
    }
}

__global__ void doubleToFloat(double* input, float* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	output[i] = (float)input[i];
    }
}

__global__ void floatToDouble(const G2G::vec_type<float,4>* input, G2G::vec_type<double,4>* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	//float x, y, z, _w;
	output[i].x = (double)(input[i].x);
	output[i].y = (double)(input[i].y);
	output[i].z = (double)(input[i].z);
	//output[i].w = (double)input[i]._w;
    }
}

__global__ void doubleToFloat(const G2G::vec_type<double,4>* input, G2G::vec_type<float,4>* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	//float x, y, z, _w;
	output[i].x = (float)(input[i].x);
	output[i].y = (float)(input[i].y);
	output[i].z = (float)(input[i].z);
	//output[i].w = (float)input[i]._w;
    }
}

////////////////////////////////////////
// Convertion from double to float TEST
void conversionTest0001(int array_size) 
{
    printf("convertionTest0001(%i)\n", array_size);
    double* input = (double*)malloc(sizeof(double)*array_size);
    float* output = (float*)malloc(sizeof(float)*array_size);

    for (int i=0; i<array_size; i++) {
	input[i] = 0.000000001*i;
	output[i] = 0;
    }

    // CUDA ALLOC
    double* input_gpu = NULL;
    float* output_gpu = NULL;

    cudaMalloc((void**)&input_gpu, sizeof(double)*array_size);
    cudaMalloc((void**)&output_gpu, sizeof(float)*array_size);

    // CUDASET ARRAYS VALUEs
    cudaMemset(input_gpu, 0, sizeof(double)*array_size);
    cudaMemset(output_gpu, 0, sizeof(float)*array_size);
    cudaMemcpy(input_gpu, input, sizeof(double)*array_size, cudaMemcpyHostToDevice);

    // KERNEL CALL
    int threadsPerBlock = 256;
    int blocksPerGrid = (array_size + threadsPerBlock - 1) / threadsPerBlock;

    printf("double->float\n");
    doubleToFloat <<<blocksPerGrid, threadsPerBlock>>>(input_gpu, output_gpu, array_size);

    // Show the results
    cudaMemcpy (output, output_gpu, sizeof(float)*array_size, cudaMemcpyDeviceToHost);

    print_array (output, array_size);

    // Free memory
    cudaFree(input_gpu);
    cudaFree(output_gpu);
    free(input);
    free(output);
}

//////////////////////////////////////////
// Convertion from float to double TEST
void conversionTest0002(int array_size) 
{
    printf("convertionTest0002(%i)\n", array_size);
    double* output = (double*)malloc(sizeof(double)*array_size);
    float* input = (float*)malloc(sizeof(float)*array_size);

    for (int i=0; i<array_size; i++) {
	input[i] = 0.000000001*i;
	output[i] = 0;
    }

    // CUDA ALLOC
    float* input_gpu = NULL;
    double* output_gpu = NULL;

    cudaMalloc((void**)&input_gpu, sizeof(float)*array_size);
    cudaMalloc((void**)&output_gpu, sizeof(double)*array_size);

    // CUDASET ARRAYS VALUES
    cudaMemset(input_gpu, 0, sizeof(float)*array_size);
    cudaMemset(output_gpu, 0, sizeof(double)*array_size);
    cudaMemcpy(input_gpu, input, sizeof(float)*array_size, cudaMemcpyHostToDevice);

    // KERNEL CALL
    int threadsPerBlock = 256;
    int blocksPerGrid = (array_size + threadsPerBlock - 1) / threadsPerBlock;

    printf("float->double\n");
    floatToDouble <<<blocksPerGrid, threadsPerBlock>>>(input_gpu, output_gpu, array_size);

    // Show the results
    cudaMemcpy (output, output_gpu, sizeof(double)*array_size, cudaMemcpyDeviceToHost);

    print_array (output, array_size);

    // Free memory
    cudaFree(input_gpu);
    cudaFree(output_gpu);
    free(input);
    free(output);
}

///////////////////////////////////////////////////////
// Convertion from float to double for vec_type TEST
void conversionTest0003(int array_size) 
{
    printf("convertionTest0003(%i)\n", array_size);
    G2G::vec_type<double,4>* output = (G2G::vec_type<double,4>*)malloc(sizeof(G2G::vec_type<double,4>)*array_size);
    G2G::vec_type<float,4>* input = (G2G::vec_type<float,4>*)malloc(sizeof(G2G::vec_type<float,4>)*array_size);

    for (int i=0; i<array_size; i++) {
	input[i].x  = 0.000000001*i;
	input[i].y  = 0.000000002*i;
	input[i].z  = 0.000000004*i;
	input[i].w = 0.000000008*i;
    	output[i].x = 0;
	output[i].y = 0;
	output[i].z = 0;
	output[i].w = 0;
    }

    // CUDA ALLOC
    G2G::vec_type<float,4>* input_gpu = NULL;
    G2G::vec_type<double,4>* output_gpu = NULL;

    cudaMalloc((void**)&input_gpu, sizeof(G2G::vec_type<float,4>)*array_size);
    cudaMalloc((void**)&output_gpu, sizeof(G2G::vec_type<double,4>)*array_size);

    // CUDASET ARRAYS VALUES
    cudaMemset(input_gpu, 0, sizeof(G2G::vec_type<float,4>)*array_size);
    cudaMemset(output_gpu, 0, sizeof(G2G::vec_type<double,4>)*array_size);
    cudaMemcpy(input_gpu, input, sizeof(G2G::vec_type<float,4>)*array_size, cudaMemcpyHostToDevice);

    // KERNEL CALL
    int threadsPerBlock = 256;
    int blocksPerGrid = (array_size + threadsPerBlock - 1) / threadsPerBlock;

    printf("float->double\n");
    floatToDouble <<<blocksPerGrid, threadsPerBlock>>>(input_gpu, output_gpu, array_size);

    // Show the results
    cudaMemcpy (output, output_gpu, sizeof(G2G::vec_type<double,4>)*array_size, cudaMemcpyDeviceToHost);

    print_vec_type (output, array_size);

    // Free memory
    cudaFree((void*)input_gpu);
    cudaFree((void*)output_gpu);
    free((void*)input);
    free((void*)output);
}

///////////////////////////////////////////////////////
// Convertion from double to float for vec_type TEST
void conversionTest0004(int array_size) 
{
    printf("convertionTest0004(%i)\n", array_size);
    G2G::vec_type<double,4>* input = (G2G::vec_type<double,4>*)malloc(sizeof(G2G::vec_type<double,4>)*array_size);
    G2G::vec_type<float,4>* output = (G2G::vec_type<float,4>*)malloc(sizeof(G2G::vec_type<float,4>)*array_size);

    for (int i=0; i<array_size; i++) {
	input[i].x  = 0.000000001*i;
	input[i].y  = 0.000000002*i;
	input[i].z  = 0.000000004*i;
	input[i].w = 0.000000008*i;
    	output[i].x = 0;
	output[i].y = 0;
	output[i].z = 0;
	output[i].w = 0;
    }

    // CUDA ALLOC
    G2G::vec_type<double,4>* input_gpu = NULL;
    G2G::vec_type<float,4>* output_gpu = NULL;

    cudaMalloc((void**)&output_gpu, sizeof(G2G::vec_type<float,4>)*array_size);
    cudaMalloc((void**)&input_gpu, sizeof(G2G::vec_type<double,4>)*array_size);

    // CUDASET ARRAYS VALUES
    cudaMemset(output_gpu, 0, sizeof(G2G::vec_type<float,4>)*array_size);
    cudaMemset(input_gpu, 0, sizeof(G2G::vec_type<double,4>)*array_size);
    cudaMemcpy(input_gpu, input, sizeof(G2G::vec_type<double,4>)*array_size, cudaMemcpyHostToDevice);

    // KERNEL CALL
    int threadsPerBlock = 256;
    int blocksPerGrid = (array_size + threadsPerBlock - 1) / threadsPerBlock;

    printf("double->float\n");
    doubleToFloat <<<blocksPerGrid, threadsPerBlock>>>(input_gpu, output_gpu, array_size);

    // Show the results
    cudaMemcpy (output, output_gpu, sizeof(G2G::vec_type<float,4>)*array_size, cudaMemcpyDeviceToHost);

    print_vec_type (output, array_size);

    // Free memory
    cudaFree((void*)input_gpu);
    cudaFree((void*)output_gpu);
    free((void*)input);
    free((void*)output);
}

void data_type_test003()
{
    printf("=== test 003 ===\n");
    float f_value = 1;
    double d_value = 2;
    int i_value = 3;
    //const std::type_info& float_type = typeid(f_value);
    //const std::type_info& double_type = typeid(d_value);

    printf("%s\n", typeid(f_value).name());
    printf("%s\n", typeid(d_value).name());
    printf("%s\n", typeid(i_value).name());
    printf("%i\n", sizeof(f_value));
    printf("%i\n", sizeof(d_value));
    printf("%i\n", sizeof(i_value));

    if (typeid(f_value).name() == "f") {
	printf("Soy un float \n");
    }
}

void pointer_proxy_test0001()
{
    printf("pointer_proxy_test0001()\n");
    LibxcProxy<double,4>* proxy_pointer = NULL;
    proxy_pointer = new LibxcProxy<double,4> (1,1,1);
    free(proxy_pointer);
}

void init_proxy_test0001()
{
    printf("init_proxy_test0001()\n");
    LibxcProxy<float,4> proxy;

    // Create the libxcproxy
    const int nspin = 1;
    const int functionalExchange = 1101;
    const int functionalCorrelation = 1130;

    proxy.init(functionalExchange, functionalCorrelation, nspin);

    proxy.printFunctionalsInformation (functionalExchange, functionalCorrelation);
    //proxy.closeProxy();
}

void printFunctionalsInformationTest001 () 
{
    int nspin = 1;
    int functionalExchange = 1101;
    int functionalCorrelation = 1130;

    LibxcProxy<double,3> aProxy;

    aProxy.printFunctionalsInformation (functionalExchange, functionalCorrelation);

}

void printFunctionalsInformationTest002 () 
{
    int nspin = 1;
    int functionalExchange = 101;
    int functionalCorrelation = 130;

    LibxcProxy<double,3> aProxy;

    aProxy.printFunctionalsInformation (functionalExchange, functionalCorrelation);

}


int main()
{
    printf("Test: Libxc Proxy GPU - BEGIN\n");
    try {
        accumulate_data_for_libxc_test0001();
	joinResultsTest0001();
#if FULL_DOUBLE
        proxyTest0001d<double>();
#else
	proxyTest0001f<float>();
#endif
        proxyTest0002();
        conversionTest0001(100);
        conversionTest0002(100);
        conversionTest0003(10);
        conversionTest0004(10);
        data_type_test003();
	pointer_proxy_test0001();
        init_proxy_test0001 ();
	//printFunctionalsInformationTest001 ();
	//printFunctionalsInformationTest002 ();
    } catch (int e) {
	printf("An exception occurred: %u \n", e);
	exit (EXIT_FAILURE);
    }
    printf("Test: Libxc Proxy GPU - END\n");
    return 0;
}

