#include <iostream>
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

//#include "../../../g2g/pointxc/calc_ggaCS.h"
//#include "../../../g2g/pointxc/calc_ggaOS.h"

//#include "../../../g2g/cuda/kernels/accumulate_point.h"

#include "../../../g2g/libxc/libxcproxy.h"
#include "../../../g2g/libxc/libxc_accumulate_point.h"


using namespace std;

using std::cout;
using std::endl;

void accumulate_data_for_libxc_test0001()
{
    printf("** accumulate_data_for_libxc_test0001 **\n");

    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in;

    // Accumulate
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_accum;

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);

    dxyz_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);

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
    cudaMemset(point_weights_gpu_in, 1, size);
    cudaMemset(partial_density_gpu_in, 1, size);
    cudaMemset(accumulated_density_gpu, 1, size);

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

}

template <class scalar_type, int width>
__global__ void funcionDeMierda(
		    double* ex, double* exchange,
		    double* ec, double* correlation,
		    double* vrho, double* vrhoC,
		    double* vsigma, double* vsigmaC,
		    double* v2rho, double* v2rhoC,
		    double* v2rhosigma, double* v2rhosigmaC,
		    double* v2sigma, double* v2sigmaC,
		    double* y2a,
		    double* sigma,
		    G2G::vec_type<double, width>* grad,
		    G2G::vec_type<double, width>* hess1,
		    G2G::vec_type<double, width>* hess2,
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

int main()
{
    cout << "Test: Libxc Proxy GPU - BEGIN" << endl;
    //accumulate_data_for_libxc_test0001();
    joinResultsTest0001();
    cout << "Test: Libxc Proxy GPU - END" << endl;
    return 0;
}

