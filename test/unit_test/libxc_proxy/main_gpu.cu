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
    int numElements = number_of_points;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    uint block_height = 1;

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


int main()
{
    cout << "Test: Libxc Proxy GPU - BEGIN" << endl;
    accumulate_data_for_libxc_test0001();
    cout << "Test: Libxc Proxy GPU - END" << endl;
    return 0;
}

