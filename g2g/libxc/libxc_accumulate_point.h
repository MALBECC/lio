#define WIDTH 4

#include "libxcproxy.h"
#include "../timer.h"

extern "C" void g2g_timer_sum_start_(const char* timer_name, unsigned int length_arg);
extern "C" void g2g_timer_sum_stop_(const char* timer_name, unsigned int length_arg);
extern "C" void g2g_timer_sum_pause_(const char* timer_name, unsigned int length_arg);

template<class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_accumulate_point_for_libxc(T* const point_weights,
                uint points, int block_height, 
		T* partial_density_in,
		G2G::vec_type<T,WIDTH>* dxyz_in,
                G2G::vec_type<T,WIDTH>* dd1_in,
		G2G::vec_type<T,WIDTH>* dd2_in,
		T* accumulated_density,
		G2G::vec_type<T,WIDTH>* dxyz_accum,
                G2G::vec_type<T,WIDTH>* dd1_accum,
		G2G::vec_type<T,WIDTH>* dd2_accum) {

  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;
  //uint point = blockIdx.x * blockDim.x + threadIdx.x;

  T _partial_density(0.0f);
  G2G::vec_type<T,WIDTH> _dxyz, _dd1, _dd2;
  _dxyz = _dd1 = _dd2 = G2G::vec_type<T,WIDTH>(0.0f,0.0f,0.0f,0.0f);

  bool valid_thread = (point < points);

  if (valid_thread) {
    for(int j =0 ; j<block_height; j++) {
      const int this_row = j*points+point;
      _partial_density += partial_density_in[this_row];
      _dxyz += dxyz_in[this_row];
      _dd1 += dd1_in[this_row];
      _dd2 += dd2_in[this_row];
    }

    // Accumulate the data for libxc.
    accumulated_density[point] = _partial_density;
    /*printf("point:%i %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", point,
	_dxyz.x, _dxyz.y, _dxyz.z,
	_dd1.x, _dd1.y, _dd1.z,
	_dd2.x, _dd2.y, _dd2.z);*/
    
    dxyz_accum[point] = _dxyz;
    dd1_accum[point] = _dd1;
    dd2_accum[point] = _dd2;
    
  }

}

template<class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_accumulate_energy_and_forces_from_libxc (T* const energy, 
		    T* const factor, 
		    T* const point_weights, 
		    uint points, 
		    T* accumulated_density)
{
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;
  //uint point = blockIdx.x * blockDim.x + threadIdx.x;

  T point_weight = 0.0f;

  bool valid_thread = (point < points);
  if (valid_thread) {
    point_weight = point_weights[point];
  }

  if (compute_energy && valid_thread && energy != NULL) {
        energy[point] *= (accumulated_density[point] * point_weight);
  }

  if (compute_factor && valid_thread && factor != NULL) {
    factor[point] *= point_weight;
  }

}

template<class T, bool compute_energy, bool compute_factor, bool lda>
    void libxc_exchange_correlation_cpu (LibxcProxy<T, WIDTH>* libxcProxy,
	T* energy_gpu,
	T* factor_gpu,
	uint points,
	T* accumulated_density_gpu,
	G2G::vec_type<T,WIDTH>* dxyz_gpu_accum,
        G2G::vec_type<T,WIDTH>* dd1_gpu_accum,
	G2G::vec_type<T,WIDTH>* dd2_gpu_accum)
{
    //printf("libxc_exchage_correlation_cpu (...) \n");

    cudaError_t err = cudaSuccess;

    // Copy to host the matrix data in gpu memory and
    // call the new libxcProxy.
    G2G::vec_type<T,WIDTH>* dxyz_cpu;
    G2G::vec_type<T,WIDTH>* dd1_cpu;
    G2G::vec_type<T,WIDTH>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = points * sizeof(G2G::vec_type<T,WIDTH>);
    dxyz_cpu = (G2G::vec_type<T,WIDTH> *)malloc(size);
    dd1_cpu = (G2G::vec_type<T,WIDTH> *)malloc(size);
    dd2_cpu = (G2G::vec_type<T,WIDTH> *)malloc(size);

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu_accum, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        //printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu_accum, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        //printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu_accum, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        //printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Allocate the host input vectors
    uint array_size = sizeof(T)*points;
    T *energy_cpu = NULL;
    if (compute_energy) 
    { 
	energy_cpu = (T *)malloc(array_size);
    }
    
    T *factor_cpu = NULL;
    if (compute_factor) 
    {
	factor_cpu = (T *)malloc(array_size);
    }

    T *accumulated_density_cpu = (T*)malloc(array_size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    if (compute_energy) 
    {
	err = cudaMemcpy(energy_cpu, energy_gpu, array_size, cudaMemcpyDeviceToHost);
        if (err != cudaSuccess)
	{
    	    //printf("Failed to copy vector energy_gpu from device to host!\n");
    	    exit(EXIT_FAILURE);
	}
    }

    if (compute_factor) 
    {
	err = cudaMemcpy(factor_cpu, factor_gpu, array_size, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess)
	{
    	    //printf("Failed to copy vector factor_gpu from device to host!\n");
    	    exit(EXIT_FAILURE);
	}
    }

    err = cudaMemcpy(accumulated_density_cpu, accumulated_density_gpu, array_size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
	//printf("Failed to copy vector accumulated_density_gpu to host!\n");
    }

    /////////////////////
    // Libxc proxy
    /////////////////////
    /* This is old. Remove later. 08/02/2018
    T y2a, exc_x, exc_c = 0;
    if (libxcProxy != NULL) {
	for (int i=0; i<points; i++) {
    	    libxcProxy->doGGA (accumulated_density_cpu[i], dxyz_cpu[i], dd1_cpu[i], dd2_cpu[i], exc_x, exc_c, y2a);
	    if (compute_energy) 
	    {
		energy_cpu[i] = exc_x + exc_c;
	    }
	    if (compute_factor)
	    {
		factor_cpu[i] = y2a;
	    }
	}
    }
    */

    // Parameters
    T* exc;
    T* corr;
    T* y2a;

    // Now alloc memory for the data
    //exc  = (T*)malloc(sizeof(T)*points);
    //corr = (T*)malloc(sizeof(T)*points);
    //y2a  = (T*)malloc(sizeof(T)*points);

    exc  = (T*)malloc(array_size);
    corr = (T*)malloc(array_size);
    y2a  = (T*)malloc(array_size);

    if (libxcProxy != NULL) {
        libxcProxy->doGGA (accumulated_density_cpu, points, dxyz_cpu, dd1_cpu, dd2_cpu, exc, corr, y2a);
	// Join the results.
	for (unsigned int i=0; i<points; i++) {
	    if (compute_energy) 
	    {
		energy_cpu[i] = exc[i] + corr[i];
	    }
	    if (compute_factor)
	    {
		factor_cpu[i] = y2a[i];
	    }
	}
    }

    // Now copy back the results to the gpu.
    if (compute_energy) {
	err = cudaMemcpy(energy_gpu, energy_cpu, array_size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
    	    //printf("Failed to copy vector energy_gpu from host to device!\n");
    	    exit(EXIT_FAILURE);
	}
    }

    if (compute_factor) {
	err = cudaMemcpy(factor_gpu, factor_cpu, array_size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
    	    //printf("Failed to copy vector factor_cpu from host to device!\n");
    	    exit(EXIT_FAILURE);
	}
    }

    // Free memory.
    if (compute_energy) {
        free(energy_cpu);
    }
    if (compute_factor) {
	free(factor_cpu);
    }

    free(accumulated_density_cpu);
    free(exc);
    free(corr);
    free(y2a);

    free(dd1_cpu);
    free(dd2_cpu);
    free(dxyz_cpu);
}

template<class T>
__global__ void calculateContractedGradient(G2G::vec_type<T,WIDTH>* grad, T* contracted_grad, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        contracted_grad[i] = (grad[i].x * grad[i].x) + (grad[i].y * grad[i].y) + (grad[i].z * grad[i].z);
    }
}

template<class T>
__global__ void vectorAdd(const T* A, const T* B, T* C, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        C[i] = A[i] + B[i];
    }
}

template<class T, bool compute_energy, bool compute_factor, bool lda>
    void libxc_exchange_correlation_gpu (LibxcProxy<T, WIDTH>* libxcProxy,
	T* energy_gpu,
	T* factor_gpu,
	uint points,
	T* accumulated_density_gpu,
	G2G::vec_type<T,WIDTH>* dxyz_gpu,
        G2G::vec_type<T,WIDTH>* dd1_gpu,
	G2G::vec_type<T,WIDTH>* dd2_gpu)
{
    //printf("libxc_exchage_correlation_gpu (...) \n");

    cudaError_t err = cudaSuccess;

    /////////////////////
    // Libxc proxy
    /////////////////////
    // Parameters
    T* exc_gpu = NULL;
    T* corr_gpu = NULL;
    T* y2a_gpu = NULL;
    T* contracted_gradient = NULL;
    int array_size = sizeof(T) * points;

    //g2g_timer_sum_start_("libxc_exchange_correlation_gpu_preparing_data",45);

    // Alloc memory for cuda variables.
    err = cudaMalloc((void **)&exc_gpu, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device exc_gpu!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&corr_gpu, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device corr_gpu!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&y2a_gpu, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device y2a_gpu!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void**)&contracted_gradient, sizeof(T)*points);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device contracted_gradient!\n");
	exit(EXIT_FAILURE);
    }

    cudaMemset(exc_gpu,0,array_size);
    cudaMemset(corr_gpu,0,array_size);
    cudaMemset(y2a_gpu,0,array_size);
    cudaMemset(contracted_gradient,0,array_size);

    // Variables for the Kernels
    int threadsPerBlock = 256;
    int blocksPerGrid = (points + threadsPerBlock - 1) / threadsPerBlock;

    if (libxcProxy != NULL) {
        //printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

	// Call the Kernel to compute the contracted gradient
        //printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
	calculateContractedGradient<<<blocksPerGrid, threadsPerBlock>>>(dxyz_gpu, contracted_gradient, points);

	//g2g_timer_sum_pause_("libxc_exchange_correlation_gpu_preparing_data",45);

	/////////////////////////////////////////
	// Compute exc_corr using LIBXC GPU.
	libxcProxy->doGGA (accumulated_density_gpu,
	    points, 
	    contracted_gradient, 
	    dxyz_gpu, 
	    dd1_gpu, 
	    dd2_gpu, 
	    exc_gpu, 
	    corr_gpu, 
	    y2a_gpu);

	///////////////////////
	// Join the results.
	//g2g_timer_sum_start_("libxc_exchange_correlation_gpu_joining_data",43);

    	if (compute_energy && energy_gpu != NULL)
	{
	    vectorAdd<<<blocksPerGrid, threadsPerBlock>>>(exc_gpu, corr_gpu, energy_gpu, points);
	}
	if (compute_factor && factor_gpu != NULL)
	{
	    err = cudaMemcpy(factor_gpu, y2a_gpu, array_size, cudaMemcpyDeviceToDevice);
	    if (err != cudaSuccess)
	    {
    		fprintf(stderr, "Failed to copy y2a_gpu to factor_gpu from DeviceToDevice for array of %i !\n", array_size);
	    }
	}
	//g2g_timer_sum_pause_("libxc_exchange_correlation_gpu_joining_data",43);
    }

    // Free memory.
    if (contracted_gradient != NULL) {
	cudaFree(contracted_gradient);
    }
    if (exc_gpu != NULL) {
	cudaFree(exc_gpu);
    }
    if (corr_gpu != NULL) {
        cudaFree(corr_gpu);
    }
    if (y2a_gpu != NULL) {
	cudaFree(y2a_gpu);
    }
}

