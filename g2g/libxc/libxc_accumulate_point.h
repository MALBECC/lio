#define WIDTH 4

#include "libxcproxy.h"

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_accumulate_point_for_libxc(scalar_type* const point_weights,
                uint points, int block_height, 
		scalar_type* partial_density_in,
		G2G::vec_type<scalar_type,WIDTH>* dxyz_in,
                G2G::vec_type<scalar_type,WIDTH>* dd1_in,
		G2G::vec_type<scalar_type,WIDTH>* dd2_in,
		scalar_type* accumulated_density,
		G2G::vec_type<scalar_type,WIDTH>* dxyz_accum,
                G2G::vec_type<scalar_type,WIDTH>* dd1_accum,
		G2G::vec_type<scalar_type,WIDTH>* dd2_accum) {

  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  scalar_type _partial_density(0.0f);
  G2G::vec_type<scalar_type,WIDTH> _dxyz, _dd1, _dd2;
  _dxyz = _dd1 = _dd2 = G2G::vec_type<scalar_type,WIDTH>(0.0f,0.0f,0.0f,0.0f);

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
    dxyz_accum[point] = _dxyz;
    dd1_accum[point] = _dd1;
    dd2_accum[point] = _dd2;

  }

}

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_accumulate_energy_and_forces_from_libxc (scalar_type* const energy, 
		    scalar_type* const factor, 
		    scalar_type* const point_weights, 
		    uint points, 
		    scalar_type* accumulated_density)
{
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  scalar_type point_weight = 0.0f;

  bool valid_thread = (point < points);
  if (valid_thread)
    point_weight = point_weights[point];

  if (compute_energy && valid_thread)
    energy[point] *= (accumulated_density[point] * point_weight);

  if (compute_factor && valid_thread)
    factor[point] *= point_weight;

}

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
    void libxc_exchange_correlation_cpu (LibxcProxy<scalar_type, WIDTH>* libxcProxy,
	scalar_type* energy_gpu,
	scalar_type* factor_gpu,
	uint points,
	scalar_type* accumulated_density_gpu,
	G2G::vec_type<scalar_type,WIDTH>* dxyz_gpu_accum,
        G2G::vec_type<scalar_type,WIDTH>* dd1_gpu_accum,
	G2G::vec_type<scalar_type,WIDTH>* dd2_gpu_accum)
{
    printf("libxc_exchage_correlation_cpu (...) \n");

    cudaError_t err = cudaSuccess;

    // Copy to host the matrix data in gpu memory and
    // call the new libxcProxy.
    G2G::vec_type<scalar_type,4>* dxyz_cpu;
    G2G::vec_type<scalar_type,4>* dd1_cpu;
    G2G::vec_type<scalar_type,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = points * sizeof(G2G::vec_type<scalar_type,4>);
    dxyz_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);
    dd1_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);
    dd2_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu_accum, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu_accum, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu_accum, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Allocate the host input vectors
    uint array_size = points * sizeof(scalar_type);
    scalar_type *energy_cpu = NULL;
    if (compute_energy) 
    { 
	energy_cpu = (scalar_type *)malloc(array_size);
    }
    
    scalar_type *factor_cpu = NULL;
    if (compute_factor) 
    {
	factor_cpu = (scalar_type *)malloc(array_size);
    }

    scalar_type *accumulated_density_cpu = (scalar_type*)malloc(array_size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    if (compute_energy) 
    {
	err = cudaMemcpy(energy_cpu, energy_gpu, array_size, cudaMemcpyDeviceToHost);
        if (err != cudaSuccess)
	{
    	    printf("Failed to copy vector energy_gpu from device to host!\n");
    	    exit(EXIT_FAILURE);
	}
    }

    if (compute_factor) 
    {
	err = cudaMemcpy(factor_cpu, factor_gpu, array_size, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess)
	{
    	    printf("Failed to copy vector factor_gpu from device to host!\n");
    	    exit(EXIT_FAILURE);
	}
    }

    err = cudaMemcpy(accumulated_density_cpu, accumulated_density_gpu, array_size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
	printf("Failed to copy vector accumulated_density_gpu to host!\n");
    }

    /////////////////////
    // Libxc proxy
    /////////////////////
    scalar_type y2a, exc_x, exc_c = 0;
    if (libxcProxy != NULL) {
	for (int i=0; i<points; i++) {
	    // TODO: try to modify libxProxy to use the 
	    //full array instead of having a cycle.
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

    // Now copy back the results to the gpu.
    if (compute_energy) {
	err = cudaMemcpy(energy_gpu, energy_cpu, array_size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
    	    printf("Failed to copy vector energy_gpu from host to device!\n");
    	    exit(EXIT_FAILURE);
	}
    }

    if (compute_factor) {
	err = cudaMemcpy(factor_gpu, factor_cpu, array_size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
    	    printf("Failed to copy vector factor_cpu from host to device!\n");
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

    free(dd1_cpu);
    free(dd2_cpu);
    free(dxyz_cpu);
}
