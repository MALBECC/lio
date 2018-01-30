#define WIDTH 4

#include "libxcproxy.h"

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
 void libxc_cpu_accumulate_point(LibxcProxy<scalar_type, WIDTH>* libxcProxy, 
	scalar_type* energy, scalar_type* factor, scalar_type* point_weights,
        uint points, int block_height, scalar_type* partial_density, 
	G2G::vec_type<scalar_type,WIDTH>* dxyz,
        G2G::vec_type<scalar_type,WIDTH>* dd1, 
	G2G::vec_type<scalar_type,WIDTH>* dd2) {

  //uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

    scalar_type point_weight = 0.0f;
    scalar_type y2a, exc_corr, exc_c, exc_x;

    scalar_type _partial_density(0.0f);
    G2G::vec_type<scalar_type,WIDTH> _dxyz, _dd1, _dd2;

    _dxyz = _dd1 = _dd2 = G2G::vec_type<scalar_type,WIDTH>(0.0f,0.0f,0.0f,0.0f);

    for(int i=0; i<points; i++) {
	point_weight = point_weights[i];
	
	for(int j=0; j<block_height; j++) {
	    //const int this_row = j*points+point;
	    const int this_row = j*points+i;
	    _partial_density += partial_density[this_row];
	    _dxyz += dxyz[this_row];
	    _dd1 += dd1[this_row];
	    _dd2 += dd2[this_row];
	}
	
	// TODO: aca va la llamada al proxy.
	//calc_ggaCS_in<scalar_type, 4>(_partial_density, _dxyz, _dd1, _dd2, exc_x, exc_c, y2a, 9);
	libxcProxy->doGGA(_partial_density, _dxyz, _dd1, _dd2, exc_x, exc_c, y2a);

	exc_corr = exc_x + exc_c;

	if (compute_energy) {
	    energy[i] = (_partial_density * point_weight) * exc_corr;
	}

	if (compute_factor) {
	    factor[i] = point_weight * y2a;
	}
    }

//  bool valid_thread = (point < points);
//  if (valid_thread)
//    point_weight = point_weights[point];

//  if (valid_thread) {
//    for(int j =0 ; j<block_height; j++) {
//      const int this_row = j*points+point;
//      _partial_density += partial_density[this_row];
//      _dxyz += dxyz[this_row];
//      _dd1 += dd1[this_row];
//      _dd2 += dd2[this_row];
//     }
//  }
  // calc_ggaCS_in<scalar_type, 4>(_partial_density, _dxyz, _dd1, _dd2, exc_x, exc_c, y2a, 9);
  // TODO: aca va la llamada al proxy.
  // calc_ggaCS_in<scalar_type, 4>(_partial_density, _dxyz, _dd1, _dd2, exc_x, exc_c, y2a, 9);
//  exc_corr = exc_x + exc_c;
//  if (compute_energy && valid_thread)
//    energy[point] = (_partial_density * point_weight) * exc_corr;
//  if (compute_factor && valid_thread)
//    factor[point] = point_weight * y2a;

}

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda> 
    void cpu_accumulate_point (LibxcProxy<scalar_type, WIDTH>* libxcProxy, scalar_type* energy_gpu,
	scalar_type* factor_gpu, scalar_type* point_weights_gpu,
        uint points, int block_height, scalar_type* partial_density_gpu, 
	G2G::vec_type<scalar_type,WIDTH>* dxyz_gpu,
        G2G::vec_type<scalar_type,WIDTH>* dd1_gpu,
	G2G::vec_type<scalar_type,WIDTH>* dd2_gpu)
{
    printf("accumulate_point_cpu (...) \n");
    cudaError_t err = cudaSuccess;

    // Copy to host the matrix data in gpu memory and
    // call the new libxc_cpu_accumulate_point.
    G2G::vec_type<scalar_type,4>* dxyz_cpu;
    G2G::vec_type<scalar_type,4>* dd1_cpu;
    G2G::vec_type<scalar_type,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = points * sizeof(G2G::vec_type<scalar_type,4>);
    dxyz_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);
    dd1_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);
    dd2_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Allocate the host input vectors
    uint array_size = points * sizeof(scalar_type);
    scalar_type *energy_cpu = NULL;
    if (compute_energy) { 
	energy_cpu = (scalar_type *)malloc(array_size);
    }
    scalar_type *factor_cpu = NULL;
    if (compute_factor) {
	factor_cpu = (scalar_type *)malloc(array_size);
    }
    scalar_type *point_weights_cpu = (scalar_type *)malloc(array_size);
    scalar_type *partial_density_cpu = (scalar_type *)malloc(array_size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    if (compute_energy) {
	err = cudaMemcpy(energy_cpu, energy_gpu, array_size, cudaMemcpyDeviceToHost);
        if (err != cudaSuccess)
	{
    	    printf("Failed to copy vector energy_gpu from device to host!\n");
    	    exit(EXIT_FAILURE);
	}
    }

    if (compute_factor) {
	err = cudaMemcpy(factor_cpu, factor_gpu, array_size, cudaMemcpyDeviceToHost);
	if (err != cudaSuccess)
	{
    	    printf("Failed to copy vector factor_gpu from device to host!\n");
    	    exit(EXIT_FAILURE);
	}
    }

    err = cudaMemcpy(point_weights_cpu, point_weights_gpu, array_size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_cpu, partial_density_gpu, array_size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    /////////////////////
    // Accumulate CPU
    /////////////////////
    libxc_cpu_accumulate_point<scalar_type, compute_energy, compute_factor, lda>(libxcProxy,
	energy_cpu, factor_cpu, point_weights_cpu,
	points, block_height, partial_density_cpu, 
	dxyz_cpu, dd1_cpu, dd2_cpu);

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

    err = cudaMemcpy(point_weights_gpu, point_weights_cpu, array_size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_gpu, partial_density_cpu, array_size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    // Copy the matrix data.
    err = cudaMemcpy(dxyz_gpu, dxyz_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dxyz_cpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_gpu, dd1_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd1_cpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_gpu, dd2_cpu, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        printf("Failed to copy dd2_cpu from host to device!\n");
        exit(EXIT_FAILURE);
    }

    // Free memory.
    free(energy_cpu);
    free(factor_cpu);
    free(point_weights_cpu);
    free(partial_density_cpu);

    free(dd1_cpu);
    free(dd2_cpu);
    free(dxyz_cpu);
}

