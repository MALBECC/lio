#include "../scalar_vector_types.h"

//////////////////////////////////////
//// PRINT UTILS

template<class scalar_type>
void print_array (scalar_type* data, int size) 
{
    printf ("[");
    if (data == NULL) {
	printf("empty");
    } else {
	for (int i=0; i<size; i++) {
	    printf("%e,", data[i]);
	}
    }
    printf("]\n");
}

template<class scalar_type>
void print_vec_type (G2G::vec_type<scalar_type,4>* data, int size) 
{
    printf ("[");
    if (data == NULL) {
	printf("empty");
    } else {
	for (int i=0; i<size; i++) {
	    printf("(%e,%e,%e),", data[i].x, data[i].y, data[i].z);
	}
    }
    printf("]\n");
}

template<class scalar_type> void print_accumulate_parameters (
    scalar_type* energy_gpu, 
    scalar_type* factor_gpu, 
    scalar_type* point_weights_gpu,
    int number_of_points,
    int block_height,
    scalar_type* partial_densities_gpu,
    G2G::vec_type<scalar_type,4>* dxyz_gpu,
    G2G::vec_type<scalar_type,4>* dd1_gpu,
    G2G::vec_type<scalar_type,4>* dd2_gpu)
{
    // Bajar la info a CPU y luego imprimirla.
    // Copy to host the matrix data in gpu memory and
    // call the new libxcProxy.
    G2G::vec_type<scalar_type,4>* dxyz_cpu;
    G2G::vec_type<scalar_type,4>* dd1_cpu;
    G2G::vec_type<scalar_type,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = number_of_points * sizeof(G2G::vec_type<scalar_type,4>);
    dxyz_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);
    dd1_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);
    dd2_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);

    // Copy data from device to host.
    if (dxyz_gpu != NULL)
        cudaMemcpy(dxyz_cpu, dxyz_gpu, size, cudaMemcpyDeviceToHost);
    if (dd1_gpu != NULL)
	cudaMemcpy(dd1_cpu, dd1_gpu, size, cudaMemcpyDeviceToHost);
    if (dd2_gpu != NULL)
        cudaMemcpy(dd2_cpu, dd2_gpu, size, cudaMemcpyDeviceToHost);

    // Allocate the host input vectors
    uint array_size = number_of_points * sizeof(scalar_type);
    scalar_type *energy_cpu = (scalar_type *)malloc(array_size);
    scalar_type *factor_cpu = (scalar_type *)malloc(array_size);
    scalar_type *partial_densities_cpu = (scalar_type*)malloc(array_size);
    scalar_type *point_weights_cpu = (scalar_type*)malloc(array_size);
    if (energy_gpu != NULL) 
        cudaMemcpy(energy_cpu, energy_gpu, array_size, cudaMemcpyDeviceToHost);
    if (factor_gpu != NULL)
        cudaMemcpy(factor_cpu, factor_gpu, array_size, cudaMemcpyDeviceToHost);
    if (partial_densities_gpu != NULL)
	cudaMemcpy(partial_densities_cpu, partial_densities_gpu, array_size, cudaMemcpyDeviceToHost);
    if (point_weights_gpu != NULL)
	cudaMemcpy(point_weights_cpu, point_weights_gpu, array_size, cudaMemcpyDeviceToHost);

    printf("=========\n");
    printf("= Data  =\n");
    printf("=========\n");
    printf("number_of_points:%i\n", number_of_points);
    printf("block_height:%i\n", block_height);
    printf("dxyz:"); print_vec_type(dxyz_cpu, number_of_points);
    printf("dd1:"); print_vec_type(dd1_cpu, number_of_points);
    printf("dd2:"); print_vec_type(dd2_cpu, number_of_points);
    printf("energy:"); print_array(energy_cpu, number_of_points);
    printf("factor:"); print_array(factor_cpu, number_of_points);
    printf("point_weights:"); print_array(point_weights_cpu, number_of_points);
    printf("partial_densities:"); print_array(partial_densities_cpu, number_of_points);
    printf("=========\n");

    // Free memory.
    free(energy_cpu);
    free(factor_cpu);
    free(partial_densities_cpu);
    free(dd1_cpu);
    free(dd2_cpu);
    free(dxyz_cpu);

}

template <class scalar_type> 
void print_proxy_input (
    scalar_type* dens,
    const int number_of_points,
    const scalar_type* contracted_grad,
    const G2G::vec_type<scalar_type, 4>* grad,
    const G2G::vec_type<scalar_type, 4>* hess1,
    const G2G::vec_type<scalar_type, 4>* hess2)
{
    // Bajar la info a CPU y luego imprimirla.
    // Copy to host the matrix data in gpu memory and
    // call the new libxcProxy.
    scalar_type* dens_cpu;
    scalar_type* contracted_grad_cpu;
    G2G::vec_type<scalar_type,4>* grad_cpu;
    G2G::vec_type<scalar_type,4>* hess1_cpu;
    G2G::vec_type<scalar_type,4>* hess2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = number_of_points * sizeof(G2G::vec_type<scalar_type,4>);
    uint array_size = number_of_points * sizeof(scalar_type);
    dens_cpu = (scalar_type*)malloc(size);
    contracted_grad_cpu = (scalar_type*)malloc(size);
    grad_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);
    hess1_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);
    hess2_cpu = (G2G::vec_type<scalar_type,4> *)malloc(size);

    // Copy data from device to host.
    if (grad != NULL)
        cudaMemcpy(grad_cpu, grad, size, cudaMemcpyDeviceToHost);
    if (hess1 != NULL)
	cudaMemcpy(hess1_cpu, hess1, size, cudaMemcpyDeviceToHost);
    if (hess2 != NULL)
        cudaMemcpy(hess2_cpu, hess2, size, cudaMemcpyDeviceToHost);
    if (contracted_grad != NULL)
	cudaMemcpy(contracted_grad_cpu, contracted_grad, array_size, cudaMemcpyDeviceToHost);
    if (dens != NULL)
	cudaMemcpy(dens_cpu, dens, array_size, cudaMemcpyDeviceToHost);

    printf("==============================\n");
    printf("= Proxy input for libxc_cuda =\n");
    printf("==============================\n");
    printf("number_of_points:%i\n", number_of_points);
    printf("grad:"); print_vec_type(grad_cpu, number_of_points);
    printf("hess1:"); print_vec_type(hess1_cpu, number_of_points);
    printf("hess2:"); print_vec_type(hess2_cpu, number_of_points);
    printf("dens:"); print_array(dens_cpu, number_of_points);
    printf("contracted_grad:"); print_array(contracted_grad_cpu, number_of_points);
    printf("==============================\n");

    // Free memory.
    free(dens_cpu);
    free(contracted_grad_cpu);
    free(hess2_cpu);
    free(hess1_cpu);
    free(grad_cpu);

}

