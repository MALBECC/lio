#ifndef __LIBXCPROXY_H_PRINT_UTILS__
#define __LIBXCPROXY_H_PRINT_UTILS__

#include "../scalar_vector_types.h"

//////////////////////////////////////
//// PRINT UTILS
template <class T>
void print_array (T* data, int size) 
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

template <class T>
void print_vec_type (G2G::vec_type<T,4>* data, int size) 
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

template <class T> 
void print_accumulate_parameters (
    T* energy_gpu, 
    T* factor_gpu, 
    T* point_weights_gpu,
    int number_of_points,
    int block_height,
    T* partial_densities_gpu,
    G2G::vec_type<T,4>* dxyz_gpu,
    G2G::vec_type<T,4>* dd1_gpu,
    G2G::vec_type<T,4>* dd2_gpu)
{
#ifdef __CUDACC__
    // Bajar la info a CPU y luego imprimirla.
    // Copy to host the matrix data in gpu memory and
    // call the new libxcProxy.
    G2G::vec_type<T,4>* dxyz_cpu;
    G2G::vec_type<T,4>* dd1_cpu;
    G2G::vec_type<T,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = number_of_points * sizeof(G2G::vec_type<T,4>);
    dxyz_cpu = (G2G::vec_type<T,4> *)malloc(size);
    dd1_cpu = (G2G::vec_type<T,4> *)malloc(size);
    dd2_cpu = (G2G::vec_type<T,4> *)malloc(size);

    // Copy data from device to host.
    if (dxyz_gpu != NULL)
        cudaMemcpy(dxyz_cpu, dxyz_gpu, size, cudaMemcpyDeviceToHost);
    if (dd1_gpu != NULL)
	cudaMemcpy(dd1_cpu, dd1_gpu, size, cudaMemcpyDeviceToHost);
    if (dd2_gpu != NULL)
        cudaMemcpy(dd2_cpu, dd2_gpu, size, cudaMemcpyDeviceToHost);

    // Allocate the host input vectors
    uint array_size = number_of_points * sizeof(T);
    T *energy_cpu = (T *)malloc(array_size);
    T *factor_cpu = (T *)malloc(array_size);
    T *partial_densities_cpu = (T*)malloc(array_size);
    T *point_weights_cpu = (T*)malloc(array_size);
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
    printf("dxyz:"); print_vec_type<T>(dxyz_cpu, number_of_points);
    printf("dd1:"); print_vec_type<T>(dd1_cpu, number_of_points);
    printf("dd2:"); print_vec_type<T>(dd2_cpu, number_of_points);
    printf("energy:"); print_array<T>(energy_cpu, number_of_points);
    printf("factor:"); print_array<T>(factor_cpu, number_of_points);
    printf("point_weights:"); print_array<T>(point_weights_cpu, number_of_points);
    printf("partial_densities:"); print_array<T>(partial_densities_cpu, number_of_points);
    printf("=========\n");

    // Free memory.
    free(energy_cpu);
    free(factor_cpu);
    free(partial_densities_cpu);
    free(dd1_cpu);
    free(dd2_cpu);
    free(dxyz_cpu);
#endif
}

template <class T> 
void print_proxy_input (
    T* dens,
    const int number_of_points,
    const T* contracted_grad,
    const G2G::vec_type<T, 4>* grad,
    const G2G::vec_type<T, 4>* hess1,
    const G2G::vec_type<T, 4>* hess2)
{
#ifdef __CUDACC__
    // Bajar la info a CPU y luego imprimirla.
    // Copy to host the matrix data in gpu memory and
    // call the new libxcProxy.
    T* dens_cpu;
    T* contracted_grad_cpu;
    G2G::vec_type<T,4>* grad_cpu;
    G2G::vec_type<T,4>* hess1_cpu;
    G2G::vec_type<T,4>* hess2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = number_of_points * sizeof(G2G::vec_type<T,4>);
    uint array_size = number_of_points * sizeof(T);
    dens_cpu = (T*)malloc(size);
    contracted_grad_cpu = (T*)malloc(size);
    grad_cpu = (G2G::vec_type<T,4> *)malloc(size);
    hess1_cpu = (G2G::vec_type<T,4> *)malloc(size);
    hess2_cpu = (G2G::vec_type<T,4> *)malloc(size);

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
#endif
}

template <class T> 
void print_libxc_exchange_correlation_gpu_input (
    T* energy_gpu,
    T* factor_gpu,
    int points,
    T* accumulated_density_gpu,
    G2G::vec_type<T, 4>* dxyz_gpu,
    G2G::vec_type<T, 4>* dd1_gpu,
    G2G::vec_type<T, 4>* dd2_gpu)
{
#ifdef __CUDACC__
    // Bajar la info a CPU y luego imprimirla.
    // Copy to host the matrix data in gpu memory and
    // call the new libxcProxy.
    T* energy_cpu;
    T* factor_cpu;
    T* accumulated_density_cpu;
    G2G::vec_type<T,4>* dxyz_cpu;
    G2G::vec_type<T,4>* dd1_cpu;
    G2G::vec_type<T,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint size = points * sizeof(G2G::vec_type<T,4>);
    uint array_size = points * sizeof(T);
    energy_cpu = (T*)malloc(size);
    factor_cpu = (T*)malloc(size);
    accumulated_density_cpu = (T*)malloc(size);
    dxyz_cpu = (G2G::vec_type<T,4> *)malloc(size);
    dd1_cpu = (G2G::vec_type<T,4> *)malloc(size);
    dd2_cpu = (G2G::vec_type<T,4> *)malloc(size);

    // Copy data from device to host.
    if (energy_gpu != NULL)
        cudaMemcpy(energy_cpu, energy_gpu, size, cudaMemcpyDeviceToHost);
    if (factor_gpu != NULL)
        cudaMemcpy(factor_cpu, factor_gpu, size, cudaMemcpyDeviceToHost);
    if (dxyz_gpu != NULL)
	cudaMemcpy(dxyz_cpu, dxyz_gpu, size, cudaMemcpyDeviceToHost);
    if (dd1_gpu != NULL)
	cudaMemcpy(dd1_cpu, dd1_gpu, size, cudaMemcpyDeviceToHost);
    if (dd2_gpu != NULL)
        cudaMemcpy(dd2_cpu, dd2_gpu, size, cudaMemcpyDeviceToHost);
    if (accumulated_density_gpu != NULL)
	cudaMemcpy(accumulated_density_cpu, accumulated_density_gpu, array_size, cudaMemcpyDeviceToHost);

    printf("=============================================\n");
    printf("= input for libxc_exchange_correlation_gpu  =\n");
    printf("=============================================\n");
    printf("points:%i\n", points);
    printf("dxyz:"); print_vec_type(dxyz_cpu, points);
    printf("dd1:"); print_vec_type(dd1_cpu, points);
    printf("dd2:"); print_vec_type(dd2_cpu, points);
    printf("energy:"); print_array(energy_cpu, points);
    printf("factor:"); print_array(factor_cpu, points);
    printf("accumulated_density:"); print_array(accumulated_density_cpu, points);
    printf("==============================\n");

    // Free memory.
    free(energy_cpu);
    free(factor_cpu);
    free(accumulated_density_cpu);
    free(dd2_cpu);
    free(dd1_cpu);
    free(dxyz_cpu);
#endif
}

#endif