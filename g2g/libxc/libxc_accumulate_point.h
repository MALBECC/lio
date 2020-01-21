#define WIDTH 4

#include "print_utils.h"
#include "libxcproxy_cuda.h"
#include "../timer.h"

extern "C" void g2g_timer_sum_start_(const char* timer_name, unsigned int length_arg);
extern "C" void g2g_timer_sum_stop_(const char* timer_name, unsigned int length_arg);
extern "C" void g2g_timer_sum_pause_(const char* timer_name, unsigned int length_arg);

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// gpu_accumulate_point_for_libx
//
// point_weights:
// points: the size of the arrays.
// block_height:
// partial_density_in:
// dxyz_in:
// dd1_in:
// dd2_in:
// accumulated_density:
// dxyz_accum:
// dd1_accum:
// dd2_accum:
// 
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// gpu_accumulate_energy_and_forces_from_libxc
// 
// energy: array of energy points,
// factor: 
// point_weights:
// points: the size of the arrays.
// accumulated_density:
//
//
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// libxc_exchange_correlation_cpu
// Use the libxc cpu version to compute the exchange-correlation values 
//
// energy_gpu: array of energy points,
// factor_gpu: array of factors points.
// points: the size of all the input arrays.
// accumulated_density_gpu: the contracted grad for libxc
// dxyz_gpu:
// dd1_gpu:
// dd2_gpu:
//
// Note: all the pointer data are pointers in CUDA memory.
//

template<class T, bool compute_energy, bool compute_factor, bool lda>
    void libxc_exchange_correlation_cpu (LibxcProxy_cuda<T, WIDTH>* libxcProxy,
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
    // Parameters
    T* exc;
    T* corr;
    T* y2a;

    // Now alloc memory for the data
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

<<<<<<< HEAD
=======
template<class T>
__global__ void all_ContractGradients(G2G::vec_type<T,WIDTH>* dxyz, G2G::vec_type<T,WIDTH>* tredxyz,
                                      T* sigma, T* cruz, uint points)
{
    uint i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < points) {

       sigma[i] = ( dxyz[i].x*dxyz[i].x ) + ( dxyz[i].y*dxyz[i].y ) + ( dxyz[i].z*dxyz[i].z ); 
       cruz[i]  = 0.5f * ( (tredxyz[i].x*dxyz[i].x)+(tredxyz[i].y*dxyz[i].y)+(tredxyz[i].z*dxyz[i].z) );

    }

}

template<class T>
__global__ void group_coefficients(G2G::vec_type<T,WIDTH>* dxyz_in, G2G::vec_type<T,WIDTH>* txyz_in,
                                   G2G::vec_type<T,WIDTH>* Dout, G2G::vec_type<T,WIDTH>* Tout,
                                   T* C_out,T* C_local, uint points)
{

    uint idx = blockDim.x * blockIdx.x + threadIdx.x;

    if ( idx < points ) {

      C_out[idx] = C_local[idx] * 0.5f;

      Dout[idx] = G2G::vec_type<T,WIDTH>(C_local[points+idx]*dxyz_in[idx].x*0.5f,
                                         C_local[points+idx]*dxyz_in[idx].y*0.5f,
                                         C_local[points+idx]*dxyz_in[idx].z*0.5f,0.0f);

      Tout[idx] = G2G::vec_type<T,WIDTH>(C_local[2*points+idx]*txyz_in[idx].x,
                                         C_local[2*points+idx]*txyz_in[idx].y,
                                         C_local[2*points+idx]*txyz_in[idx].z,0.0f);


   } // end if valid thread

}


// Call from gpu to libxc on gpu in LR
template<class T, bool compute_energy, bool compute_factor, bool lda>
void libxc_gpu_coefLR(LibxcProxy_cuda<T, WIDTH>* libxcProxy,
                      uint points, T* dens_gpu, T* tred_gpu,
                      G2G::vec_type<T,WIDTH>* dxyz_gpu, G2G::vec_type<T,WIDTH>* tredxyz_gpu,
                      // outputs
                      G2G::vec_type<T,WIDTH>* Dxyz, G2G::vec_type<T,WIDTH>* Txyz,
                      T* coef_gpu)
{
   cudaError_t err = cudaSuccess;

// CALCULATE GRADIENTS
   T* sigma = NULL;
   T* cruz  = NULL;
   T* coef_local  = NULL;
   int size = sizeof(T) * points;

   err = cudaMalloc((void**)&sigma, size);
   if ( err != cudaSuccess ) cout << "Falied to allocate sigma on GPU" << endl;

   err = cudaMalloc((void**)&cruz, size);
   if ( err != cudaSuccess ) cout << "Falied to allocate cruz on GPU" << endl;

   err = cudaMalloc((void**)&coef_local, 3*size);
   if ( err != cudaSuccess ) cout << "Falied to allocate coef_local on GPU" << endl;

   int threadsPerBlock = 256;
   int blocksPerGrid = (points + threadsPerBlock - 1) / threadsPerBlock;
   all_ContractGradients<<<blocksPerGrid,threadsPerBlock>>>(dxyz_gpu,tredxyz_gpu,
                           sigma,cruz,points);

// LIBXC CALCULATE ALL POINTS
   libxcProxy->coefLR(points,dens_gpu,sigma,tred_gpu,cruz,coef_local);
   group_coefficients<<<blocksPerGrid,threadsPerBlock>>>(dxyz_gpu,tredxyz_gpu,
                      Dxyz, Txyz, coef_gpu, coef_local, points);

// Free memory
   cudaFree(sigma);      sigma = NULL;
   cudaFree(cruz);       cruz = NULL;
   cudaFree(coef_local); coef_local = NULL;
}
<<<<<<< HEAD
>>>>>>> ce4ae22... save GS rho values on GPU
=======

// Call from gpu to libxc on gpu in Zvector
template<class T, bool compute_energy, bool compute_factor, bool lda>
void libxc_gpu_coefZv(LibxcProxy_cuda<T, WIDTH>* libxcProxy,
                      uint points, T* dens_gpu, T* tred_gpu,
                      G2G::vec_type<T,WIDTH>* dxyz_gpu, G2G::vec_type<T,WIDTH>* tredxyz_gpu,
                      // outputs
                      G2G::vec_type<T,WIDTH>* Dxyz, G2G::vec_type<T,WIDTH>* Txyz,
                      T* coef_gpu)
{
   cudaError_t err = cudaSuccess;

// CALCULATE CONTRACTED GRADIENTS
   T* sigma = NULL;
   T* cruz  = NULL;
   T* coef_local  = NULL;
   int size = sizeof(T) * points;

   err = cudaMalloc((void**)&sigma, size);
   if ( err != cudaSuccess ) cout << "Falied to allocate sigma on GPU" << endl;

   err = cudaMalloc((void**)&cruz, size);
   if ( err != cudaSuccess ) cout << "Falied to allocate cruz on GPU" << endl;

   err = cudaMalloc((void**)&coef_local, 3*size);
   if ( err != cudaSuccess ) cout << "Falied to allocate coef_local on GPU" << endl;

   int threadsPerBlock = 256;
   int blocksPerGrid = (points + threadsPerBlock - 1) / threadsPerBlock;
   all_ContractGradients<<<blocksPerGrid,threadsPerBlock>>>(dxyz_gpu,tredxyz_gpu,
                           sigma,cruz,points);

   cudaFree(cruz); cruz = NULL; // In this case, cruz is not referenced
   libxcProxy->coefZv(points,dens_gpu,sigma,tred_gpu,dxyz_gpu,tredxyz_gpu,coef_local);
   group_coefficients<<<blocksPerGrid,threadsPerBlock>>>(dxyz_gpu,tredxyz_gpu,
                      Dxyz, Txyz, coef_gpu, coef_local, points);

   cudaFree(sigma); sigma = NULL;
   cudaFree(coef_local); coef_local = NULL;
}
>>>>>>> 8a77770... Relaxed Density Matrix works on cuda
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// calculateContractedGradient
// Calculate a contracted gradient for a given gradient.
// grad: array of vec_type.
// contracted_grad: the contracted gradient as a result of 
//                  contract the variable grad.
// numElements: the size of the arrays.
//
// Note: all the pointer data are pointers in CUDA memory.
//

template<class T>
__global__ void calculateContractedGradient(G2G::vec_type<T,WIDTH>* grad, T* contracted_grad, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        contracted_grad[i] = (grad[i].x * grad[i].x) + (grad[i].y * grad[i].y) + (grad[i].z * grad[i].z);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// vectorAdd
// Adds two arrays and puts the result in a third array.
// A simple array addition.
// T* A: array of T
// T* B: array of T
// T* C: the result of adding the arrays A+B.
// numElements: the size of the arrays.
//
// Note: all the pointer data are pointers in CUDA memory.
//
template<class T>
__global__ void vectorAdd(const T* A, const T* B, T* C, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        C[i] = A[i] + B[i];
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// libxc_exchange_correlation_gpu
// Use the libxc gpu version to compute the exchange-correlation values.
// energy_gpu: array of energy points,
// factor_gpu: array of factors points.
// points: the size of all the input arrays.
// accumulated_density_gpu: the contracted grad for libxc
// dxyz_gpu:
// dd1_gpu:
// dd2_gpu:
//
// Note: all the pointer data are pointers in CUDA memory.
//
template<class T, bool compute_energy, bool compute_factor, bool lda>
    void libxc_exchange_correlation_gpu (LibxcProxy_cuda<T, WIDTH>* libxcProxy,
	T* energy_gpu,
	T* factor_gpu,
	uint points,
	T* accumulated_density_gpu,
	G2G::vec_type<T,WIDTH>* dxyz_gpu,
        G2G::vec_type<T,WIDTH>* dd1_gpu,
	G2G::vec_type<T,WIDTH>* dd2_gpu)
{
    //printf("libxc_exchage_correlation_gpu<%u,%d,%d,%d>(%u) \n", sizeof(T), compute_energy, compute_factor, lda, points);
    //print_libxc_exchange_correlation_gpu_input (energy_gpu, factor_gpu, points, accumulated_density_gpu, dxyz_gpu, dd1_gpu, dd2_gpu);

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
    //printf ("array_size: %u \n", array_size);

    // Alloc memory for cuda variables.
    err = cudaMalloc((void **)&exc_gpu, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device exc_gpu! %u \n", err);
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

    //err = cudaMalloc((void**)&contracted_gradient, sizeof(T)*points);
    err = cudaMalloc((void**)&contracted_gradient, array_size);
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
        //printf("Calling calculateContractedGradient() \n");
	// Call the Kernel to compute the contracted gradient
        //printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
	calculateContractedGradient<<<blocksPerGrid, threadsPerBlock>>>(dxyz_gpu, contracted_gradient, points);

	/////////////////////////////////////////
	// Compute exc_corr using LIBXC GPU.
        //printf("Calling libxcProxy->doGGA() \n");

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

