/*-*- mode: c -*-*/

#include <cstdio>
#include "cuda_extra.h"
#include "../matrix.h"
using namespace G2G;
using namespace std;

// Potencia de 2 / 2 * BLOCK_SIZE * k = input.width * input.height
#define BLOCK_SIZE 16

__global__ void accum_kernel(const float* input, float* output);

extern "C" void calc_accum(const HostMatrix& input, HostMatrix& output)
{	
	// por el momento asume que dimSize.x es potencia de 2
	CudaMatrix gpu_input(input);

	// define sizes
	dim3 dimBlock(BLOCK_SIZE);
	dim3 dimSize((input.width * input.height) / (BLOCK_SIZE * 2));

	CudaMatrix gpu_output(dimSize.x);
	
	// accumulate with blocks of BLOCK_SIZE, do this until there are no more block to process
	float* gpu_input_data = gpu_input.data;
	
	do {
		printf("Pasada\n");
		accum_kernel<<<dimSize, dimBlock>>>(gpu_input_data, gpu_output.data);

		/*float test;
		cudaMemcpy(&test, gpu_output.data, sizeof(float), cudaMemcpyDeviceToHost);
		printf("----- %f\n", test);*/
		
		dimSize.x /= (BLOCK_SIZE * 2);
		gpu_input_data = gpu_output.data;
	} while (dimSize.x >= 1);
	
	cudaThreadSynchronize();

	output.copy_submatrix(gpu_output, 1);
}

/** TODO: probar sin cargar a shared **/
__global__ void accum_kernel(const float* input, float* output)
{
	
	__shared__ float shared_data[BLOCK_SIZE * 2];
	const uint pos = (blockDim.x * blockIdx.x + threadIdx.x) * 2;
	const uint shared_pos = threadIdx.x * 2;
	
	shared_data[shared_pos] = input[pos];
	shared_data[shared_pos + 1] = input[pos + 1];	// this could be optional
		
	uint offset = 1;
	
	for (unsigned int d = BLOCK_SIZE; d > 0; d /= 2) {
		if (threadIdx.x % offset == 0)
			shared_data[shared_pos] += shared_data[shared_pos + offset];
		
		offset *= 2;

		__syncthreads();		
	}

	// printf("t: %i, b: %i\n", threadIdx.x, blockIdx.x);
	if (threadIdx.x == 0) output[blockIdx.x] = shared_data[shared_pos];
}

