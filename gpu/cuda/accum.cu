/*-*- mode: c -*-*/

#include <cstdio>
#include "cuda_extra.h"
#include "../matrix.h"
using namespace G2G;
using namespace std;

// Potencia de 2 / 2 * BLOCK_SIZE * k = input.width * input.height
#define BLOCK_SIZE 16

__global__ void accum_kernel(const float* input, float* output);

extern "C" void calc_accum(const Matrix& input, Matrix& output)
{
	float* gpu_input;
	float* gpu_output;
	
	cudaMalloc((void**)&gpu_input, input.bytes());
	cudaMalloc((void**)&gpu_output, output.bytes());
	
	cudaMemcpy(gpu_input, input.data, input.bytes(), cudaMemcpyHostToDevice);	

	dim3 dimBlock(BLOCK_SIZE);
	dim3 dimSize((input.width * input.height) / (BLOCK_SIZE * 2));
	accum_kernel<<<dimSize, dimBlock>>>(gpu_input, gpu_output);
	
	cudaMemcpy(output.data, gpu_output, output.bytes(), cudaMemcpyDeviceToHost);
	
	cudaThreadSynchronize();
		
	cudaFree(gpu_input);
	cudaFree(gpu_output);	
}

/** TODO: probar sin cargar a shared **/
__global__ void accum_kernel(const float* input, float* output)
{
	
	__shared__ float shared_data[BLOCK_SIZE * 2];
	const uint pos = index(blockDim, blockIdx, threadIdx).x * 2;
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
	
	if (threadIdx.x == 0) *output = shared_data[shared_pos];
}

