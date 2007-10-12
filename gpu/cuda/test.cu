/*-*- mode: c -*-*/

#include "cuda_extra.h"

#define BLOCK_SIZE 16

__global__ void dists_kernel(float* coords, float* dists, int atoms);
__global__ void dists_cached_kernel(float* coords, float* dists, int atoms);

extern "C" void run_test(float* coords_cpu, float* dists_cpu, int atoms) {
	float* coords;
	float* dists;
	
	cudaMalloc((void**)&coords, atoms * sizeof(float));
	cudaMemcpy(coords, coords_cpu, atoms * sizeof(float), cudaMemcpyHostToDevice);
	
	//int pitch;	
	//cudaMallocPitch(&dists, &pitch, atoms * atoms * sizeof(float), atoms);
	cudaMalloc((void**)&dists, atoms * atoms * sizeof(float));

	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimSize(atoms / dimBlock.x, atoms / dimBlock.y);
	//dists_kernel<<<dimSize, dimBlock>>>(coords, dists, atoms);
	dists_cached_kernel<<<dimSize, dimBlock>>>(coords, dists, atoms);	
	
	cudaMemcpy(dists_cpu, dists, atoms * atoms * sizeof(float), cudaMemcpyDeviceToHost);	
		
	cudaFree(coords);
	cudaFree(dists);	
}

__global__ void dists_cached_kernel(float* coords, float* dists, int atoms)
{
	const uint3 pos = index(blockDim, blockIdx, threadIdx);

	__shared__ float coords_block[BLOCK_SIZE * BLOCK_SIZE];

	coords_block[threadIdx.x] = coords[pos.x];
	coords_block[threadIdx.y + BLOCK_SIZE] = coords[pos.y];
	
	__syncthreads();
	
	dists[pos.x + atoms * pos.y] = coords_block[threadIdx.x] - coords_block[threadIdx.y + BLOCK_SIZE];

	__syncthreads();
}

__global__ void dists_kernel(float* coords, float* dists, int atoms)
{
	const uint3 pos = index(blockDim, blockIdx, threadIdx);
	dists[pos.x + atoms * pos.y] = coords[pos.x] - coords[pos.y];
}
