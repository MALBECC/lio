/*-*- mode: c -*-*/

#include <cstdio>
#include "cuda_extra.h"
#include "../matrix.h"
using namespace std;
using namespace G2G;

#define BLOCK_SIZE 32

__global__ void dists_kernel(float* coords, float* dists, int atoms);
__global__ void dists_cached_kernel(float* coords, float* dists, int atoms);
__global__ void dists_cached2_kernel(float* coords, float* dists, int atoms);

void calc_dists(const HostMatrixFloat& coords, HostMatrixFloat& dists) {
	CudaMatrixFloat gpu_coords(coords);
	CudaMatrixFloat gpu_dists(coords.width, coords.width);
	
	dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE);
	dim3 threads(coords.width, coords.width);
	dim3 gridSize = divUp(threads, blockSize);
	
	dists_kernel<<<gridSize, blockSize>>>(gpu_coords.data, gpu_dists.data, gpu_coords.width);
	
	dists = gpu_dists;
}

__global__ void dists_cached_kernel(float* coords, float* dists, int atoms)
{
	const uint3 pos = index(blockDim, blockIdx, threadIdx);
	

	__shared__ float coords_block[BLOCK_SIZE * BLOCK_SIZE];
	unsigned int shared_x = threadIdx.x * 2;
	unsigned int shared_y = threadIdx.y * 2 + 1;
	//unsigned int shared_x = threadIdx.x;
	//unsigned int shared_y = threadIdx.y + BLOCK_SIZE;

	coords_block[shared_x] = coords[pos.x];
	coords_block[shared_y] = coords[pos.y];
	
	__syncthreads();
	
	dists[pos.x + atoms * pos.y] = coords_block[shared_x] - coords_block[shared_y];

	__syncthreads();
}

#define NUM_BANKS 16

/**
 * Esta tiene sentido probarla solo para BLOCK_SIZE mayor a 16
 */
__global__ void dists_cached2_kernel(float* coords, float* dists, int atoms)
{
	const uint3 pos = index(blockDim, blockIdx, threadIdx);

	__syncthreads();
	_EMU(printf("t: %i, b: %i\n", threadIdx.x, blockIdx.x));
	
	__shared__ float coords_block[2 * (BLOCK_SIZE + BLOCK_SIZE / NUM_BANKS)];
	unsigned int shared_x = threadIdx.x * 2;
	unsigned int shared_y = threadIdx.y * 2 + 1;
	unsigned int padding_x = shared_x / NUM_BANKS;
	unsigned int padding_y = shared_y / NUM_BANKS;

	coords_block[shared_x + padding_x] = coords[pos.x];
	coords_block[shared_y + padding_y] = coords[pos.y];
	
	__syncthreads();
	
	dists[pos.x + atoms * pos.y] = coords_block[shared_x + padding_x] - coords_block[shared_y + padding_y];

	__syncthreads();
}

__global__ void dists_kernel(float* coords, float* dists, int atoms)
{
	const uint3 pos = index(blockDim, blockIdx, threadIdx);
	_EMU(printf("Hola %i\n", pos.x));
	dists[pos.x + atoms * pos.y] = coords[pos.x] - coords[pos.y];
}
