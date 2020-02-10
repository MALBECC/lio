#ifndef _G2G_GPU_VARIABLES_H
#define _G2G_GPU_VARIABLES_H

extern __device__ __constant__ uint gpu_atoms;
#if FULL_DOUBLE
extern __device__ __constant__ double3 gpu_atom_positions[MAX_ATOMS];
extern __device__ __constant__ double gpu_normalization_factor;
#else
extern __device__ __constant__ float3 gpu_atom_positions[MAX_ATOMS];
extern __device__ __constant__ float gpu_normalization_factor;
#endif

#endif
