#ifndef _GPU_VARIABLES_H
#define _GPU_VARIABLES_H

#include "../common.h"

//#define GAMMA_INC 1
//#define GAMMA_LENGTH 100

//cudaArray* gammaArray;

/* grid */
__device__ __constant__ uint gpu_atoms;
__device__ __constant__ uint gpu_clatoms;
__device__ __constant__ uint gpu_Iexch;
__device__ __constant__ uint gpu_m;
#if FULL_DOUBLE
__device__ __constant__ double3 gpu_atom_positions[MAX_ATOMS];
__device__ __constant__ double3 gpu_clatom_positions[MAX_CLATOMS];
__device__ __constant__ double gpu_clatom_charges[MAX_CLATOMS];
__device__ __constant__ double gpu_normalization_factor;

__device__ double* gpu_str;
__device__ __constant__ double gpu_fac[17];
#else
__device__ __constant__ float3 gpu_atom_positions[MAX_ATOMS];
__device__ __constant__ float3 gpu_clatom_positions[MAX_CLATOMS];
__device__ __constant__ float gpu_clatom_charges[MAX_CLATOMS];
__device__ __constant__ float gpu_normalization_factor;

__device__ float* gpu_str;
__device__ __constant__ float gpu_fac[17];
#endif

#endif
