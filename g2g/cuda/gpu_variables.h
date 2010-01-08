#ifndef _GPU_VARIABLES_H
#define _GPU_VARIABLES_H

#include "../init.h"

/* grid */
__device__ __constant__ float3 gpu_atom_positions[MAX_ATOMS];
__device__ __constant__ float gpu_rm[MAX_ATOMS];  // TODO: merge this into a float4 with positions
__device__ __constant__ uint gpu_atoms;
__device__ __constant__ uint gpu_Iexch = 0;
__device__ __constant__ uint gpu_nco = 0;
__device__ __constant__ float gpu_normalization_factor = 1.0f;

#endif
