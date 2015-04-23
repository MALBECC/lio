#ifndef _QMMM_COUL_GPU_VARIABLES_H
#define _QMMM_COUL_GPU_VARIABLES_H

#include "../common.h"

#define NUM_COULOMB_TERM_TYPES 6

cudaArray* gammaArray;

/* grid */
extern __device__ __constant__ uint gpu_atoms;
__device__ __constant__ uint gpu_clatoms;
__device__ __constant__ uint gpu_m;
__device__ __constant__ uint gpu_d_offset; // Needed for d-d in QM/MM
__device__ __constant__ uint gpu_dens_gauss;
__device__ __constant__ uint gpu_dens_s_gauss;
__device__ __constant__ uint gpu_dens_p_gauss;
__device__ __constant__ uint gpu_out_offsets[NUM_COULOMB_TERM_TYPES];
#if FULL_DOUBLE
extern __device__ __constant__ double3 gpu_atom_positions[MAX_ATOMS];
extern __device__ __constant__ double gpu_normalization_factor;
#else
extern __device__ __constant__ float3 gpu_atom_positions[MAX_ATOMS];
extern __device__ __constant__ float gpu_normalization_factor;
#endif

#if !QMMM_MP || !COULOMB_MP || FULL_DOUBLE
//__device__ double* gpu_str;
__device__ __constant__ double gpu_fac[17];
#else
//__device__ float* gpu_str;
__device__ __constant__ float gpu_fac[17];
#endif

#endif
