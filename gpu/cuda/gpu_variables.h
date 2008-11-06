#ifndef _GPU_VARIABLES_H
#define _GPU_VARIABLES_H

#include "../init.h"

/* grid */
__device__ __constant__ float3 gpu_atom_positions[MAX_ATOMS];
__device__ __constant__ uint gpu_atoms;

/* pot_kernel constants */
__device__ __constant__ float pot_alpha = -0.738558766382022447;
__device__ __constant__ float pot_gl = 0.620350490899400087;
	
__device__ __constant__ float pot_vosko_a1 = 0.03109205;
__device__ __constant__ float pot_vosko_b1 = 3.72744;
__device__ __constant__ float pot_vosko_c1 = 12.9352;
__device__ __constant__ float pot_vosko_x0 = -0.10498;
__device__ __constant__ float pot_vosko_q = 6.15199066246304849;
__device__ __constant__ float pot_vosko_a16 = 0.005182008333;
__device__ __constant__ float pot_vosko_a2 = 0.015546025;
__device__ __constant__ float pot_vosko_b2 = 7.06042;
__device__ __constant__ float pot_vosko_c2 = 18.0578;
__device__ __constant__ float pot_vosko_x02 = -0.32500;
__device__ __constant__ float pot_vosko_q2 = 4.7309269;
__device__ __constant__ float pot_vosko_a26 = 0.0025910042;

__device__ __constant__ uint gpu_Iexch = 0;
__device__ __constant__ uint gpu_nco = 0;
__device__ __constant__ float gpu_normalization_factor = 1.0f;


#endif
