#ifndef _COULOMB_GPU_VARIABLES_H
#define _COULOMB_GPU_VARIABLES_H

__device__ __constant__ uint gpu_out_offsets[NUM_TERM_TYPES];
__device__ __constant__ uint gpu_s_end;
__device__ __constant__ uint gpu_p_offset;
__device__ __constant__ uint gpu_p_end;
__device__ __constant__ uint gpu_d_offset;
__device__ __constant__ uint gpu_d_end;

#endif
