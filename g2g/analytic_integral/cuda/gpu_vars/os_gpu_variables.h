#ifndef _OS_GPU_VARIABLES_H
#define _OS_GPU_VARIABLES_H

extern cudaArray* gammaArray;

extern __device__ __constant__ uint gpu_m;

#if !AINT_MP || FULL_DOUBLE
extern __device__ __constant__ double gpu_fac[17];
extern texture<int2, cudaTextureType2D, cudaReadModeElementType>
    str_tex;  // Texture for STR array (used in F(m,U))
#else
extern __device__ __constant__ float gpu_fac[17];
extern texture<float, cudaTextureType2D, cudaReadModeElementType> str_tex;
#endif

extern __device__ __constant__ uint TERM_TYPE_GAUSSIANS[6];  // How many
                                                             // individual force
                                                             // terms each type
                                                             // (s-s,etc) is
                                                             // calculating
extern __device__ __constant__ uint gpu_atom_types[MAX_ATOMS];
extern __device__ __constant__ uint gpu_atom_Z[MAX_ATOMS];

#endif
