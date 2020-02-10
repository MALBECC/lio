#ifndef __FAKECOMPILE_H__
#define __FAKECOMPILE_H__

#if !GPU_KERNELS
#define __host__
#define __device__
#else
#include <cuda_runtime.h>
#endif

#ifdef __CUDACC__
    // Nothing
#else
    // Fake definitions for the kernels so they'll
    // compile under gcc or cc
    struct fake_blockDim {
	int x;
    };

    struct fake_blockIdx {
	int x;
    };

    struct fake_threadIdx {
	int x;
    };

    extern fake_blockDim blockDim;
    extern fake_blockIdx blockIdx;
    extern fake_threadIdx threadIdx;
#endif

#endif