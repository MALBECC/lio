#if !GPU_KERNELS
#   define __host__
#   define __device__
#else
#   include <cuda_runtime.h>
#endif

#if !INTELCOMP
#endif

