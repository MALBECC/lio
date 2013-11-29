#ifndef __GLOBAL_MEMORY_POOL
#define __GLOBAL_MEMORY_POOL
#include <cassert>
#include <cstdlib>
#include "cuda/cuda_extra.h"

class globalMemoryPool
{
    public:
        static void init(float free_factor=0.3)
        {
            //Aca tenemos que leer el GPU_global y restar un factor de tolerancia (1/4?)
#if !CPU_KERNELS
            size_t free_memory, total_memory;

	        cudaGetMemoryInfo(free_memory, total_memory);

            _totalGlobalMemory=total_memory;
            if(free_factor>1.0f) free_factor=1.0f;
            if(free_factor<0.0f) free_factor=0.0f;
            _freeFactor=free_factor;

            _freeGlobalMemory=static_cast<size_t>(static_cast<float>(free_memory)*_freeFactor);
#else
            _totalGlobalMemory=0;
            _freeGlobalMemory=0;
            _freeFactor=0.0f;
#endif
        }
        static bool tryAlloc(size_t size);
        static void dealloc(size_t size);
        static size_t getFreeMemory()
        {
            return _freeGlobalMemory;
        }
    private:
        static size_t _totalGlobalMemory;
        static size_t _freeGlobalMemory;
        static float _freeFactor;
};
#endif
