#ifndef __GLOBAL_MEMORY_POOL
#define __GLOBAL_MEMORY_POOL

#include <cassert>
#include <cstdlib>
#include "cuda/cuda_extra.h"

class globalMemoryPool
{
    public:
        static void init (double free_global_memory = 0.0);
        static bool tryAlloc (size_t size);
        static void dealloc (size_t size);
        static size_t getFreeMemory()
        {
            if(!_init) init();
            return _freeGlobalMemory;
        }

    private:
        static size_t _totalGlobalMemory;
        static size_t _freeGlobalMemory;
        static float _freeFactor;
        static bool  _init;
};
#endif
