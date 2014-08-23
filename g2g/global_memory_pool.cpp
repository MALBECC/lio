#include "global_memory_pool.h"
#include <cassert>

// Only include it to cudaGetMemoryInfo to accurately measure the device memory.
#include "cuda_includes.h"

//TryAlloc return 1 if error, 0 if success
bool globalMemoryPool::tryAlloc(size_t size)
{
    if (!_init) init();
    if (_freeGlobalMemory < size)
        return true;
    _freeGlobalMemory -= size;
    return false;
}

void globalMemoryPool::dealloc(size_t size)
{
    if(!_init) init();
    assert (_freeGlobalMemory>size);
    _freeGlobalMemory += size;
}

void globalMemoryPool::init(double free_global_memory)
{
    //Aca tenemos que leer el GPU_global y restar un factor de tolerancia (1/4?)
#if !CPU_KERNELS
    double free_factor = free_global_memory;
    size_t free_memory, total_memory;

    cudaGetMemoryInfo(free_memory, total_memory);

    _totalGlobalMemory=total_memory;
    if (free_factor > 1.0f) free_factor = 1.0f;
    if (free_factor < 0.0f) free_factor = 0.0f;
    _freeFactor=free_factor;

    _freeGlobalMemory=static_cast<size_t>(static_cast<double>(free_memory)*_freeFactor);
#else
    _totalGlobalMemory = 0;
    _freeGlobalMemory = 0;
    _freeFactor = 0.0f;
#endif
    _init = true;
}

size_t globalMemoryPool::_freeGlobalMemory = 0;
size_t globalMemoryPool::_totalGlobalMemory = 0;
float  globalMemoryPool::_freeFactor = 0.8;
bool   globalMemoryPool::_init = false;
