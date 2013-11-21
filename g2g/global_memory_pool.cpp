#include "global_memory_pool.h"
#include <cassert>
//TryAlloc return 1 if error, 0 if success
bool globalMemoryPool::tryAlloc(size_t size)
{
            if(_freeGlobalMemory<size)
                return true; 
            _freeGlobalMemory-=size; 
            return false;
}

void globalMemoryPool::dealloc(size_t size)
{ 
    assert(_freeGlobalMemory>size); 
    _freeGlobalMemory+=size;
}
size_t globalMemoryPool::_freeGlobalMemory=0;
size_t globalMemoryPool::_totalGlobalMemory=0;
float  globalMemoryPool::_freeFactor=0.8;
