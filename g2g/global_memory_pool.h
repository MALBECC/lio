#ifndef __GLOBAL_MEMORY_POOL
#define __GLOBAL_MEMORY_POOL

#include <cassert>
#include <cstdlib>
#include <vector>
#include "cuda/cuda_extra.h"

class GlobalMemoryPool {
 public:
  static void init(double free_global_memory = 0.0);
  static int tryAlloc(size_t size);
  static void dealloc(size_t size, int device);
  static void dealloc(size_t size);
  static size_t getFreeMemory();

 private:
  static std::vector<size_t> _totalGlobalMemory;
  static std::vector<size_t> _freeGlobalMemory;
  static float _freeFactor;
  static bool _init;
};
#endif
