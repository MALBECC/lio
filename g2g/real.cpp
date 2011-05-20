#include "common.h"
#include "init.h"



real3 to_real3(double3 f3) {
  return make_real3(f3.x, f3.y, f3.z);
}

creal3 to_creal3(double3 d3) {
  #if FULL_DOUBLE
  return d3;
  #else
  return creal3(d3);
  #endif
}

