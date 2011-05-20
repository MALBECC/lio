#ifndef __REAL_H__
#define __REAL_H__

#include <float.h>
#include "cuda/cuda_extra.h"

#if FULL_DOUBLE
  #define real double
  #define creal3 double3
  #define make_creal3 make_double3
  #define real3 double3
  #define make_real3 make_double3
  #define exp exp
  #define pow pow
  #define log log
  #define REAL_MIN DBL_MIN
  #define gsl_matrix_real_view gsl_matrix_view
  #define gsl_matrix_real_view_array gsl_matrix_view_array
  #define gsl_vector_real_const_view gsl_vector_const_view
  #define gsl_vector_real_const_view_array gsl_vector_const_view_array
#else
  #define real float
  #define creal3 cfloat3
  #define make_creal3 cfloat3
  #define real3 float3
  #define make_real3 make_float3
  #define exp expf
  #define pow powf
  #define log logf
  #define REAL_MIN FLT_MIN
  #define gsl_matrix_real_view gsl_matrix_float_view
  #define gsl_matrix_real_view_array gsl_matrix_float_view_array
  #define gsl_vector_real_const_view gsl_vector_float_const_view
  #define gsl_vector_real_const_view_array gsl_vector_float_const_view_array
#endif

real3 to_real3(double3 f3);
creal3 to_creal3(double3 d3);

#endif
