#ifndef _SCALAR_VECTOR_TYPES_GPU_H
#define	_SCALAR_VECTOR_TYPES_GPU_H

#include "datatypes/cpu_primitives.h"

namespace G2G {
  template<class T, unsigned int n> class vec_type_gpu { };

  template<> class vec_type_gpu<float, 2> : public float2 {
    public:
      __device__ __host__ vec_type_gpu(void) {}
      __device__ __host__ vec_type_gpu(float _x, float _y) { this->x = _x; this->y = _y; }
      typedef float2 base_type;
  };
  template<> class vec_type_gpu<float, 3> : public float3 {
    public:
      __device__ __host__ vec_type_gpu(void) {}
      __device__ __host__ vec_type_gpu(const float4& other) { float3::x = other.x; float3::y = other.y; float3::z = other.z; }
      __device__ __host__ vec_type_gpu(const float3& other) { float3::x = other.x; float3::y = other.y; float3::z = other.z; }
      __device__ __host__ explicit vec_type_gpu(const double3& other) { float3::x = other.x; float3::y = other.y; float3::z = other.z; }
      __device__ __host__ vec_type_gpu(float _x, float _y, float _z) { float3::x = _x; float3::y = _y; float3::z = _z; }
      __device__ __host__ float length2() const { return float3::x * float3::x + float3::y * float3::y + float3::z * float3::z; }
      typedef float3 base_type;

  };

  template<> class vec_type_gpu<float, 4> : public float4 {
    public:
      __device__ __host__ vec_type_gpu(void) {}
      __device__ __host__ vec_type_gpu(const float3& other) { float4::x = other.x; float4::y = other.y; float4::z = other.z; float4::w = 0.0f; }
      __device__ __host__ vec_type_gpu(const float4& other) { this->x = other.x; this->y = other.y; this->z = other.z; this->w = other.w; }
      __device__ __host__ vec_type_gpu(float _x, float _y, float _z, float _w) { this->x = _x; this->y = _y; this->z = _z; this->w = _w; }
      typedef float4 base_type;
  };

  template<> class vec_type_gpu<double, 2> : public double2 {
    public:
      __device__ __host__ vec_type_gpu(void) {}
      __device__ __host__ vec_type_gpu(double _x, double _y) { this->x = _x; this->y = _y; }
      typedef double2 base_type;
  };


  template<> class vec_type_gpu<double, 3> : public double3 {
    public:
      __device__ __host__ vec_type_gpu(void) {}
      __device__ __host__ vec_type_gpu(const double3& other) : double3(other) { }
      __device__ __host__ explicit vec_type_gpu(const float3& other) { double3::x = other.x; double3::y = other.y; double3::z = other.z; }
      __device__ __host__ explicit vec_type_gpu(const double4& other) { double3::x = other.x; double3::y = other.y; double3::z = other.z; }
      __device__ __host__ vec_type_gpu(double _x, double _y, double _z) { double3::x = _x; double3::y = _y; double3::z = _z; }
      typedef double3 base_type;
      __device__ __host__ double length2(void) const { return double3::x * double3::x + double3::y * double3::y + double3::z * double3::z; }
  };

  template<> class vec_type_gpu<double, 4> : public double4 {
    public:
      __device__ __host__ vec_type_gpu(void) {}
      __device__ __host__ explicit vec_type_gpu(const double3& other) { this->x = other.x; this->y = other.y; this->z = other.z; }
      __device__ __host__ explicit vec_type_gpu(const double4& other) { this->x = other.x; this->y = other.y; this->z = other.z; }
      __device__ __host__ vec_type_gpu(double _x, double _y, double _z, double _w) { this->x = _x; this->y = _y; this->z = _z; this->w = _w; }
  };

  template<> class vec_type_gpu<int2, 4>{
    public:
      __device__ __host__ vec_type_gpu(void) {}
      __device__ __host__ vec_type_gpu(int2 _x, int2 _y, int2 _z, int2 _w) { this->x = _x; this->y = _y; this->z = _z; this->w = _w; }
      int2 x, y, z, w;
  };
}

#endif	/* _SCALAR_VECTOR_TYPES_GPU_H */

