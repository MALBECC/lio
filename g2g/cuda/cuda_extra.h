#ifndef __CUTOOLS_H__
#define __CUTOOLS_H__

#if GPU_KERNELS
#include <cuda_runtime.h>
#include <cuda.h>
#endif

#define WARP_SIZE 32
#define BANKS 32

#include <stdexcept>
#include <iostream>
#include <string>

#include <cmath>

#if CPU_KERNELS && !GPU_KERNELS
#define __device__
#define __host__
#include "datatypes/cpu_primitives.h"
#endif

// TODO: usar cutil sdk para todos estos operadores, o usar classes de C++
// templatizadas
/** operators **/
inline __device__ __host__ float2 operator-(const float2 a) {
  return make_float2(-a.x, -a.y);
}

inline __device__ __host__ double3 operator-(const double3& a,
                                             const double3& b) {
  return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __device__ __host__ double3 operator-(const double3& a, double b) {
  return make_double3(a.x - b, a.y - b, a.z - b);
}

inline __device__ __host__ double3 operator+(const double3& a,
                                             const double3& b) {
  return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __device__ __host__ double4 operator+(const double4& a,
                                             const double4& b) {
  return make_double4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

inline __device__ __host__ double4 operator*(const double4& a,
                                             const double4& b) {
  return make_double4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

inline __device__ __host__ float4 operator+(const float4& a, const float4& b) {
  return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

inline __device__ __host__ float4 operator*(const float4& a, const float4& b) {
  return make_float4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

inline __device__ __host__ double3 operator/(const double3& a, double b) {
  return make_double3(a.x / b, a.y / b, a.z / b);
}

inline __device__ __host__ double3 operator/(const double3& a, const uint& b) {
  return make_double3(a.x / b, a.y / b, a.z / b);
}

inline __device__ __host__ uint3 ceil_uint3(const double3& a) {
  return make_uint3(static_cast<uint>(ceil(a.x)), static_cast<uint>(ceil(a.y)),
                    static_cast<uint>(ceil(a.z)));
}

inline __device__ __host__ uint3 floor_uint3(const double3& a) {
  return make_uint3(static_cast<uint>(floor(a.x)),
                    static_cast<uint>(floor(a.y)),
                    static_cast<uint>(floor(a.z)));
}

inline __device__ __host__ float3 operator*(const float3& a, const float3& b) {
  return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __device__ __host__ double3 operator*(const double3& a, double b) {
  return make_double3(a.x * b, a.y * b, a.z * b);
}

inline __device__ __host__ float3 operator/(const float3& a, const float3& b) {
  return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}

inline __device__ __host__ float3 operator+(const float3& a, const float3& b) {
  return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __device__ __host__ float3 operator+(const float3& a, const float& b) {
  return make_float3(a.x + b, a.y + b, a.z + b);
}

inline __device__ __host__ float3 operator-(const float3 a, const float3 b) {
  return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __device__ __host__ float3 operator-(const float3 a, const float b) {
  return make_float3(a.x - b, a.y - b, a.z - b);
}

inline __device__ __host__ float3 operator+(const uint3 a, const float3 b) {
  return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __device__ __host__ uint3 operator*(const uint3 a, const uint3 b) {
  return make_uint3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __device__ __host__ uint3 operator*(const dim3 a, const uint3 b) {
  return make_uint3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __device__ __host__ uint3 operator+(const uint3 a, const uint3 b) {
  return make_uint3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __device__ __host__ bool operator<(const uint3 i, const uint3 r) {
  return i.x < r.x && i.y < r.y && i.z < r.z;
}

inline __device__ __host__ float3 operator*(const float3& a, float b) {
  return make_float3(a.x * b, a.y * b, a.z * b);
}

inline __device__ __host__ double3 operator*(const double3& a,
                                             const double3& b) {
  return make_double3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __device__ __host__ float3 operator*(float b, const float3& a) {
  return make_float3(a.x * b, a.y * b, a.z * b);
}

inline __device__ __host__ float4 operator*(const float4& a, float b) {
  return make_float4(a.x * b, a.y * b, a.z * b, a.w * b);
}

inline __device__ __host__ float1 operator*(const float1& a, float b) {
  return make_float1(a.x * b);
}

inline __device__ __host__ uint3 index(const dim3 bd, const uint3 bi,
                                       const uint3 ti) {
  return bd * bi + ti;
}

inline __device__ __host__ uint index_x(const dim3 bd, const uint3 bi,
                                        const uint3 ti) {
  return bd.x * bi.x + ti.x;
}

inline __device__ __host__ uint reference(const uint3 i, const uint3 g) {
  return i.x + i.y * g.x + i.z * g.x * g.y;
}

inline __device__ __host__ uint len2(const float3 a) {
  return (uint)a.x * (uint)a.x + (uint)a.y * (uint)a.y + (uint)a.z * (uint)a.z;
}

inline __device__ __host__ void swap(float4& a, float4& b) {
  float4 tmp = a;
  a = b;
  b = tmp;
}

inline __device__ __host__ void negate(float2& a) {
  a.x = -a.x;
  a.y = -a.y;
}

inline __device__ __host__ const float4& min(const float4& a, const float4& b) {
  return (a.x < b.x) ? a : b;
}

inline __device__ __host__ const float4& max(const float4& a, const float4& b) {
  return (a.x > b.x) ? a : b;
}

inline __device__ __host__ float length2(const float3& a) {
  return (a.x * a.x + a.y * a.y + a.z * a.z);
}

inline double& elem(double3& a, uint i) {
  switch (i) {
    case 0:
      return a.x;
    case 1:
      return a.y;
    case 2:
      return a.z;
    default:
      return a.x;
  }
}

inline const double& elem(const double3& a, uint i) {
  switch (i) {
    case 0:
      return a.x;
    case 1:
      return a.y;
    case 2:
      return a.z;
    default:
      return a.x;
  }
}

inline __device__ __host__ double length(const double3& a) {
  return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

inline __device__ __host__ double length2(const double3& a) {
  return (a.x * a.x + a.y * a.y + a.z * a.z);
}

inline __device__ __host__ float distance2(const float3& a, const float3& b) {
  float x = a.x - b.x;
  float y = a.y - b.y;
  float z = a.z - b.z;
  return x * x + y * y + z * z;
}

inline __device__ __host__ double distance2(const double3& a,
                                            const double3& b) {
  double x = a.x - b.x;
  double y = a.y - b.y;
  double z = a.z - b.z;
  return x * x + y * y + z * z;
}

inline __device__ __host__ float distance(const float3& a, const float3& b) {
  return sqrtf(distance2(a, b));
}

inline __device__ __host__ double distance(const double3& a, const double3& b) {
  return sqrt(distance2(a, b));
}

inline __device__ __host__ uint index_from3d(const dim3& size,
                                             const dim3& pos) {
  return size.z * (pos.x * size.y + pos.y) + pos.z;
}

inline __device__ __host__ uint index_from4d(const uint4& size,
                                             const uint4& pos) {
  return size.w * (size.z * (size.y * pos.x + pos.y) + pos.z) + pos.w;
}

inline __device__ __host__ float sum(const float3& a) {
  return (a.x + a.y + a.z);
}

inline __device__ __host__ uint sum(const uint3& a) {
  return (a.x + a.y + a.z);
}

inline __device__ __host__ float3& operator-=(float3& a, const float& b) {
  a.x -= b;
  a.y -= b;
  a.z -= b;
  return a;
}

inline __device__ __host__ float3& operator-=(float3& a, const float3& b) {
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
  return a;
}

inline __device__ __host__ float4& operator-=(float4& a, const float4& b) {
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
  a.w -= b.w;
  return a;
}

inline __device__ __host__ float3& operator+=(float3& a, const float3& b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}

inline __device__ __host__ double4& operator+=(double4& a, const double4& b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  a.w += b.w;
  return a;
}

inline __device__ __host__ double4& operator-=(double4& a, const double4& b) {
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
  a.w -= b.w;
  return a;
}

inline __device__ __host__ double3& operator+=(double3& a, const double3& b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}

inline __device__ __host__ double3& operator-=(double3& a, const double3& b) {
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
  return a;
}

inline __device__ __host__ float3& operator+=(float3& a, const float& b) {
  a.x += b;
  a.y += b;
  a.z += b;
  return a;
}

inline __device__ __host__ float4& operator+=(float4& a, const float4& b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  a.w += b.w;
  return a;
}

inline __device__ __host__ dim3 operator/(const dim3& a, const uint b) {
  return dim3(a.x / b, a.y / b, a.z / b);
}

inline __device__ __host__ dim3 operator/(const dim3& a, const dim3& b) {
  return dim3(a.x / b.x, a.y / b.y, a.z / b.z);
}

inline __device__ __host__ dim3 operator%(const dim3& a, const dim3& b) {
  return dim3(a.x % b.x, a.y % b.y, a.z % b.z);
}

inline __device__ __host__ int divUp(int a, int b) {
  return ((a % b) != 0) ? (a / b + 1) : (a / b);
}

inline __device__ __host__ dim3 divUp(const dim3& a, const dim3& b) {
  return dim3(divUp(a.x, b.x), divUp(a.y, b.y), divUp(a.z, b.z));
}

inline __device__ __host__ uint4 operator+(const dim3& a, const uint4& b) {
  return make_uint4(a.x + b.x, a.y + b.y, a.z + b.z, b.w);
}

inline __device__ __host__ double3 operator*(const uint3& a, double b) {
  return make_double3(a.x * b, a.y * b, a.z * b);
}

inline __device__ __host__ double4 operator*(const double4& a, double b) {
  return make_double4(a.x * b, a.y * b, a.z * b, a.w * b);
}

inline __device__ __host__ float4 to_float4(const float3& f) {
  return make_float4(f.x, f.y, f.z, 0.0f);
}

inline __device__ __host__ float3 to_float3(const float4& f) {
  return make_float3(f.x, f.y, f.z);
}

inline __device__ __host__ double3 to_double3(const float3& f) {
  return make_double3(f.x, f.y, f.z);
}

inline __device__ __host__ float3 to_float3(const double3& f) {
  return make_float3(f.x, f.y, f.z);
}

inline __device__ __host__ float3 to_float3(const float3& f) { return f; }

#ifdef _DEBUG
inline void cudaAssertNoError(const char* msg = NULL) {
#if GPU_KERNELS

  cudaThreadSynchronize();
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    std::cerr << "CUDA ERROR: " << cudaGetErrorString(error) << " ["
              << (msg ? msg : "??") << "]" << std::endl;
    abort();
  }
#endif
}
#else
#define cudaAssertNoError(s)
#endif

inline int cudaGetGPUCount() {
  int devices = 0;
#if GPU_KERNELS
  if (cudaGetDeviceCount(&devices) != cudaSuccess)
    throw std::runtime_error("cudaGetDeviceCount failed");
#endif
  return devices;
}

inline void cudaGetMemoryInfo(size_t& free, size_t& total) {
#if GPU_KERNELS
  if (cuMemGetInfo(&free, &total) != CUDA_SUCCESS)
    throw std::runtime_error("cuMemGetInfo failed");
#endif
}

inline void cudaPrintMemoryInfo(void) {
#if GPU_KERNELS
  size_t free = 0, total = 0;
  cudaGetMemoryInfo(free, total);
  std::cout << "mem_used: " << (total - free) / (1024.0 * 1024.0)
            << "MB | mem_perc: "
            << ((double)(total - free) / (double)total) * 100.0 << "%"
            << std::endl;
#endif
}

template <typename T>
void to_constant(const std::string& name, const T* ptr) {
#if GPU_KERNELS
  cudaMemcpyToSymbol(name.c_str(), ptr, sizeof(T), 0, cudaMemcpyHostToDevice);
#endif
}

inline std::ostream& operator<<(std::ostream& o, const double2& a) {
  o << "(" << a.x << "," << a.y << ")";
  return o;
}
inline std::ostream& operator<<(std::ostream& o, const double3& a) {
  o << "(" << a.x << "," << a.y << "," << a.z << ")";
  return o;
}
inline std::ostream& operator<<(std::ostream& o, const double4& a) {
  o << "(" << a.x << "," << a.y << "," << a.z << ")";
  return o;
}
inline std::ostream& operator<<(std::ostream& o, const uint1& a) {
  o << "(" << a.x << ")";
  return o;
}
inline std::ostream& operator<<(std::ostream& o, const uint2& a) {
  o << "(" << a.x << "," << a.y << ")";
  return o;
}

inline std::ostream& operator<<(std::ostream& o, const float1& a) {
  o << "(" << a.x << ")";
  return o;
}
inline std::ostream& operator<<(std::ostream& o, const float2& a) {
  o << "(" << a.x << "," << a.y << ")";
  return o;
}
inline std::ostream& operator<<(std::ostream& o, const float3& a) {
  o << "(" << a.x << "," << a.y << "," << a.z << ")";
  return o;
}
inline std::ostream& operator<<(std::ostream& o, const float4& a) {
  o << "(" << a.x << "," << a.y << "," << a.z << "," << a.w << ")";
  return o;
}

#endif /* __CUTOOLS_H__ */
