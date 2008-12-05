#ifndef __CUTOOLS_H__
#define __CUTOOLS_H__

#include <cuda_runtime.h>
#include <cuda.h>

#ifdef __DEVICE_EMULATION__
#define _EMU(code) code
#else
#define _EMU(code)
#endif

#define WARP_SIZE 32
#define BANKS 16

#include <stdexcept>
#include <iostream>
#include "double.h"

/** operators **/
inline __device__ __host__ float2 operator -(const float2 a)
{ return make_float2(-a.x, -a.y); }

inline __device__ __host__ float3 operator *(const float3& a, const float3& b)
{ return make_float3(a.x*b.x, a.y*b.y, a.z*b.z); }

inline __device__ __host__ float3 operator /(const float3& a, const float3& b)
{ return make_float3(a.x/b.x, a.y/b.y, a.z/b.z); }

inline __device__ __host__ float3 operator +(const float3& a, const float3& b)
{ return make_float3(a.x+b.x, a.y+b.y, a.z+b.z); }

inline __device__ __host__ float3 operator +(const float3& a, const float& b)
{ return make_float3(a.x+b, a.y+b, a.z+b); }

inline __device__ __host__ float3 operator -(const float3 a, const float3 b)
{ return make_float3(a.x-b.x, a.y-b.y, a.z-b.z); }

inline __device__ __host__ float3 operator +(const uint3 a, const float3 b)
{ return make_float3(a.x+b.x, a.y+b.y, a.z+b.z); }

inline __device__ __host__ uint3 operator *(const uint3 a, const uint3 b)
{ return make_uint3(a.x*b.x, a.y*b.y, a.z*b.z); }

inline __device__ __host__ uint3 operator *(const dim3 a, const uint3 b)
{ return make_uint3(a.x*b.x, a.y*b.y, a.z*b.z); }

inline __device__ __host__ uint3 operator +(const uint3 a, const uint3 b)
{ return make_uint3(a.x+b.x, a.y+b.y, a.z+b.z); }

inline __device__ __host__ bool operator < (const uint3 i, const uint3 r)
{ return i.x < r.x && i.y < r.y && i.z < r.z; }

inline __device__ __host__ float3 operator*(const float3& a, float b)
{ return make_float3(a.x * b, a.y * b, a.z * b); }

inline __device__ __host__ float1 operator*(const float1 &a, float b)
{ return make_float1(a.x * b); }

inline __device__ __host__ uint3 index(const dim3 bd, const uint3 bi, const uint3 ti)
{ return bd * bi + ti; }

inline __device__ __host__ uint reference(const uint3 i, const uint3 g)
{ return i.x + i.y * g.x + i.z * g.x * g.y; }

inline __device__ __host__ uint len2(const float3 a)
{ return (uint)a.x * (uint)a.x + (uint)a.y * (uint)a.y + (uint)a.z * (uint)a.z; }

inline __device__ __host__ void swap(float4 &a, float4 &b)
{ float4 tmp = a; a = b; b = tmp; }

inline __device__ __host__ void negate(float2 &a)
{ a.x = -a.x; a.y = -a.y; }

inline __device__ __host__ const float4& min(const float4 &a, const float4 &b) 
{ return (a.x < b.x)? a: b; }

inline __device__ __host__ const float4& max(const float4 &a, const float4 &b) 
{ return (a.x > b.x)? a: b; }

inline __device__ __host__ float length2(const float3& a) {
	return (a.x * a.x + a.y * a.y + a.z * a.z);
}

inline __device__ __host__ float distance2(const float3& a, const float3& b)
{
	float x = a.x - b.x;
	float y = a.y - b.y;
	float z = a.z - b.z;
	return x * x + y * y + z * z;
}

inline __device__ __host__ float distance(const float3& a, const float3& b)
{ return sqrtf(distance2(a, b)); }

inline __device__ __host__ uint index_from3d(const dim3& size, const dim3& pos)
{
	return size.z * (pos.x * size.y + pos.y) + pos.z;
}

inline __device__ __host__ uint index_from4d(const uint4& size, const uint4& pos)
{
	return size.w * (size.z * (size.y * pos.x + pos.y) + pos.z) + pos.w;
}

inline __device__ __host__ const float& float3_elem(const float3& a, uint i)
{
	switch(i) {
		case 0: return a.x;
		case 1: return a.y;
		case 2: return a.z;
		default: return a.x;
	}
}

inline __device__ __host__ float& float3_elem(float3& a, uint i)
{
	switch(i) {
		case 0: return a.x;
		case 1: return a.y;
		case 2: return a.z;
		default: return a.x;
	}
}

inline __device__ __host__ float sum(const float3& a)
{ return (a.x + a.y + a.z); }

inline __device__ __host__ uint sum(const uint3& a)
{ return (a.x + a.y + a.z); }

inline __device__ __host__ float3& operator-=(float3& a, const float& b)
{ a.x -= b; a.y -= b; a.z -= b; return a; }

inline __device__ __host__ float3& operator-=(float3& a, const float3& b)
{ a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }

inline __device__ __host__ float3& operator+=(float3& a, const float3& b)
{ a.x += b.x; a.y += b.y; a.z += b.z; return a; }

inline __device__ __host__ float3& operator+=(float3& a, const float& b)
{ a.x += b; a.y += b; a.z += b; return a; }

inline __device__ __host__ dim3 operator/(const dim3& a, const uint b)
{ return dim3(a.x / b, a.y / b, a.z / b); }

inline __device__ __host__ dim3 operator/(const dim3& a, const dim3& b)
{ return dim3(a.x / b.x, a.y / b.y, a.z / b.z); }

inline __device__ __host__ dim3 operator%(const dim3& a, const dim3& b)
{ return dim3(a.x % b.x, a.y % b.y, a.z % b.z); }

inline __device__ __host__ int divUp(int a, int b)
{ return ((a % b) != 0) ? (a / b + 1) : (a / b); }

inline __device__ __host__ dim3 divUp(const dim3& a, const dim3& b)
{ return dim3(divUp(a.x, b.x), divUp(a.y, b.y), divUp(a.z, b.z)); }

inline __device__ __host__ uint4 operator +(const dim3& a, const uint4& b)
{ return make_uint4(a.x + b.x, a.y + b.y, a.z + b.z, b.w); }

inline void cudaAssertNoError(const char* msg = NULL) {
#ifdef _DEBUG
	cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) {
		std::cerr << "CUDA ERROR: " << cudaGetErrorString(error) << " " << (msg ? msg : "??") << std::endl;
		abort();
	}
#endif
}

inline void cudaPrintMemoryInfo(void) {
	uint free = 0, total = 0;
	CUresult res = cuMemGetInfo(&free, &total);
	if (res != CUDA_SUCCESS) {
		throw std::runtime_error("cuMemGetInfo");
	}
	std::cout << "free: " << (free / (1024.0 * 1024.0)) << "/" << (total / (1024.0 * 1024.0)) << " used: " << (total - free) / (1024.0 * 1024.0) << " (" << ((double)free / (double)total) * 100.0 << "%)" << std::endl;
}

#endif /* __CUTOOLS_H__ */

