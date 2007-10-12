#ifndef __CUTOOLS_H__
#define __CUTOOLS_H__

#include <cuda.h>

extern CUcontext CUDAContext;

#ifdef _DEBUG
#define MESSAGE(call) printf("[" #call "] @ %s %i\n", __FILE__, __LINE__);
#else
#define MESSAGE(call)
#endif

#ifdef NOCUDA
#define CUDA_SAFE_CALL( call, exc, mess, RETURN ) { } while(0)
#define CU_SAFE_CALL( call, exc, mess, RETURN ) { } while(0)
#else
#define CUDA_SAFE_CALL( call, exc, mess, RETURN ) do {                       \
    MESSAGE(#call);                                                          \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        PyErr_Format( exc, mess                                              \
	              ":CUDA error in file '%s' in line %i : %s.",           \
                      __FILE__, __LINE__, cudaGetErrorString( err ) );       \
        RETURN;                                                              \
    } } while (0)
#define CU_SAFE_CALL( call, exc, mess, RETURN ) do {                         \
    MESSAGE(#call);                                                          \
    CUresult err = call;                                                     \
    if( CUDA_SUCCESS != err) {                                               \
        PyErr_Format( exc, mess                                              \
	              ":CUDA driver error in file '%s' in line %i.\n",       \
                      __FILE__, __LINE__ );                                  \
        RETURN;                                                              \
    } } while (0)
#endif

inline __device__ __host__ float dot(const float4 a, const float4 b)
{ return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

inline __device__ __host__ float4 dot(const float4 a[4], const float4 b)
{ return make_float4( dot(a[0], b), dot(a[1], b), dot(a[2], b), dot(a[3], b) ); }

inline __device__ __host__ float2 Cmul(float2 a, float2 b)
{ return make_float2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }

inline __device__ __host__ float2 Csum(float2 a, float2 b)
{ return make_float2(a.x + b.x, a.y + a.x); }

inline __device__ __host__ float2 Cscale(float2 a, float s)
{ return make_float2(a.x * s, a.y * s); }

inline __device__ __host__ float2 operator -(const float2 a)
{ return make_float2(-a.x, -a.y); }

inline __device__ __host__ float3 operator *(const float3 a, const float3 b)
{ return make_float3(a.x*b.x, a.y*b.y, a.z*b.z); }

inline __device__ __host__ float3 operator /(const float3 a, const float3 b)
{ return make_float3(a.x/b.x, a.y/b.y, a.z/b.z); }

inline __device__ __host__ float3 operator +(const float3 a, const float3 b)
{ return make_float3(a.x+b.x, a.y+b.y, a.z+b.z); }

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

inline __device__ __host__ uint3 index(const dim3 bd, const uint3 bi, const uint3 ti)
{ return bd * bi + ti; }

inline __device__ __host__ bool operator < (const uint3 i, const uint3 r)
{ return i.x < r.x && i.y < r.y && i.z < r.z; }

inline __device__ __host__ uint reference(const uint3 i, const uint3 g)
{ return i.x + i.y * g.x + i.z * g.x * g.y; }

inline __device__ __host__ uint len2(const float3 a)
{ return a.x*a.x + a.y*a.y + a.z*a.z; }

inline __device__ __host__ void swap(float4 &a, float4 &b)
{ float4 tmp = a; a = b; b = tmp; }

inline __device__ __host__ void negate(float2 &a)
{ a.x = -a.x; a.y = -a.y; }

inline int iDivUp(int a, int b)
{ return ((a % b) != 0) ? (a / b + 1) : (a / b); }

inline __device__ __host__ const float4& min(const float4 &a, const float4 &b) 
{ return (a.x < b.x)? a: b; }

inline __device__ __host__ const float4& max(const float4 &a, const float4 &b) 
{ return (a.x > b.x)? a: b; }

#endif /* __CUTOOLS_H__ */

