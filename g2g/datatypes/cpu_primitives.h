#ifndef _CPU_PRIMITIVES
#define _CPU_PRIMITIVES
#if (CPU_KERNELS && !GPU_KERNELS)

typedef unsigned int uint;
typedef struct uint1 { uint x; } uint1;
typedef struct uint2 { uint x, y; } uint2;
typedef struct uint3 { uint x, y, z; } uint3;
typedef struct uint4 { uint x, y, z, w; } uint4;
typedef struct int1 { int x; } int1;
typedef struct int2 { int x, y; } int2;
typedef struct int3 { int x, y, z; } int3;
typedef struct int4 { int x, y, z, w; } int4;
typedef struct float1 { float x; } float1;
typedef struct float2 { float x, y; } float2;
typedef struct float3 { float x, y, z; } float3;
typedef struct float4 { float x, y, z, w; } float4;
typedef struct double1 { double x; } double1;
typedef struct double2 { double x, y; } double2;
typedef struct double3 { double x, y, z; } double3;
typedef struct double4 { double x, y, z, w; } double4;

static int1 make_int1(int x) {
  int1 t;
  t.x = x;
  return t;
}

static uint1 make_uint1(unsigned int x) {
  uint1 t;
  t.x = x;
  return t;
}

static int2 make_int2(int x, int y) {
  int2 t;
  t.x = x;
  t.y = y;
  return t;
}

static uint2 make_uint2(unsigned int x, unsigned int y) {
  uint2 t;
  t.x = x;
  t.y = y;
  return t;
}

static int3 make_int3(int x, int y, int z) {
  int3 t;
  t.x = x;
  t.y = y;
  t.z = z;
  return t;
}

static uint3 make_uint3(unsigned int x, unsigned int y, unsigned int z) {
  uint3 t;
  t.x = x;
  t.y = y;
  t.z = z;
  return t;
}

static int4 make_int4(int x, int y, int z, int w) {
  int4 t;
  t.x = x;
  t.y = y;
  t.z = z;
  t.w = w;
  return t;
}

static uint4 make_uint4(unsigned int x, unsigned int y, unsigned int z,
                        unsigned int w) {
  uint4 t;
  t.x = x;
  t.y = y;
  t.z = z;
  t.w = w;
  return t;
}

static float1 make_float1(float x) {
  float1 t;
  t.x = x;
  return t;
}

static float2 make_float2(float x, float y) {
  float2 t;
  t.x = x;
  t.y = y;
  return t;
}

static float3 make_float3(float x, float y, float z) {
  float3 t;
  t.x = x;
  t.y = y;
  t.z = z;
  return t;
}

static float4 make_float4(float x, float y, float z, float w) {
  float4 t;
  t.x = x;
  t.y = y;
  t.z = z;
  t.w = w;
  return t;
}

static double1 make_double1(double x) {
  double1 t;
  t.x = x;
  return t;
}

static double2 make_double2(double x, double y) {
  double2 t;
  t.x = x;
  t.y = y;
  return t;
}

static double3 make_double3(double x, double y, double z) {
  double3 t;
  t.x = x;
  t.y = y;
  t.z = z;
  return t;
}

static double4 make_double4(double x, double y, double z, double w) {
  double4 t;
  t.x = x;
  t.y = y;
  t.z = z;
  t.w = w;
  return t;
}

typedef struct dim3 {
  unsigned int x, y, z;
  dim3(unsigned int vx = 1, unsigned int vy = 1, unsigned int vz = 1)
      : x(vx), y(vy), z(vz) {}
  dim3(uint3 v) : x(v.x), y(v.y), z(v.z) {}
  operator uint3(void) {
    uint3 t;
    t.x = x;
    t.y = y;
    t.z = z;
    return t;
  }
} dim3;

#endif
#endif
