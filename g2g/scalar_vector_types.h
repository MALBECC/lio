#ifndef _SCALAR_VECTOR_TYPES_H
#define _SCALAR_VECTOR_TYPES_H

#include "datatypes/cpu_primitives.h"

namespace G2G {
template <class T, unsigned int n>
class vec_type {};

#if !GPU_KERNELS

template <>
class vec_type<float, 2> : public float2 {
 public:
  vec_type(void) {}
  vec_type(float _x, float _y) {
    this->x = _x;
    this->y = _y;
  }
  // operator float2 () { return (float2)(*this); }
};

template <>
class vec_type<float, 3> {
 public:
  float x, y, z;

  vec_type(void) {}
  vec_type(const float3& other) : x(other.x), y(other.y), z(other.z) {}
  vec_type(const double3& other) : x(other.x), y(other.y), z(other.z) {}
  vec_type(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}

  inline float length2() const { return x * x + y * y + z * z; }

  friend vec_type operator*(const vec_type& lo, const vec_type& ro) {
    return vec_type(lo.x * ro.x, lo.y * ro.y, lo.z * ro.z);
  }

  friend vec_type operator-(const vec_type& lo, const vec_type& ro) {
    return vec_type(lo.x - ro.x, lo.y - ro.y, lo.z - ro.z);
  }

  friend vec_type operator+(const vec_type& lo, const vec_type& ro) {
    return vec_type(lo.x + ro.x, lo.y + ro.y, lo.z + ro.z);
  }

  friend std::ostream& operator<<(std::ostream& o, const vec_type& v) {
    o << v.x << " " << v.y << " " << v.z;
    return o;
  }

  void operator+=(const vec_type& lo) { x += lo.x, y += lo.y, z += lo.z; }

  void operator-=(const vec_type& lo) { x -= lo.x, y -= lo.y, z -= lo.z; }

  friend vec_type operator*(const vec_type& lo, float f) {
    return vec_type(lo.x * f, lo.y * f, lo.z * f);
  }

  friend vec_type operator-(const vec_type& lo, float f) {
    return vec_type(lo.x - f, lo.y - f, lo.z - f);
  }
};

template <>
class vec_type<float, 4> {
 private:
  //float x, y, z, _w;

 public:
  /** La integracion a Libxc necesita acceder a estos componentes **/
  float x, y, z, _w;

  vec_type(void) {}
  vec_type(float x_in, float y_in, float z_in, float w_in)
      : x(x_in), y(y_in), z(z_in), _w(w_in) {}

  friend std::ostream& operator<<(std::ostream& o, const vec_type& v) {
    o << v.x << " " << v.y << " " << v.z << " " << v._w;
    return o;
  }
};
#else

template <>
class vec_type<float, 2> : public float2 {
 public:
  __device__ __host__ vec_type(void) {}
  __device__ __host__ vec_type(float _x, float _y) {
    this->x = _x;
    this->y = _y;
  }
  typedef float2 base_type;
};
template <>
class vec_type<float, 3> : public float3 {
 public:
  __device__ __host__ vec_type(void) {}
  __device__ __host__ vec_type(const float4& other) {
    float3::x = other.x;
    float3::y = other.y;
    float3::z = other.z;
  }
  __device__ __host__ vec_type(const float3& other) {
    float3::x = other.x;
    float3::y = other.y;
    float3::z = other.z;
  }
  __device__ __host__ explicit vec_type(const double3& other) {
    float3::x = other.x;
    float3::y = other.y;
    float3::z = other.z;
  }
  __device__ __host__ vec_type(float _x, float _y, float _z) {
    float3::x = _x;
    float3::y = _y;
    float3::z = _z;
  }
  __device__ __host__ float length2() const {
    return float3::x * float3::x + float3::y * float3::y +
           float3::z * float3::z;
  }
  typedef float3 base_type;
};

template <>
class vec_type<float, 4> : public float4 {
 public:
  __device__ __host__ vec_type(void) {}
  __device__ __host__ vec_type(const float3& other) {
    float4::x = other.x;
    float4::y = other.y;
    float4::z = other.z;
    float4::w = 0.0f;
  }
  __device__ __host__ vec_type(const float4& other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->w = other.w;
  }
  __device__ __host__ vec_type(float _x, float _y, float _z, float _w) {
    this->x = _x;
    this->y = _y;
    this->z = _z;
    this->w = _w;
  }
  typedef float4 base_type;
};

#endif

#if !GPU_KERNELS

template <>
class vec_type<double, 2> : public double2 {
 public:
  vec_type(void) {}
  vec_type(double _x, double _y) {
    this->x = _x;
    this->y = _y;
  }
  typedef double2 base_type;
};

template <>
class vec_type<double, 3> : public double3 {
 public:
  vec_type(void) {}
  vec_type(const double3& other) : double3(other) {}
  explicit vec_type(const float3& other) {
    double3::x = other.x;
    double3::y = other.y;
    double3::z = other.z;
  }
  explicit vec_type(const double4& other) {
    double3::x = other.x;
    double3::y = other.y;
    double3::z = other.z;
  }
  vec_type(double _x, double _y, double _z) {
    double3::x = _x;
    double3::y = _y;
    double3::z = _z;
  }
  typedef double3 base_type;

  double length2(void) const {
    return double3::x * double3::x + double3::y * double3::y +
           double3::z * double3::z;
  }
};

template <>
class vec_type<double, 4> : public double4 {
 public:
  vec_type(void) {}
  explicit vec_type(const double3& other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
  }
  explicit vec_type(const double4& other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
  }
  vec_type(double _x, double _y, double _z, double _w) {
    this->x = _x;
    this->y = _y;
    this->z = _z;
    this->w = _w;
  }
};

#else

template <>
class vec_type<double, 2> : public double2 {
 public:
  __device__ __host__ vec_type(void) {}
  __device__ __host__ vec_type(double _x, double _y) {
    this->x = _x;
    this->y = _y;
  }
  typedef double2 base_type;
};

template <>
class vec_type<double, 3> : public double3 {
 public:
  __device__ __host__ vec_type(void) {}
  __device__ __host__ vec_type(const double3& other) : double3(other) {}
  __device__ __host__ explicit vec_type(const float3& other) {
    double3::x = other.x;
    double3::y = other.y;
    double3::z = other.z;
  }
  __device__ __host__ explicit vec_type(const double4& other) {
    double3::x = other.x;
    double3::y = other.y;
    double3::z = other.z;
  }
  __device__ __host__ vec_type(double _x, double _y, double _z) {
    double3::x = _x;
    double3::y = _y;
    double3::z = _z;
  }
  typedef double3 base_type;
  __device__ __host__ double length2(void) const {
    return double3::x * double3::x + double3::y * double3::y +
           double3::z * double3::z;
  }
};

template <>
class vec_type<double, 4> : public double4 {
 public:
  __device__ __host__ vec_type(void) {}
  __device__ __host__ explicit vec_type(const double3& other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
  }
  __device__ __host__ explicit vec_type(const double4& other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
  }
  __device__ __host__ vec_type(double _x, double _y, double _z, double _w) {
    this->x = _x;
    this->y = _y;
    this->z = _z;
    this->w = _w;
  }
};

template <>
class vec_type<int2, 4> {
 public:
  __device__ __host__ vec_type(void) {}
  __device__ __host__ vec_type(int2 _x, int2 _y, int2 _z, int2 _w) {
    this->x = _x;
    this->y = _y;
    this->z = _z;
    this->w = _w;
  }
  int2 x, y, z, w;
};
#endif
}

#endif /* _SCALAR_VECTOR_TYPES_H */
