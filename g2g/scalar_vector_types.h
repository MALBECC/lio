#ifndef _SCALAR_VECTOR_TYPES_H
#define	_SCALAR_VECTOR_TYPES_H

#include "datatypes/cpu_primitives.h"

namespace G2G {
  template<class T, unsigned int n> class vec_type { };

#if CPU_KERNELS

  template<> class vec_type<float, 2> : public float2 {
    public:
      vec_type(void) {}
      vec_type(float _x, float _y) { this->x = _x; this->y = _y; }
      typedef float2 base_type;
      //operator float2 () { return (float2)(*this); }
  };

  template<> class vec_type<float, 3> {
  	private:
  	  float _x, _y, _z;

    public:
      vec_type(void) {}
      vec_type(const float3& other) : _x(other.x),_y(other.y),_z(other.z) { }
      vec_type(const double3& other) : _x(other.x),_y(other.y),_z(other.z) { }
      vec_type(float x, float y, float z) : _x(x), _y(y), _z(z) { }
      
      inline float x() const { return _x; }
      inline float y() const { return _y; }
      inline float z() const { return _z; }
      inline float length2(){ return _x * _x + _y * _y + _z * _z; }

      friend bool isinf(const vec_type & e){ return isinf(e._x) || isinf(e._y) || isinf(e._z); }
      friend bool isnan(const vec_type & e){ return isnan(e._x) || isnan(e._y) || isnan(e._z); }

      friend vec_type operator*(const vec_type & lo, const vec_type & ro) {
      	return vec_type(lo.x() * ro.x(), lo.y() * ro.y(), lo.z() * ro.z());
      }

      friend vec_type operator-(const vec_type & lo, const vec_type & ro) {
      	return vec_type(lo.x() - ro.x(), lo.y() - ro.y(), lo.z() - ro.z());
      }

      friend vec_type operator+(const vec_type & lo, const vec_type & ro) {
      	return vec_type(lo.x() + ro.x(), lo.y() + ro.y(), lo.z() + ro.z());
      }

      friend std::ostream& operator<<(std::ostream& o, const vec_type & v) {
          o << v.x() << " " << v.y() << " " << v.z();
      }

      void operator+=(const vec_type & lo){
      	_x += lo.x(), _y += lo.y(), _z += lo.z();
      }
    
      void operator-=(const vec_type & lo){
      	_x -= lo.x(), _y -= lo.y(), _z -= lo.z();
      }

      friend vec_type operator*(const vec_type & lo, float f) {
      	return vec_type(lo.x()*f,lo.y()*f,lo.z()*f);
      }

      friend vec_type operator-(const vec_type & lo, float f) {
      	return vec_type(lo.x()-f,lo.y()-f,lo.z()-f);
      }

      typedef float3 base_type;
  };

  template<> class vec_type<float, 4> {
    private:
      float _x, _y, _z, _w;

    public:
      vec_type(void) {}
      vec_type(float x, float y, float z, float w) : _x(_x), _y(y), _z(z), _w(w) { }
      typedef float4 base_type;

      friend bool isinf(const vec_type & e){ return false; }
      friend bool isnan(const vec_type & e){ return false; }

      friend std::ostream& operator<<(std::ostream& o, const vec_type & v) {
          o << v._x << " " << v._y << " " << v._z << " " << v._w;
      }

  };
#else

  template<> class vec_type<float, 2> : public float2 {
    public:
      __device__ __host__ vec_type(void) {}
      __device__ __host__ vec_type(float _x, float _y) { this->x = _x; this->y = _y; }
      typedef float2 base_type;
      //operator float2 () { return (float2)(*this); }
  };
  template<> class vec_type<float, 3> : public float3 {
    public:
      __device__ __host__ vec_type(void) {}
      __device__ __host__ vec_type(const float4& other) { float3::x = other.x; float3::y = other.y; float3::z = other.z; }
      __device__ __host__ vec_type(const float3& other) { float3::x = other.x; float3::y = other.y; float3::z = other.z; }
      __device__ __host__ explicit vec_type(const double3& other) { float3::x = other.x; float3::y = other.y; float3::z = other.z; }
      __device__ __host__ vec_type(float _x, float _y, float _z) { float3::x = _x; float3::y = _y; float3::z = _z; }
      typedef float3 base_type;
     /*
      __device__ __host__ inline double x(void) const { return float3::x; }
      __device__ __host__ inline double y(void) const { return float3::y; }
      __device__ __host__ inline double z(void) const { return float3::z; }
      */

  };

  template<> class vec_type<float, 4> : public float4 {
    public:
      __device__ __host__ vec_type(void) {}
      __device__ __host__ vec_type(const float3& other) { float4::x = other.x; float4::y = other.y; float4::z = other.z; float4::w = 0.0f; }
      __device__ __host__ vec_type(const float4& other) { this->x = other.x; this->y = other.y; this->z = other.z; this->w = other.w; }
      __device__ __host__ vec_type(float _x, float _y, float _z, float _w) { this->x = _x; this->y = _y; this->z = _z; this->w = _w; }
      typedef float4 base_type;
  };

  #endif

  #if CPU_KERNELS

  template<> class vec_type<double, 2> : public double2 {
    public:
      vec_type(void) {}
      vec_type(double _x, double _y) { this->x = _x; this->y = _y; }
      typedef double2 base_type;
      // operator double2 () { return (double2)(*this); }
  };

  template<> class vec_type<double, 3> : public double3 {
    public:
      vec_type(void) {}
      vec_type(const double3& other) : double3(other) { }
      explicit vec_type(const float3& other) { double3::x = other.x; double3::y = other.y; double3::z = other.z; }
      explicit vec_type(const double4& other) { double3::x = other.x; double3::y = other.y; double3::z = other.z; }
      vec_type(double _x, double _y, double _z) { double3::x = _x; double3::y = _y; double3::z = _z; }
      typedef double3 base_type;

      double length2(void) { return double3::x * double3::x + double3::y * double3::y + double3::z * double3::z; }
      inline double x(void) const { return double3::x; }
      inline double y(void) const { return double3::y; }
      inline double z(void) const { return double3::z; }
      //operator double3 () { return (double3)(*this); }
  };

  template<> class vec_type<double, 4> : public double4 {
    public:
      vec_type(void) {}
      explicit vec_type(const double3& other) { this->x = other.x; this->y = other.y; this->z = other.z; }
      explicit vec_type(const double4& other) { this->x = other.x; this->y = other.y; this->z = other.z; }
      vec_type(double _x, double _y, double _z, double _w) { this->x = _x; this->y = _y; this->z = _z; this->w = _w; }
      //operator double4 () { return (double4)(*this); }
  };


#else

  template<> class vec_type<double, 2> : public double2 {
    public:
      __device__ __host__ vec_type(void) {}
      __device__ __host__ vec_type(double _x, double _y) { this->x = _x; this->y = _y; }
      typedef double2 base_type;
      // operator double2 () { return (double2)(*this); }
  };


  template<> class vec_type<double, 3> : public double3 {
    public:
      __device__ __host__ vec_type(void) {}
      __device__ __host__ vec_type(const double3& other) : double3(other) { }
      __device__ __host__ explicit vec_type(const float3& other) { double3::x = other.x; double3::y = other.y; double3::z = other.z; }
      __device__ __host__ explicit vec_type(const double4& other) { double3::x = other.x; double3::y = other.y; double3::z = other.z; }
      __device__ __host__ vec_type(double _x, double _y, double _z) { double3::x = _x; double3::y = _y; double3::z = _z; }
      typedef double3 base_type;

      __device__ __host__ double length2(void) { return double3::x * double3::x + double3::y * double3::y + double3::z * double3::z; }
      /*
      __device__ __host__ inline double x(void) const { return double3::x; }
      __device__ __host__ inline double y(void) const { return double3::y; }
      __device__ __host__ inline double z(void) const { return double3::z; }
      */
      //operator double3 () { return (double3)(*this); }
  };

  template<> class vec_type<double, 4> : public double4 {
    public:
      __device__ __host__ vec_type(void) {}
      __device__ __host__ explicit vec_type(const double3& other) { this->x = other.x; this->y = other.y; this->z = other.z; }
      __device__ __host__ explicit vec_type(const double4& other) { this->x = other.x; this->y = other.y; this->z = other.z; }
      __device__ __host__ vec_type(double _x, double _y, double _z, double _w) { this->x = _x; this->y = _y; this->z = _z; this->w = _w; }
      //operator double4 () { return (double4)(*this); }
  };

  template<> class vec_type<int2, 4>{
    public:
      __device__ __host__ vec_type(void) {}
      __device__ __host__ vec_type(int2 _x, int2 _y, int2 _z, int2 _w) { this->x = _x; this->y = _y; this->z = _z; this->w = _w; }
      int2 x, y, z, w;
  };
#endif
}

#endif	/* _SCALAR_VECTOR_TYPES_H */

