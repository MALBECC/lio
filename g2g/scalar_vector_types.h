#ifndef _SCALAR_VECTOR_TYPES_H
#define	_SCALAR_VECTOR_TYPES_H

namespace G2G {
  template<class T, unsigned int n> class vec_type { };
  template<> class vec_type<float, 2> : public float2 {
    public:
      __device__ __host__ vec_type(void) {}      
      __device__ __host__ vec_type(float _x, float _y) { this->x = _x; this->y = _y; }
      typedef float2 base_type;
      //operator float2 () { return (float2)(*this); }
  };
  
#if CPU_KERNELS
  template<> class vec_type<float, 3> : public cfloat3 {
    public:
      vec_type(void) {}
      vec_type(const cfloat3& other) : cfloat3(other) { }
      vec_type(const double3& other) : cfloat3(other.x, other.y, other.z) { }
      vec_type(float _x, float _y, float _z) : cfloat3(_x, _y, _z) { }
      typedef cfloat3 base_type;
      //operator cfloat3 () { return (cfloat3)(*this); }
  };
  
  template<> class vec_type<float, 4> : public cfloat4 {
    public:
      vec_type(void) {}      
      vec_type(float _x, float _y, float _z, float _w) : cfloat4(_w, _y, _z, _w) { }
      typedef cfloat4 base_type;      
      //operator cfloat4 () { return (cfloat4)(*this); }
  };
#else
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
  
  template<> class vec_type<double, 2> : public double2 {
    public:
      __device__ __host__ vec_type(void) {}      
      __device__ __host__ vec_type(double _x, double _y) { this->x = _x; this->y = _y; }
      typedef double2 base_type;      
      // operator double2 () { return (double2)(*this); }
  };
  
#if CPU_KERNELS
  template<> class vec_type<double, 3> : public double3 {
    public:
      __device__ __host__ vec_type(void) {}
      __device__ __host__ vec_type(const double3& other) : double3(other) { }
      __device__ __host__ explicit vec_type(const float3& other) { double3::x = other.x; double3::y = other.y; double3::z = other.z; }
      __device__ __host__ explicit vec_type(const double4& other) { double3::x = other.x; double3::y = other.y; double3::z = other.z; }
      __device__ __host__ vec_type(double _x, double _y, double _z) { double3::x = _x; double3::y = _y; double3::z = _z; }
      typedef double3 base_type;      
      
      __device__ __host__ double length2(void) { return double3::x * double3::x + double3::y * double3::y + double3::z * double3::z; }    
      __device__ __host__ inline double x(void) const { return double3::x; }
      __device__ __host__ inline double y(void) const { return double3::y; }
      __device__ __host__ inline double z(void) const { return double3::z; }
      //operator double3 () { return (double3)(*this); }
  };
#else

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
#endif  
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
}

#endif	/* _SCALAR_VECTOR_TYPES_H */

