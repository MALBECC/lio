#ifndef __G2G_CPU_VECTOR_TYPES_H__
#define __G2G_CPU_VECTOR_TYPES_H__

#ifndef __CUDACC__

#include <math.h>
#include <sys/types.h>
#include <fvec.h>
#include <xmmintrin.h>

#include "cuda/cuda_extra.h"

namespace G2G {
	class cfloat4 : public F32vec4 {
	  public:
      cfloat4(void) : F32vec4() { }
      explicit cfloat4(float x) : F32vec4(x) { }
    	cfloat4(F32vec4 a) : F32vec4(a) {  }
      explicit cfloat4(float4 a) : F32vec4(a.w, a.z, a.y, a.x) { }
  		explicit cfloat4(float x, float y, float z, float w) : F32vec4(w, z, y, x) { }

  		inline const float& x(void) const { return (*this)[0]; }
  		inline const float& y(void) const { return (*this)[1]; }
  		inline const float& z(void) const { return (*this)[2]; }
  		inline const float& w(void) const { return (*this)[3]; }

      inline float& x(void) { return (*this)[0]; }
  		inline float& y(void) { return (*this)[1]; }
  		inline float& z(void) { return (*this)[2]; }
  		inline float& w(void) { return (*this)[3]; }
	
    	friend std::ostream& operator<<(std::ostream & os, const cfloat4& a)
  	  {                                                                                                                                           
  		  float *fp = (float*)&a;
  		  os << "(" << *fp << "," << *(fp+1) << "," << *(fp+2) << "," << *(fp+3) << ")";
        return os;
    	}

      inline operator float4() { return make_float4(x(), y(), z(), w()); }
  };

  inline bool isinf(cfloat4 v) { return isinff(v.x()) || isinff(v.y()) || isinff(v.z()) || isinff(v.w()); }
  inline bool isnan(cfloat4 v) { return isnanf(v.x()) || isnanf(v.y()) || isnanf(v.z()) || isnanf(v.w()); }
	
	class cfloat3 : public cfloat4 {
		public:
      cfloat3(void) : cfloat4() { }
      explicit cfloat3(float x) : cfloat4(x) { }
    	cfloat3(F32vec4 a) : cfloat4(a) {  }
      explicit cfloat3(float3 a) : cfloat4(a.x, a.y, a.z, 0.0f) { }
		  explicit cfloat3(float x, float y, float z) : cfloat4(x, y, z, 0.0f) { }

      inline operator float3() { return make_float3(x(), y(), z()); }

      friend cfloat3 operator*(const cfloat3& a, float b) { return a * cfloat3(b, b, b); }

      friend std::ostream& operator<<(std::ostream & os, const cfloat3& a)
  	  {
  		  float *fp = (float*)&a;
  		  os << "(" << *fp << "," << *(fp+1) << "," << *(fp+2) << ")";
        return os;
    	}
		
		private:
		  cfloat3(float x, float y, float z, float w);
	};

  inline bool isinf(cfloat3 v) { return isinff(v.x()) || isinff(v.y()) || isinff(v.z()); }
  inline bool isnan(cfloat3 v) { return isnanf(v.x()) || isnanf(v.y()) || isnanf(v.z()); }
}

#endif

#endif
