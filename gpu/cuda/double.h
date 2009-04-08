#ifndef __DOUBLE_H__
#define __DOUBLE_H__

#include <cuda_runtime.h>
#include <cmath>
#include <stdexcept>

struct double3 {
  union {
    struct {
      double x, y, z;
    };
    double elements[3];
  };
  
  double3(void) : x(0),y(0),z(0) {  }

	double3(double _x, double _y, double _z) : x(_x),y(_y),z(_z) { }
	
	inline double length(void) const {
		return sqrt(x * x + y * y + z * z);
	}
	
	inline double length_square(void) const {
		return x * x + y * y + z * z;
	}

  inline double& operator[](unsigned int axis) {
    if (axis > 3) throw std::runtime_error("Axis number out of bounds");
    return elements[axis];
  }

  inline const double& operator[](unsigned int axis) const {
    if (axis > 3) throw std::runtime_error("Axis number out of bounds");
    return elements[axis];
  }

  inline double3& operator=(const double3& other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& o, const double3& v);
};

inline std::ostream& operator<<(std::ostream& o, const double3& v) {
  o << "(" << v.x << "," << v.y << "," << v.z << ")";
  return o;
}

inline double3 operator-(const double3& a, const double3& b) {
	return double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline double3 operator+(const double3& a, const double3& b) {
	return double3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline double3 operator*(const double3& a, const double3& b) {
	return double3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline double3 operator*(const double3& a, const double& b) {
	return double3(a.x * b, a.y * b, a.z * b);
}

inline double3 operator/(const double3& a, const uint b) {
	return double3(a.x / b, a.y / b, a.z / b);
}

inline double3 operator/(const double3& a, const double b) {
	return double3(a.x / b, a.y / b, a.z / b);
}


inline uint3 ceil_uint3(const double3& a) {
	return make_uint3(static_cast<uint>(ceil(a.x)),
										static_cast<uint>(ceil(a.y)),
										static_cast<uint>(ceil(a.z)));
}

inline uint3 floor_uint3(const double3& a) {
	return make_uint3(static_cast<uint>(floor(a.x)),
										static_cast<uint>(floor(a.y)),
										static_cast<uint>(floor(a.z)));
}

inline double3 operator*(const uint3& a, double b) {
	return double3(a.x * b, a.y * b, a.z * b);
}

#endif
