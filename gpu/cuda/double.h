#ifndef __DOUBLE_H__
#define __DOUBLE_H__

#include <cuda_runtime.h>
#include <cmath>

struct double3 {
	double3(void) : x(0.0), y(0.0), z(0.0) {}
	double3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {} 
	double x; double y; double z;
	
	inline double length(void) const {
		return sqrt(x * x + y * y + z * z);
	}
	
	inline double length_square(void) const {
		return x * x + y * y + z * z;
	}
};

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

inline double elem(const double3& a, uint i)
{
	switch(i) {
		case 0: return a.x;
		case 1: return a.y;
		case 2: return a.z;
		default: return a.x;
	}
}

inline double3 operator*(const uint3& a, double b) {
	return double3(a.x * b, a.y * b, a.z * b);
}

#endif
