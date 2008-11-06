#ifndef __CUBES_H__
#define __CUBES_H__

#include <list>
#include <set>
#include "cuda/double.h"

struct Point {
	Point(uint _atom, uint _shell, uint _point, double3 _position, double _weight) :
		atom(_atom), shell(_shell), point(_point), position(_position), weight(_weight) {}
	
	uint atom, shell, point;
	double3 position;
	double weight;
};

struct LittleCube {
	LittleCube(void) : number_of_points(0), s_functions(0), p_functions(0), d_functions(0) {}
	std::list<Point> points;
	uint number_of_points;
	uint s_functions, p_functions, d_functions;
	std::set<uint> functions;
	std::set<uint> nucleii;
	G2G::CudaMatrixFloat function_values;
};

extern std::list<LittleCube> final_cube;

void regenerate_cubes(void);

const uint max_function_exponent = 8;
const double little_cube_size = 5.0/*2.0*/; // largo de una arista de los cubos chiquitos [Angstrom]
const uint min_points_per_cube = 5;

#endif
