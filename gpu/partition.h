#ifndef __CUBES_H__
#define __CUBES_H__

#include <list>
#include <set>
#include <vector>
#include "cuda/double.h"

/********************
 * Point information
 ********************/
struct Point {
	Point(uint _atom, uint _shell, uint _point, double3 _position, double _weight) :
		atom(_atom), shell(_shell), point(_point), position(_position), weight(_weight) {}
	
	uint atom, shell, point;
	double3 position;
	double weight;
};

class PointGroup {
  public:
    PointGroup(void) : number_of_points(0), s_functions(0), p_functions(0), d_functions(0) {}
    std::list<Point> points;
    uint number_of_points;

  	uint s_functions, p_functions, d_functions;
  	std::set<uint> functions;
  	std::set<uint> nucleii;
  	#ifdef COMPUTE_FUNCTIONS_CPU
    G2G::HostMatrixFloat function_values;
    G2G::HostMatrixFloat4 gradient_values;
    #else
    G2G::CudaMatrixFloat function_values;
  	G2G::CudaMatrixFloat4 gradient_values;
    #endif

    void add_point(const Point& p);
    uint total_functions(void);
    void compute_weights(void);

    bool is_sphere; // for debugging
    bool is_cube;
};

class Sphere : public PointGroup {
  public:
    Sphere(void);
    Sphere(uint _atom, double _radius);

    void assign_significative_functions(const std::vector<double>& min_exps);

    uint atom;
    double radius;
};

class Cube : public PointGroup {
  public:
		Cube(void);
    void assign_significative_functions(const double3& cube_coord, const std::vector<double>& min_exps);
};

extern std::list<PointGroup> final_partition;
void regenerate_partition(void);

#endif
