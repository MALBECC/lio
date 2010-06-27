#ifndef __CUBES_H__
#define __CUBES_H__

#include <list>
#include <set>
#include <vector>

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
    std::vector<uint> functions;
    std::set<uint> nucleii;
    std::vector<uint> nuc_map;
    #if CPU_KERNELS
    G2G::HostMatrixFloat function_values;
    G2G::HostMatrixFloat3 gradient_values;
    G2G::HostMatrixFloat3 hessian_values;
    #else
    G2G::CudaMatrixFloat function_values;
    G2G::CudaMatrixFloat4 gradient_values;
    #endif

    void add_point(const Point& p);
    void compute_weights(void);

    enum FunctionType { FUNCTION_S = 1, FUNCTION_P = 3, FUNCTION_D = 6 };
    inline FunctionType small_function_type(uint f) const {
      if (f < s_functions) return FUNCTION_S;
      else if (f < s_functions + p_functions) return FUNCTION_P;
      else return FUNCTION_D;
    }

    inline uint total_functions(void) const { return s_functions + p_functions * 3 + d_functions * 6; }
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
