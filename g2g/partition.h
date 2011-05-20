#ifndef __CUBES_H__
#define __CUBES_H__

#include <list>
#include <set>
#include <vector>
#include <algorithm>

namespace G2G {
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
    PointGroup(void) : number_of_points(0), s_functions(0), p_functions(0), d_functions(0) {  }
    virtual ~PointGroup(void) { }
    std::list<Point> points;
    uint number_of_points;
    uint s_functions, p_functions, d_functions;

    G2G::HostMatrixUInt func2global_nuc; // size == total_functions_simple()
    G2G::HostMatrixUInt func2local_nuc; // size == total_functions()

    std::vector<uint> local2global_func; // size == total_functions_simple()
    std::vector<uint> local2global_nuc;  // size == total_nucleii()
    #if CPU_KERNELS
    G2G::HostMatrix<real> function_values;
    G2G::HostMatrix<creal3> gradient_values;
    G2G::HostMatrix<creal3> hessian_values;
    #else
    G2G::CudaMatrixFloat function_values;
    G2G::CudaMatrixFloat4 gradient_values;
    G2G::CudaMatrixFloat4 hessian_values;
    #endif

    void add_point(const Point& p);
    void compute_weights(void);
    
    template<bool forces, bool gga> void compute_functions(void);

    enum FunctionType { FUNCTION_S = 1, FUNCTION_P = 3, FUNCTION_D = 6 };
    inline FunctionType small_function_type(uint f) const {
      if (f < s_functions) return FUNCTION_S;
      else if (f < s_functions + p_functions) return FUNCTION_P;
      else return FUNCTION_D;
    }

    inline uint total_functions(void) const { return s_functions + p_functions * 3 + d_functions * 6; }
    inline uint total_functions_simple(void) const { return local2global_func.size(); } // == s_functions + p_functions + d_functions
    inline uint total_nucleii(void) const { return local2global_nuc.size(); }
    inline bool has_nucleii(uint atom) const { return (std::find(local2global_nuc.begin(), local2global_nuc.end(), atom) != local2global_nuc.end()); }

    void get_rmm_input(G2G::HostMatrix<real>& rmm_input) const;
    void add_rmm_output(const G2G::HostMatrix<real>& rmm_output) const;

    void compute_nucleii_maps(void);
    virtual bool is_sphere(void) = 0;
    virtual bool is_cube(void) = 0;

    bool is_significative(FunctionType, double exponent, double coeff, double d2);
};

class Sphere : public PointGroup {
  public:
    Sphere(void);
    Sphere(uint _atom, double _radius);

    void assign_significative_functions(const std::vector<double>& min_exps, const std::vector<double>& min_coeff);
    bool is_sphere(void) { return true; }
    bool is_cube(void) { return false; }

    uint atom;
    double radius;
};

class Cube : public PointGroup {
  public:
    Cube(void);
    void assign_significative_functions(const double3& cube_coord, const std::vector<double>& min_exps, const std::vector<double>& min_coeff);
    bool is_sphere(void) { return false; }
    bool is_cube(void) { return true; }

};

class Partition {
  public:
    void clear(void) {
      for (std::list<PointGroup*>::iterator it = group_list.begin(); it != group_list.end(); ++it) { delete *it; }
      group_list.clear();
    }

    void regenerate(void);
    template <bool forces, bool gga> void compute_functions(void);
    std::list<PointGroup*> group_list;
};

extern Partition partition;
}

#endif
