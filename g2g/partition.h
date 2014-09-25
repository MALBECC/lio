#ifndef __CUBES_H__
#define __CUBES_H__

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <cstdio>

#include "scalar_vector_types.h"
#include "timer.h"

#include "global_memory_pool.h"
#include "buffer_pool.h"

using std::cout;
using std::endl;
using std::pair;

namespace G2G {
  struct Timers {
    Timer memcpy, trmms, density_calcs, total, ciclos, rmm, density, forces, resto, pot, functions, density_derivs;
  };

  std::ostream& operator<<(std::ostream& io, const Timers& t);

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

enum FunctionType { FUNCTION_S = 1, FUNCTION_P = 3, FUNCTION_D = 6 };

template<class scalar_type>
class PointGroup {
  public:
    PointGroup(void) : number_of_points(0), s_functions(0), p_functions(0), d_functions(0), inGlobal(false) {  }
    virtual ~PointGroup(void);
    std::vector<Point> points;
    uint number_of_points;
    uint s_functions, p_functions, d_functions;

    G2G::HostMatrixUInt func2global_nuc; // size == total_functions_simple()
    G2G::HostMatrixUInt func2local_nuc; // size == total_functions()

    std::vector<uint> local2global_func; // size == total_functions_simple()
    std::vector<uint> local2global_nuc;  // size == total_nucleii()

    typedef vec_type<scalar_type,2> vec_type2;
    typedef vec_type<scalar_type,3> vec_type3;
    typedef vec_type<scalar_type,4> vec_type4;

    #if CPU_KERNELS
    G2G::HostMatrix<scalar_type> function_values;
    G2G::HostMatrix<vec_type3> gradient_values;
    G2G::HostMatrix<vec_type3> hessian_values;
    G2G::HostMatrix<scalar_type> gX, gY, gZ;
    G2G::HostMatrix<scalar_type> hIX, hIY, hIZ;
    G2G::HostMatrix<scalar_type> hPX, hPY, hPZ;
    #else
    G2G::CudaMatrix<scalar_type> function_values;
    G2G::CudaMatrix<vec_type4> gradient_values;
    G2G::CudaMatrix<vec_type4> hessian_values;
    #endif

    long long cost() const;
    inline FunctionType small_function_type(uint f) const {
      if (f < s_functions) return FUNCTION_S;
      else if (f < s_functions + p_functions) return FUNCTION_P;
      else return FUNCTION_D;
    }
    //Las funciones totales, son totales del grupo, no las totales para todos los grupos.
    inline uint total_functions(void) const { return s_functions + p_functions * 3 + d_functions * 6; }
    inline uint total_functions_simple(void) const { return local2global_func.size(); } // == s_functions + p_functions + d_functions
    inline uint total_nucleii(void) const { return local2global_nuc.size(); }
    inline bool has_nucleii(uint atom) const { return (std::find(local2global_nuc.begin(), local2global_nuc.end(), atom) != local2global_nuc.end()); }

    void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input) const;
    void add_rmm_output(const G2G::HostMatrix<scalar_type>& rmm_output) const;

    void compute_nucleii_maps(void);

    void add_point(const Point& p);
    void compute_weights(void);

    void compute_functions(bool forces, bool gga);
    void solve(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, 
        bool compute_energy, double& energy, HostMatrix<double> &, ThreadBufferPool<scalar_type> &, int, HostMatrix<scalar_type> &);

    bool is_significative(FunctionType, double exponent, double coeff, double d2);
    bool operator<(const PointGroup<scalar_type>& T) const;
    size_t size_in_gpu() const;
    int pool_elements() const;

    bool inGlobal;

};

#if FULL_DOUBLE
typedef double base_scalar_type;
#else
typedef float base_scalar_type;
#endif

class Sphere : public PointGroup<base_scalar_type> {
  public:
    Sphere(void);
    Sphere(uint _atom, double _radius);

    void assign_significative_functions(const std::vector<double>& min_exps, const std::vector<double>& min_coeff);

    uint atom;
    double radius;
};

class Cube : public PointGroup<base_scalar_type> {
  public:
    void assign_significative_functions(const double3& cube_coord, const std::vector<double>& min_exps, const std::vector<double>& min_coeff);
};

class Partition {
  public:
    void clear(void) {
      cubes.clear(); spheres.clear();
      work.clear(); pool_sizes.clear();
    }

    void solve(Timers& timers, bool compute_rmm,bool lda,bool compute_forces, bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr);

    void regenerate(void);

    void compute_functions(bool forces, bool gga)
    {
      Timer t1;
      t1.start_and_sync();

      #pragma omp parallel for
      for(int i = 0; i < cubes.size(); i++){
        cubes[i].compute_functions(forces, gga);
      }

      #pragma omp parallel for
      for(int i = 0; i < spheres.size(); i++){
        spheres[i].compute_functions(forces, gga);
      }

      t1.stop_and_sync();
    }

    std::vector<Cube> cubes;
    std::vector<Sphere> spheres;

    std::vector< std::vector< int > > work;
    std::vector< int > pool_sizes;
    int outer_threads, inner_threads;
};

}

#define ALIGN(x) (((x) % 16) ? ((x) + 16 - ((x) % 16)) : (x))

#endif
