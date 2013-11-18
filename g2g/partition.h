#ifndef __CUBES_H__
#define __CUBES_H__

#include <list>
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include "scalar_vector_types.h"
#include "timer.h"

namespace G2G {
  struct Timers {
    Timer total, ciclos, rmm, density, forces, resto, pot, functions, density_derivs;
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
    PointGroup(void) : number_of_points(0), s_functions(0), p_functions(0), d_functions(0) {  }
    virtual ~PointGroup(void) { }
    std::list<Point> points;
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
    #else
    G2G::CudaMatrix<scalar_type> function_values;
    G2G::CudaMatrix<vec_type4> gradient_values;
    G2G::CudaMatrix<vec_type4> hessian_values;
    #endif

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
    void solve(Timers& timers, bool compute_rmm, bool lda, bool compute_forces, bool compute_energy, double& energy, double* fort_forces_ptr);

    bool is_significative(FunctionType, double exponent, double coeff, double d2);
    bool operator<(const PointGroup<scalar_type>& T) const;
    size_t size_in_gpu() const;

    virtual bool is_sphere(void) = 0;
    virtual bool is_cube(void) = 0;    
};

#if FULL_DOUBLE
class Sphere : public PointGroup<double> {
#else
class Sphere : public PointGroup<float> {
#endif
  public:
    Sphere(void);
    Sphere(uint _atom, double _radius);

    void assign_significative_functions(const std::vector<double>& min_exps, const std::vector<double>& min_coeff);
    bool is_sphere(void) { return true; }
    bool is_cube(void) { return false; }

    uint atom;
    double radius;
};

#if FULL_DOUBLE
class Cube : public PointGroup<double> {
#else
class Cube : public PointGroup<float> {
#endif
  public:
    void assign_significative_functions(const double3& cube_coord, const std::vector<double>& min_exps, const std::vector<double>& min_coeff);
    bool is_sphere(void) { return false; }
    bool is_cube(void) { return true; }

};

class Partition {
  public:
    void clear(void) {
      for (std::list<Cube*>::const_iterator it = cubes.begin(); it != cubes.end(); ++it) delete *it;
      for (std::list<Sphere*>::const_iterator it = spheres.begin(); it != spheres.end(); ++it) delete *it;
      cubes.clear(); spheres.clear();
    }

    void solve(Timers& timers, bool compute_rmm,bool lda,bool compute_forces, bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr)
    {
      double cubes_energy = 0, spheres_energy = 0;

      long long int accumulated_size=0;
      for (std::list<Cube*>::const_iterator it = cubes.begin(); it != cubes.end(); ++it)
      {
        (*it)->solve(timers, compute_rmm,lda,compute_forces, compute_energy, cubes_energy, fort_forces_ptr);
        accumulated_size+=(*it)->size_in_gpu();
//        printf("\t\t\t\t So far %luKb\n",accumulated_size/1024);
      }

      for (std::list<Sphere*>::const_iterator it = spheres.begin(); it != spheres.end(); ++it)
        (*it)->solve(timers, compute_rmm,lda,compute_forces, compute_energy, spheres_energy, fort_forces_ptr);
        
//      std::cout << "cubes XC energy: " << cubes_energy << std::endl;
//      std::cout << "spheres XC energy: " << spheres_energy << std::endl;
      *fort_energy_ptr = cubes_energy + spheres_energy;
    }

    void regenerate(void);

    void compute_functions(bool forces, bool gga)
    {
      Timer t1;
      t1.start_and_sync();
      for (std::list<Cube*>::const_iterator it = cubes.begin(); it != cubes.end(); ++it)
        (*it)->compute_functions(forces, gga);
      for (std::list<Sphere*>::const_iterator it = spheres.begin(); it != spheres.end(); ++it)
        (*it)->compute_functions(forces, gga);
      t1.stop_and_sync();
//      std::cout << "TIMER: funcs: " << t1 << std::endl;
    }

    std::list<Cube*> cubes;
    std::list<Sphere*> spheres;
};

extern Partition partition;
}

#endif
