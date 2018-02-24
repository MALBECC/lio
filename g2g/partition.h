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

using std::cout;
using std::endl;
using std::pair;

namespace G2G {
struct Timers {
  Timer rmm, density, forces, functions, density_derivs;
};

std::ostream& operator<<(std::ostream& io, const Timers& t);

/********************
 * Point information
 ********************/
struct Point {
  Point(uint _atom, uint _shell, uint _point, double3 _position, double _weight)
      : atom(_atom),
        shell(_shell),
        point(_point),
        position(_position),
        weight(_weight) {}

  uint atom, shell, point;
  double3 position;
  double weight;
};

enum FunctionType { FUNCTION_S = 1, FUNCTION_P = 3, FUNCTION_D = 6 };

template <class scalar_type>
class PointGroup {
 public:
  PointGroup(void)
      : number_of_points(0),
        s_functions(0),
        p_functions(0),
        d_functions(0),
        inGlobal(false) {}
  virtual ~PointGroup(void);
  virtual void deallocate() = 0;
  std::vector<Point> points;
  uint number_of_points;
  uint s_functions, p_functions, d_functions;

  G2G::HostMatrixUInt func2global_nuc;  // size == total_functions_simple()
  G2G::HostMatrixUInt func2local_nuc;   // size == total_functions()

  std::vector<uint> local2global_func;  // size == total_functions_simple()
  std::vector<uint> local2global_nuc;   // size == total_nucleii()

  typedef vec_type<scalar_type, 2> vec_type2;
  typedef vec_type<scalar_type, 3> vec_type3;
  typedef vec_type<scalar_type, 4> vec_type4;
  long long cost() const;
  inline FunctionType small_function_type(uint f) const {
    if (f < s_functions)
      return FUNCTION_S;
    else if (f < s_functions + p_functions)
      return FUNCTION_P;
    else
      return FUNCTION_D;
  }
  // Las funciones totales, son totales del grupo, no las totales para todos los
  // grupos.
  inline uint total_functions(void) const {
    int v = s_functions + p_functions * 3 + d_functions * 6;
    return v;
  }
  inline uint total_functions_simple(void) const {
    return local2global_func.size();
  }  // == s_functions + p_functions + d_functions
  inline uint total_nucleii(void) const { return local2global_nuc.size(); }
  inline bool has_nucleii(uint atom) const {
    return (std::find(local2global_nuc.begin(), local2global_nuc.end(), atom) !=
            local2global_nuc.end());
  }

  virtual void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input,
                             FortranMatrix<double>& source) const = 0;
  virtual void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input) const = 0;
  virtual void get_rmm_input(
      G2G::HostMatrix<scalar_type>& rmm_input_a,
      G2G::HostMatrix<scalar_type>& rmm_input_b) const = 0;

  void add_rmm_output(const G2G::HostMatrix<scalar_type>& rmm_output,
                      FortranMatrix<double>& target) const;
  void add_rmm_output(const G2G::HostMatrix<scalar_type>& rmm_output) const;
  void add_rmm_output(const G2G::HostMatrix<scalar_type>& rmm_output,
                      G2G::HostMatrix<double>& rmm_destination) const;

  void add_rmm_output_a(const G2G::HostMatrix<scalar_type>& rmm_output) const;
  void add_rmm_output_b(const G2G::HostMatrix<scalar_type>& rmm_output) const;
  void add_rmm_open_output(
      const G2G::HostMatrix<scalar_type>& rmm_output_a,
      const G2G::HostMatrix<scalar_type>& rmm_output_b) const;

  void compute_nucleii_maps(void);

  void add_point(const Point& p);
  virtual void compute_weights(void) = 0;

  virtual bool is_big_group() const = 0;
  void compute_indexes();
  std::vector<uint> rmm_rows;
  std::vector<uint> rmm_cols;
  std::vector<uint> rmm_bigs;
  virtual void compute_functions(bool forces, bool gga) = 0;

  virtual void solve_opened(Timers& timers, bool compute_rmm, bool lda,
                            bool compute_forces, bool compute_energy,
                            double& energy, double&, double&, double&, double&,
                            HostMatrix<double>&, HostMatrix<double>&,
                            HostMatrix<double>&) = 0;

  virtual void solve_closed(Timers& timers, bool compute_rmm, bool lda,
                            bool compute_forces, bool compute_energy,
                            double& energy, HostMatrix<double>&, int,
                            HostMatrix<double>&) = 0;

  virtual void solve(Timers& timers, bool compute_rmm, bool lda,
                     bool compute_forces, bool compute_energy, double& energy,
                     double&, double&, double&, double&, HostMatrix<double>&,
                     int, HostMatrix<double>&, bool) = 0;

  bool is_significative(FunctionType, double exponent, double coeff, double d2);

  void assign_functions_as_sphere(uint, double, const std::vector<double>&,
                                  const std::vector<double>&);
  void assign_functions_as_cube(const double3&, const std::vector<double>&,
                                const std::vector<double>&);
  void assign_functions(HostMatrix<double>, const std::vector<double>&,
                        const std::vector<double>&);

  bool operator<(const PointGroup<scalar_type>& T) const;
  size_t size_in_gpu() const;
  int elements() const;

  bool inGlobal;
};

template <class scalar_type>
class PointGroupCPU : public PointGroup<scalar_type> {
 public:
  virtual ~PointGroupCPU(void);
  virtual void deallocate();
  virtual void compute_functions(bool, bool);
  virtual void compute_weights(void);
  void output_cost() const;
  bool is_big_group() const;
  virtual void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input,
                             FortranMatrix<double>& source) const;
  virtual void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input) const;
  virtual void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input_a,
                             G2G::HostMatrix<scalar_type>& rmm_input_b) const;
  virtual void solve_opened(Timers& timers, bool compute_rmm, bool lda,
                            bool compute_forces, bool compute_energy,
                            double& energy, double&, double&, double&, double&,
                            HostMatrix<double>&, HostMatrix<double>&,
                            HostMatrix<double>&);

  virtual void solve_closed(Timers& timers, bool compute_rmm, bool lda,
                            bool compute_forces, bool compute_energy,
                            double& energy, HostMatrix<double>&, int,
                            HostMatrix<double>&);

  virtual void solve(Timers& timers, bool compute_rmm, bool lda,
                     bool compute_forces, bool compute_energy, double& energy,
                     double&, double&, double&, double&, HostMatrix<double>&,
                     int, HostMatrix<double>&, bool);

  typedef vec_type<scalar_type, 2> vec_type2;
  typedef vec_type<scalar_type, 3> vec_type3;
  typedef vec_type<scalar_type, 4> vec_type4;
  G2G::HostMatrix<scalar_type> function_values;
  G2G::HostMatrix<scalar_type> gX, gY, gZ;
  G2G::HostMatrix<scalar_type> hIX, hIY, hIZ;
  G2G::HostMatrix<scalar_type> hPX, hPY, hPZ;
  G2G::HostMatrix<scalar_type> function_values_transposed;
};

template<class scalar_type>
class PointGroupGPU: public PointGroup<scalar_type> {
  public:
    virtual ~PointGroupGPU(void);
    virtual void deallocate();
    virtual void compute_functions(bool, bool);
    virtual void compute_weights(void);
    bool is_big_group() const;
    virtual void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input, FortranMatrix<double>& source) const;
    virtual void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input) const;
    virtual void get_rmm_input(G2G::HostMatrix<scalar_type>& rmm_input_a, G2G::HostMatrix<scalar_type>& rmm_input_b) const;
    virtual void solve_opened(Timers& timers, bool compute_rmm, bool lda, bool compute_forces,
        bool compute_energy, double& energy, double &, double &, double &, double &,
        HostMatrix<double> &, HostMatrix<double> &, HostMatrix<double> &);
    virtual void solve_closed(Timers& timers, bool compute_rmm, bool lda, bool compute_forces,
        bool compute_energy, double& energy, HostMatrix<double> &, int, HostMatrix<double> &);

    virtual void solve(Timers& timers, bool compute_rmm, bool lda, bool compute_forces,
        bool compute_energy, double& energy, double &, double &, double &, double &,
        HostMatrix<double> &, int, HostMatrix<double> &, bool);

    typedef vec_type<scalar_type,2> vec_type2;
    typedef vec_type<scalar_type,3> vec_type3;
    typedef vec_type<scalar_type,4> vec_type4;
    G2G::CudaMatrix<scalar_type> function_values;
    G2G::CudaMatrix<vec_type4> gradient_values;
    G2G::CudaMatrix<vec_type4> hessian_values_transposed;
    int current_device;

};

#if FULL_DOUBLE
typedef double base_scalar_type;
#else
typedef float base_scalar_type;
#endif

// =======Partition Class ========//

class Partition {
  public:
    void compute_work_partition();
    void clear(void);
    void regenerate(void);

    void solve(Timers& timers, bool compute_rmm,bool lda,bool compute_forces, bool compute_energy,
               double* fort_energy_ptr, double* fort_forces_ptr, bool OPEN);
    void compute_functions(bool forces, bool gga);
    void rebalance(std::vector<double> &, std::vector<double> &);

    std::vector<PointGroup<base_scalar_type>*> cubes;
    std::vector<PointGroup<base_scalar_type>*> spheres;

    std::vector< HostMatrix<double> > fort_forces_ms;
    std::vector< HostMatrix<double> > rmm_outputs;
    std::vector< HostMatrix<double> > rmm_outputs_a;
    std::vector< HostMatrix<double> > rmm_outputs_b;

    std::vector< std::vector< int > > work;
    std::vector< double > next;
    std::vector< double > timeforgroup;

};

extern int MINCOST, THRESHOLD, SPLITPOINTS;
extern int cpu_threads, gpu_threads;
}

#endif
