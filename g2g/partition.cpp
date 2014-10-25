#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "timer.h"
#include "mkl.h"
using namespace std;

namespace G2G {
Partition partition;

ostream& operator<<(ostream& io, const Timers& t) {
#ifdef TIMINGS
    ostringstream ss;
    ss << "memcpys: " << t.memcpy << "trmms: " << t.trmms << "density_calcs: " << t.density_calcs << "rmm: " << t.rmm << " density: " 
       << t.density << " pot: " << t.pot << " forces: " << t.forces << " resto: " << t.resto << " functions: " << t.functions
       << " rmm_input = " << t.rmm_input << " rmm_ssyr = " << t.rmm_calcs << " rmm_update = " << t.rmm_update;
    io << ss.str() << endl;
#endif
  return io;
}

/********************
 * PointGroup
 ********************/

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input,
    FortranMatrix<double>& source) const {
  rmm_input.zero();
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);

          rmm_input(ii, jj) = (scalar_type)source.data[big_index];
          rmm_input(jj, ii) = rmm_input(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input) const {
  get_rmm_input(rmm_input, fortran_vars.rmm_input_ndens1);
}

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input_a, HostMatrix<scalar_type>& rmm_input_b) const {
  get_rmm_input(rmm_input_a, fortran_vars.rmm_dens_a);
  get_rmm_input(rmm_input_b, fortran_vars.rmm_dens_b);
}

template<class scalar_type>
void PointGroup<scalar_type>::compute_indexes()
{
  rmm_bigs.clear(); rmm_rows.clear(); rmm_cols.clear(); 
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          rmm_bigs.push_back(big_index);
          rmm_rows.push_back(ii); rmm_cols.push_back(jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(const HostMatrix<scalar_type>& rmm_output,
    FortranMatrix<double>& target ) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          target(big_index) += (double)rmm_output(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output_a(const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output_a);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output_b(const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output_b);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_open_output(const HostMatrix<scalar_type>& rmm_output_a,
    const HostMatrix<scalar_type>& rmm_output_b) const {
  add_rmm_output(rmm_output_a, fortran_vars.rmm_output_a);
  add_rmm_output(rmm_output_b, fortran_vars.rmm_output_b);
}

template<class scalar_type>
void PointGroup<scalar_type>::compute_nucleii_maps(void)
{
  if (total_functions_simple() != 0) {
    func2global_nuc.resize(total_functions_simple());
    for (uint i = 0; i < total_functions_simple(); i++) {
      func2global_nuc(i) = fortran_vars.nucleii(local2global_func[i]) - 1;
    }

    func2local_nuc.resize(total_functions());
    uint ii = 0;
    for (uint i = 0; i < total_functions_simple(); i++) {
      uint global_atom = func2global_nuc(i);
      uint local_atom = std::distance(local2global_nuc.begin(),
          std::find(local2global_nuc.begin(), local2global_nuc.end(), global_atom));
      uint inc = small_function_type(i);
      for (uint k = 0; k < inc; k++, ii++) func2local_nuc(ii) = local_atom;
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_point(const Point& p) {
  points.push_back(p);
  number_of_points++;
}

#define EXP_PREFACTOR 1.01057089636005 // (2 * pow(4, 1/3.0)) / M_PI

template<class scalar_type>
bool PointGroup<scalar_type>::is_significative(FunctionType type, double exponent, double coeff, double d2) {
  switch(type) {
    case FUNCTION_S:
      return (exponent * d2 < max_function_exponent-log(pow((2.*exponent/M_PI),3))/4);
    break;
    default:
    {
      double x = 1;
      double delta;
      double e = 0.1;
      double factor = pow((2.0*exponent/M_PI),3);
      factor = sqrt(factor*4.0*exponent) ;
      double norm = (type == FUNCTION_P ? sqrt(factor) : abs(factor)) ;
      do {
        double div = (type == FUNCTION_P ? log(x) : 2 * log(x));
        double x1 = sqrt((max_function_exponent - log(norm) + div) / exponent);
        delta = abs(x-x1);
        x = x1;
      } while (delta > e);
      return (sqrt(d2) < x);
    }
    break;
  }
}

template<class scalar_type>
long long PointGroup<scalar_type>::cost() const {
  long long MINCOST = 60000;
  long long np = number_of_points, gm = total_functions();
  return 10*((np * gm * (1+gm)) / 2) + 2 * (gm * gm * np);
}

template<class scalar_type>
bool PointGroup<scalar_type>::operator<(const PointGroup<scalar_type>& T) const{
    return cost() < T.cost();
}

template<class scalar_type>
int PointGroup<scalar_type>::pool_elements() const {
    int t = total_functions(), n = number_of_points;
    return t * n;
}

template<class scalar_type>
size_t PointGroup<scalar_type>::size_in_gpu() const
{
  uint total_cost=0;
  uint single_matrix_cost = COALESCED_DIMENSION(number_of_points) * total_functions();

  total_cost += single_matrix_cost;       //1 scalar_type functions
  if (fortran_vars.do_forces || fortran_vars.gga)
    total_cost += (single_matrix_cost*4); //4 vec_type gradient
  if (fortran_vars.gga)
    total_cost+= (single_matrix_cost*8);  //2*4 vec_type hessian
  return total_cost*sizeof(scalar_type);  // size in bytes according to precision
}

template<class scalar_type>
PointGroup<scalar_type>::~PointGroup<scalar_type>()
{
#if !CPU_KERNELS
  if(inGlobal) {
    globalMemoryPool::dealloc(size_in_gpu());
    function_values.deallocate();
    gradient_values.deallocate();
    hessian_values.deallocate();
  }
#else
  function_values.deallocate();
  gradient_values.deallocate();
  hessian_values.deallocate();
  gX.deallocate(); gY.deallocate(); gZ.deallocate();
  hPX.deallocate(); hPY.deallocate(); hPZ.deallocate();
  hIX.deallocate(); hIY.deallocate(); hIZ.deallocate();
  function_values_transposed.deallocate();
#endif
}

void Partition::compute_functions(bool forces, bool gga) { 
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

void Partition::clear() {
  cubes.clear(); spheres.clear(); work.clear(); 
}

void Partition::solve(Timers& timers, bool compute_rmm,bool lda,bool compute_forces, 
                      bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr, bool OPEN){
  double energy = 0.0;

  double cubes_energy = 0, spheres_energy = 0;
  double cubes_energy_i = 0, spheres_energy_i = 0;
  double cubes_energy_c = 0, spheres_energy_c = 0;
  double cubes_energy_c1 = 0, spheres_energy_c1 = 0;
  double cubes_energy_c2 = 0, spheres_energy_c2 = 0;

  HostMatrix<double> fort_forces_ms[outer_threads];
  HostMatrix<base_scalar_type> rmm_outputs[outer_threads];

  int order[outer_threads]; int next = 0;

  #pragma omp parallel for reduction(+:energy) num_threads(outer_threads)
  for(int i = 0; i< work.size(); i++) {
    double local_energy = 0; Timers ts; Timer t;
    int id = omp_get_thread_num();

    t.start();
    long long cost = 0;

    if(compute_rmm) rmm_outputs[i].resize(fortran_vars.rmm_output.width, fortran_vars.rmm_output.height);
    if(compute_forces) fort_forces_ms[i].resize(fortran_vars.max_atoms, 3);

    for(int j = 0; j < work[i].size(); j++) {
      int ind = work[i][j];
      if (ind >= cubes.size()) {
        spheres[ind-cubes.size()].solve(ts, compute_rmm,lda,compute_forces, compute_energy, 
          local_energy, spheres_energy_i, spheres_energy_c, spheres_energy_c1, spheres_energy_c2,
          fort_forces_ms[i], inner_threads, rmm_outputs[i], OPEN);
        cost += spheres[ind-cubes.size()].cost();
      } else {
        cubes[ind].solve(ts, compute_rmm,lda,compute_forces, compute_energy, 
          local_energy, cubes_energy_i, cubes_energy_c, cubes_energy_c1, cubes_energy_c2,
          fort_forces_ms[i], inner_threads, rmm_outputs[i], OPEN);
        cost += cubes[ind].cost();
      }
    }

    t.stop();
    printf("Workload %d took %lus %lums and it has %lu elements (%lld nanounits) (%d)\n", i, 
      t.getSec(), t.getMicrosec(), work[i].size(), cost, id);
    cout << "breakdown: " << ts << endl;

    energy += local_energy;

    order[i] = next;
    #pragma omp atomic
    next++;
  }

  // Work steal
  int first = 0, last = 0;
  for(int i = 0; i < outer_threads; i++) {
    if (order[i] == 0) first = i;
    if (order[i] == outer_threads-1) last = i;
  }

  if(first != last && work[last].size() >= 1) {
    int group = work[last].back();
    work[first].push_back(group);
    work[last].pop_back();
  }

  if (compute_forces) {
    FortranMatrix<double> fort_forces_out(fort_forces_ptr, fortran_vars.atoms, 3, fortran_vars.max_atoms);
    for(int k = 0; k < outer_threads; k++) {
      for(int i = 0; i < fortran_vars.atoms; i++) {
        for(int j = 0; j < 3; j++) {
          fort_forces_out(i,j) += fort_forces_ms[k](i,j);
        }
      }
    }
  }

  if (compute_rmm) {
    for(int k = 0; k < outer_threads; k++) {
      for(int i = 0; i < rmm_outputs[k].width; i++) {
        for(int j = 0; j < rmm_outputs[k].height; j++) {
          fortran_vars.rmm_output(i,j) += rmm_outputs[k](i,j);
        }
      }
    }
  }
    
  if(OPEN && compute_energy) {
    std::cout << "Ei: " << cubes_energy_i+spheres_energy_i;
    std::cout << " Ec: " << cubes_energy_c+spheres_energy_c;
    std::cout << " Ec1: " << cubes_energy_c1+spheres_energy_c1;
    std::cout << " Ec2: " << cubes_energy_c2+spheres_energy_c2 << std::endl;
  }

  *fort_energy_ptr = energy;
  if(*fort_energy_ptr != *fort_energy_ptr) {
      std::cout << "I see dead peaple " << std::endl;
#ifndef CPU_KERNELS
    cudaDeviceReset();
#endif
     exit(1);
   }
}

/**********************
 * Sphere
 **********************/
Sphere::Sphere(void) : atom(0), radius(0) { }
Sphere::Sphere(uint _atom, double _radius) : atom(_atom), radius(_radius) { }

/**********************
 * Cube
 **********************/

template class PointGroup<double>;
template class PointGroup<float>;
}
