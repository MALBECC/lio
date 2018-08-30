#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <climits>
#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "timer.h"
using namespace std;

namespace G2G {

int MINCOST, THRESHOLD, SPLITPOINTS;
Partition partition;

ostream& operator<<(ostream& io, const Timers& t) {
  ostringstream ss;
  ss << "density = " << t.density << " rmm = " << t.rmm
     << " forces = " << t.forces << " functions = " << t.functions;
  io << ss.str() << endl;
  return io;
}

/********************
 * PointGroup
 ********************/

template <class scalar_type>
void PointGroupCPU<scalar_type>::output_cost() const {
  // printf("%d %d %d %lld %lld\n", number_of_points, total_functions(),
  // rmm_bigs.size(), cost(), size_in_gpu());
}

template <class scalar_type>
bool PointGroupCPU<scalar_type>::is_big_group() const {
  return false;
}

template <class scalar_type>
bool PointGroupGPU<scalar_type>::is_big_group() const {
  return true;
}

template <class scalar_type>
void PointGroupGPU<scalar_type>::get_rmm_input(
    HostMatrix<scalar_type>& rmm_input, FortranMatrix<double>& source) const {
  rmm_input.zero();
  for (uint i = 0, ii = 0; i < this->total_functions_simple(); i++) {
    uint inc_i = this->small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = this->local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < this->total_functions_simple(); j++) {
        uint inc_j = this->small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = this->local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index =
              (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) +
              (big_j - big_i);

          rmm_input(ii, jj) = (scalar_type)source.data[big_index];

          rmm_input(jj, ii) = rmm_input(ii, jj);
        }
      }
    }
  }
}

template <class scalar_type>
void PointGroupCPU<scalar_type>::get_rmm_input(
    HostMatrix<scalar_type>& rmm_input, FortranMatrix<double>& source) const {
  rmm_input.zero();
  const int indexes = this->rmm_bigs.size();
  for (int i = 0; i < indexes; i++) {
    int ii = this->rmm_rows[i], jj = this->rmm_cols[i], bi = this->rmm_bigs[i];
    rmm_input(ii, jj) = rmm_input(jj, ii) = (scalar_type)source.data[bi];
  }
}

template <class scalar_type>
void PointGroupCPU<scalar_type>::get_rmm_input(
    HostMatrix<scalar_type>& rmm_input) const {
  get_rmm_input(rmm_input, fortran_vars.rmm_input_ndens1);
}

template <class scalar_type>
void PointGroupGPU<scalar_type>::get_rmm_input(
    HostMatrix<scalar_type>& rmm_input) const {
  get_rmm_input(rmm_input, fortran_vars.rmm_input_ndens1);
}

template <class scalar_type>
void PointGroupCPU<scalar_type>::get_rmm_input(
    HostMatrix<scalar_type>& rmm_input_a,
    HostMatrix<scalar_type>& rmm_input_b) const {
  get_rmm_input(rmm_input_a, fortran_vars.rmm_dens_a);
  get_rmm_input(rmm_input_b, fortran_vars.rmm_dens_b);
}

template <class scalar_type>
void PointGroupGPU<scalar_type>::get_rmm_input(
    HostMatrix<scalar_type>& rmm_input_a,
    HostMatrix<scalar_type>& rmm_input_b) const {
  get_rmm_input(rmm_input_a, fortran_vars.rmm_dens_a);
  get_rmm_input(rmm_input_b, fortran_vars.rmm_dens_b);
}

template <class scalar_type>
void PointGroup<scalar_type>::compute_indexes() {
  rmm_bigs.clear();
  rmm_cols.clear();
  rmm_rows.clear();
  for (uint i = 0, ii = 0; i < this->total_functions_simple(); i++) {
    uint inc_i = this->small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = this->local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < this->total_functions_simple(); j++) {
        uint inc_j = this->small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = this->local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index =
              (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) +
              (big_j - big_i);
          if (ii > jj) swap(ii, jj);
          rmm_rows.push_back(ii);
          rmm_cols.push_back(jj);
          rmm_bigs.push_back(big_index);
        }
      }
    }
  }
}

template <class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(
    const HostMatrix<scalar_type>& rmm_output,
    FortranMatrix<double>& target) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index =
              (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) +
              (big_j - big_i);
          target(big_index) += (double)rmm_output(ii, jj);
        }
      }
    }
  }
}

template <class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(
    const HostMatrix<scalar_type>& rmm_output,
    HostMatrix<double>& target) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index =
              (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) +
              (big_j - big_i);
          target(big_index) += (double)rmm_output(ii, jj);
        }
      }
    }
  }
}

template <class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(
    const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output);
}

template <class scalar_type>
void PointGroup<scalar_type>::add_rmm_output_a(
    const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output_a);
}

template <class scalar_type>
void PointGroup<scalar_type>::add_rmm_output_b(
    const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output_b);
}

template <class scalar_type>
void PointGroup<scalar_type>::add_rmm_open_output(
    const HostMatrix<scalar_type>& rmm_output_a,
    const HostMatrix<scalar_type>& rmm_output_b) const {
  add_rmm_output(rmm_output_a, fortran_vars.rmm_output_a);
  add_rmm_output(rmm_output_b, fortran_vars.rmm_output_b);
}

template <class scalar_type>
void PointGroup<scalar_type>::add_point(const Point& p) {
  points.push_back(p);
  number_of_points++;
}

#define EXP_PREFACTOR 1.01057089636005  // (2 * pow(4, 1/3.0)) / M_PI

template <class scalar_type>
bool PointGroup<scalar_type>::is_significative(FunctionType type,
                                               double exponent, double coeff,
                                               double d2) {
  switch (type) {
    case FUNCTION_S:
      return (exponent * d2 <
              max_function_exponent - log(pow((2. * exponent / M_PI), 3)) / 4);
      break;
    default: {
      double x = 1;
      double delta;
      double e = 0.1;
      double factor = pow((2.0 * exponent / M_PI), 3);
      factor = sqrt(factor * 4.0 * exponent);
      double norm = (type == FUNCTION_P ? sqrt(factor) : abs(factor));
      do {
        double div = (type == FUNCTION_P ? log(x) : 2 * log(x));
        double x1 = sqrt((max_function_exponent - log(norm) + div) / exponent);
        delta = abs(x - x1);
        x = x1;
      } while (delta > e);
      return (sqrt(d2) < x);
    } break;
  }
}

template <class scalar_type>
long long PointGroup<scalar_type>::cost() const {
  long long np = number_of_points, gm = total_functions();
  return 10 * ((np * gm * (1 + gm)) / 2) + MINCOST;
}

template <class scalar_type>
bool PointGroup<scalar_type>::operator<(
    const PointGroup<scalar_type>& T) const {
  return cost() < T.cost();
}

template <class scalar_type>
int PointGroup<scalar_type>::elements() const {
  int t = total_functions(), n = number_of_points;
  return t * n;
}

template <class scalar_type>
size_t PointGroup<scalar_type>::size_in_gpu() const {
  uint total_cost = 0;
  uint single_matrix_cost =
      COALESCED_DIMENSION(number_of_points) * total_functions();

  total_cost += single_matrix_cost;  // 1 scalar_type functions
  if (fortran_vars.do_forces || fortran_vars.gga)
    total_cost += (single_matrix_cost * 4);  // 4 vec_type gradient
  if (fortran_vars.gga)
    total_cost += (single_matrix_cost * 8);  // 2*4 vec_type hessian
  return total_cost *
         sizeof(scalar_type);  // size in bytes according to precision
}

template<class scalar_type>
PointGroup<scalar_type>::~PointGroup<scalar_type>() {
}

template <class scalar_type>
PointGroupCPU<scalar_type>::~PointGroupCPU<scalar_type>() {
  deallocate();
}

template <class scalar_type>
void PointGroupCPU<scalar_type>::deallocate() {
  function_values.deallocate();
  gX.deallocate();
  gY.deallocate();
  gZ.deallocate();
  hPX.deallocate();
  hPY.deallocate();
  hPZ.deallocate();
  hIX.deallocate();
  hIY.deallocate();
  hIZ.deallocate();
  function_values_transposed.deallocate();
}

template <class scalar_type>
void PointGroupGPU<scalar_type>::deallocate() {
  if (this->inGlobal) {
    GlobalMemoryPool::dealloc(this->size_in_gpu(), current_device);
    function_values.deallocate();
    gradient_values.deallocate();
    hessian_values_transposed.deallocate();
    this->inGlobal = false;
  }
}

template <class scalar_type>
PointGroupGPU<scalar_type>::~PointGroupGPU<scalar_type>() {
  deallocate();
}

void Partition::compute_functions(bool forces, bool gga) {
  Timer t1;
  t1.start_and_sync();

#pragma omp parallel for schedule(guided, 8)
  for (uint i = 0; i < cubes.size(); i++) {
    if (!cubes[i]->is_big_group()) cubes[i]->compute_functions(forces, gga);
  }

#pragma omp parallel for schedule(guided, 8)
  for (uint i = 0; i < spheres.size(); i++) {
    if (!spheres[i]->is_big_group()) spheres[i]->compute_functions(forces, gga);
  }

  t1.stop_and_sync();
  if (timer_single) cout << "Functions: " << t1 << endl;
}

void Partition::clear() {
  for (uint i = 0; i < cubes.size(); i++) delete cubes[i];
  for (uint i = 0; i < spheres.size(); i++) delete spheres[i];
  cubes.clear();
  spheres.clear();
  work.clear();
}

void Partition::rebalance(vector<double>& times, vector<double>& finishes) {
  for (int device = 0; device < 2; device++) {
    for (int rondas = 0; rondas < 5; rondas++) {
      int largest = 0;
      int smallest = 0;

      if (device == 0) {
        largest =
            std::max_element(finishes.begin(), finishes.end() - gpu_threads) -
            finishes.begin();
        smallest =
            std::min_element(finishes.begin(), finishes.end() - gpu_threads) -
            finishes.begin();
      } else {
        largest =
            std::max_element(finishes.begin() + cpu_threads, finishes.end()) -
            finishes.begin();
        smallest =
            std::min_element(finishes.begin() + cpu_threads, finishes.end()) -
            finishes.begin();
      }

      double diff = finishes[largest] - finishes[smallest];

      if (largest != smallest && work[largest].size() > 1) {
        double lt = finishes[largest];
        double moved = 0;
        while (diff / lt >= 0.02) {
          int mini = -1;
          double currentmini = diff;
          for (uint i = 0; i < work[largest].size(); i++) {
            int ind = work[largest][i];
            if (times[ind] > diff / 2) continue;
            double cost = times[ind];
            if (currentmini > diff - 2 * cost) {
              currentmini = diff - 2 * cost;
              mini = i;
            }
          }

          if (mini == -1) {
            // printf("Nothing more to swap!\n");
            break;
          }

          int topass = mini;
          int workindex = work[largest][topass];

          // printf("Swapping %d from %d to %d\n", work[largest][topass],
          // largest, smallest);

          if (device == 1) {
            if (workindex < (int)cubes.size())
              cubes[workindex]->deallocate();
            else
              spheres[workindex - cubes.size()]->deallocate();
          }

          work[smallest].push_back(work[largest][topass]);
          work[largest].erase(work[largest].begin() + topass);

          diff -= 2 * times[workindex];
          moved += times[workindex];
        }
        finishes[largest] -= moved;
        finishes[smallest] += moved;
      }
    }
  }
}

void Partition::solve(Timers& timers, bool compute_rmm, bool lda,
                      bool compute_forces, bool compute_energy,
                      double* fort_energy_ptr, double* fort_forces_ptr,
                      bool OPEN) {
  double energy = 0.0;

  double cubes_energy = 0, spheres_energy = 0;
  double cubes_energy_i = 0, spheres_energy_i = 0;
  double cubes_energy_c = 0, spheres_energy_c = 0;
  double cubes_energy_c1 = 0, spheres_energy_c1 = 0;
  double cubes_energy_c2 = 0, spheres_energy_c2 = 0;

  Timer smallgroups, biggroups;

// Verificar si anda reduction (+:energy) FF
#pragma omp parallel for num_threads(cpu_threads + gpu_threads) schedule( \
    static) reduction(+ : energy)
  for (uint i = 0; i < work.size(); i++) {
#if GPU_KERNELS
    bool gpu_thread = false;
    if (i >= cpu_threads) {
      gpu_thread = true;
      cudaSetDevice(i - cpu_threads);
    }
#endif
    double local_energy = 0;

    Timers ts;
    Timer t;
    t.start_and_sync();

    if (compute_forces) fort_forces_ms[i].zero();
    if (compute_rmm) {
      if (OPEN) {
        rmm_outputs_a[i].zero();
        rmm_outputs_b[i].zero();
      } else {
        rmm_outputs[i].zero();
      }
    }

    for (uint j = 0; j < work[i].size(); j++) {
      int ind = work[i][j];
      Timer element;
      element.start_and_sync();
      if (OPEN) {
        if (ind >= cubes.size()) {
          spheres[ind - cubes.size()]->solve_opened(
              ts, compute_rmm, lda, compute_forces, compute_energy,
              local_energy, spheres_energy_i, spheres_energy_c,
              spheres_energy_c1, spheres_energy_c2, fort_forces_ms[i],
              rmm_outputs_a[i], rmm_outputs_b[i]);
        } else {
          cubes[ind]->solve_opened(
              ts, compute_rmm, lda, compute_forces, compute_energy,
              local_energy, spheres_energy_i, spheres_energy_c,
              spheres_energy_c1, spheres_energy_c2, fort_forces_ms[i],
              rmm_outputs_a[i], rmm_outputs_b[i]);
        }
      } else {
        if (ind >= cubes.size()) {
          spheres[ind - cubes.size()]->solve_closed(
              ts, compute_rmm, lda, compute_forces, compute_energy,
              local_energy, fort_forces_ms[i], 1, rmm_outputs[i]);
        } else {
          cubes[ind]->solve_closed(ts, compute_rmm, lda, compute_forces,
                                   compute_energy, local_energy,
                                   fort_forces_ms[i], 1, rmm_outputs[i]);
        }
      }

#if GPU_KERNELS
      if (gpu_thread) {
        cudaDeviceSynchronize();
      }
#endif

      element.stop_and_sync();
      timeforgroup[ind] = element.getTotal();
    }
    t.stop_and_sync();

    next[i] = t.getTotal();

    energy += local_energy;
  }

  Timer enditer;
  enditer.start();
  if (work.size() > 1) rebalance(timeforgroup, next);
  if (compute_forces) {
    FortranMatrix<double> fort_forces_out(fort_forces_ptr, fortran_vars.atoms,
                                          3, fortran_vars.max_atoms);
    for (uint k = 0; k < fort_forces_ms.size(); k++) {
      for (int i = 0; i < fortran_vars.atoms; i++) {
        for (int j = 0; j < 3; j++) {
          fort_forces_out(i, j) += fort_forces_ms[k](i, j);
        }
      }
    }
  }

  if (compute_rmm) {
    if (fortran_vars.OPEN) {
      double* dst_a = fortran_vars.rmm_output_a.data;
      double* dst_b = fortran_vars.rmm_output_b.data;
      const int elements =
          fortran_vars.rmm_output_a.width * fortran_vars.rmm_output_a.height;
      const int elements_b =
          fortran_vars.rmm_output_b.width * fortran_vars.rmm_output_b.height;

      if (!(rmm_outputs_a.size() == rmm_outputs_b.size())) {
        std::cout << "ERROR in partition.solve: outputs A and B of different "
                     "size. \n";
      }

      if (!(elements == elements_b)) {
        std::cout << "ERROR in partition.solve: different number of elements A "
                     "and B.\n";
      }

      for (uint k = 0; k < rmm_outputs_a.size(); k++) {
        const double* src_a = rmm_outputs_a[k].asArray();
        const double* src_b = rmm_outputs_b[k].asArray();
#if INTEL_COMP
#pragma ivdep
#pragma vector always
#endif
        for (int i = 0; i < elements; i++) {
          dst_a[i] += src_a[i];
          dst_b[i] += src_b[i];
        }
      }
    } else {
      double* dst = fortran_vars.rmm_output.data;
      const int elements =
          fortran_vars.rmm_output.width * fortran_vars.rmm_output.height;
      for (uint k = 0; k < rmm_outputs.size(); k++) {
        const double* src = rmm_outputs[k].asArray();
#if INTEL_COMP
#pragma ivdep
#pragma vector always
#endif
        for (int i = 0; i < elements; i++) {
          dst[i] += src[i];
        }
      }
    }
  }

/*  if (OPEN && compute_energy) {
    std::cout << " Ei:  " << cubes_energy_i + spheres_energy_i << std::endl;
    std::cout << " Ec:  " << cubes_energy_c + spheres_energy_c << std::endl;
    std::cout << " Ec1: " << cubes_energy_c1 + spheres_energy_c1 << std::endl;
    std::cout << " Ec2: " << cubes_energy_c2 + spheres_energy_c2 << std::endl;
  }*/

  *fort_energy_ptr = energy;
  if (*fort_energy_ptr != *fort_energy_ptr) {
    std::cout << "I see dead peaple " << std::endl;
#if GPU_KERNELS
    cudaDeviceReset();
#endif
    exit(1);
  }
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#if GPU_KERNELS
template class PointGroupGPU<double>;
#endif
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#if GPU_KERNELS
template class PointGroupGPU<float>;
#endif
#endif
}
