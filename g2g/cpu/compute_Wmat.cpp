#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cstring>
#include <omp.h>
#include <cstdio>
#include "../common.h"
#include "../init.h"
#include "../cuda_includes.h"
#include "../matrix.h"
#include "../timer.h"
#include "../partition.h"

namespace G2G {
template <class scalar_type>
void PointGroupCPU<scalar_type>::calc_W_mat(HostMatrix<double>& W_output_local){

  const uint group_m = this->total_functions();
  const int npoints = this->points.size();

#if CPU_RECOMPUTE or !GPU_KERNELS
  /** Compute functions **/
  compute_functions(false, true);
#endif

  HostMatrix<scalar_type> factors_cdft;
  factors_cdft.resize(this->points.size(), cdft_vars.regions);

  /** Calculates w factors **/
  for (int point = 0; point < npoints; point++) {
    const scalar_type wp = this->points[point].weight;

    for (int i = 0; i < cdft_vars.regions; i++) {
    for (int j = 0; j < cdft_vars.natom(i); j++) {
      factors_cdft(point,i) = wp * (this->points[point].atom_weights(cdft_vars.atoms(j,i)));
    }
    }
  }

  const int indexes = this->rmm_bigs.size();
  for (int i = 0; i < indexes; i++) {
      int bi = this->rmm_bigs[i],
          row = this->rmm_rows[i],
          col = this->rmm_cols[i];

      double res = 0.0;
      const scalar_type* fvr = function_values_transposed.row(row);
      const scalar_type* fvc = function_values_transposed.row(col);

      for (int point = 0; point < npoints; point++) {
      for (int j = 0; j < cdft_vars.regions; j++) {
        res += fvr[point] * fvc[point] * factors_cdft(point,j) * cdft_vars.Vc(j);
      }
      }
      
      W_output_local(bi) += res;
  }
  
#if CPU_RECOMPUTE or !GPU_KERNELS
  /* clear functions */
  function_values_transposed.deallocate();
#endif
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif

}