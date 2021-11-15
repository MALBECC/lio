#ifndef __OS_H__
#define __OS_H__

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include "../scalar_vector_types.h"
#include "../timer.h"

#include "aint_common.h"

namespace AINT {

template <class scalar_type>
class OSIntegral {
 public:
  //
  // On host
  //
  bool apply_cutoff;

  // Details of the distribution of the valid primitive-primitive and density
  // terms
  uint term_type_offsets[MAX_TERM_TYPES];
  uint dens_offsets[MAX_TERM_TYPES];
  uint out_offsets[MAX_TERM_TYPES];
  uint energies_offsets[MAX_TERM_TYPES];
  uint term_type_counts[MAX_TERM_TYPES];
  uint dens_counts[MAX_TERM_TYPES];

  // Input
  std::vector<uint> func_code;
  std::vector<uint> local_dens;
  std::vector<scalar_type> dens_values;

  std::vector<uint> local2globaldens;

#if GPU_KERNELS
  //
  // On device
  //
  // Common input for O-S integral evaluation
  G2G::CudaMatrix<G2G::vec_type<scalar_type, 2> > factor_ac_dev;
  G2G::CudaMatrixUInt nuc_dev;
  G2G::CudaMatrixUInt func_code_dev;
  G2G::CudaMatrixUInt local_dens_dev;
  G2G::CudaMatrix<scalar_type> dens_values_dev;

  // Common output for O-S integral evaluation (energy or gradient)
  G2G::CudaMatrix<double> partial_fock_dev;
  G2G::CudaMatrix<double> partial_energies_dev;
  G2G::CudaMatrix<G2G::vec_type<scalar_type, 3> > partial_qm_forces_dev;

  int my_device;
#endif

  //
  // Functions
  //
  OSIntegral(void) : apply_cutoff(true) {}
  OSIntegral(bool _apply_cutoff) : apply_cutoff(_apply_cutoff) {}

  // Apply cutoff, get per-thread information
  void new_cutoff(void);

  // Allocate common input/output on GPU
  void load_params(void);

  bool load_input(void);
  bool alloc_output(void);

  void reload_density(void);

  // Get common results from GPU and send to final output
  void get_fock_output(double& Es, G2G::FortranMatrix<double>& fock_out);
  void get_gradient_output(double* qm_forces, uint partial_size);

  void clear(void);
  void deinit(void);
};
}  // namespace AINT

#endif
