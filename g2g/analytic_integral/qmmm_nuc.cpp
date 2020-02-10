#include <fstream>
#include <iostream>
#include <vector>

#include "../common.h"
#include "../init.h"
#include "../matrix.h"
#include "../timer.h"
#include "../scalar_vector_types.h"

#include "aint_init.h"
#include "qmmm_integral.h"

using std::cout;
using std::vector;
using std::endl;

namespace AINT {

template <class scalar_type>
void QMMMIntegral<scalar_type>::calc_nuc_gradient(double* qm_forces,
                                                  double* mm_forces) {
  for (uint i = 0; i < G2G::fortran_vars.atoms; i++) {
    double3 qm_pos = G2G::fortran_vars.atom_positions(i);
    for (uint j = 0; j < integral_vars.clatoms; j++) {
      double3 mm_pos = integral_vars.clatom_positions(j);
      double3 diff = qm_pos - mm_pos;
      double dist = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
      dist = sqrt(dist);

      double prefactor = -integral_vars.clatom_charges(j) *
                          integral_vars.atom_Z(i) / pow(dist, 3.0);
      qm_forces[i + 0 * G2G::fortran_vars.atoms] += prefactor * diff.x;
      qm_forces[i + 1 * G2G::fortran_vars.atoms] += prefactor * diff.y;
      qm_forces[i + 2 * G2G::fortran_vars.atoms] += prefactor * diff.z;
      mm_forces[j + 0 * integral_vars.clatoms] -= prefactor * diff.x;
      mm_forces[j + 1 * integral_vars.clatoms] -= prefactor * diff.y;
      mm_forces[j + 2 * integral_vars.clatoms] -= prefactor * diff.z;
    }
  }
}

template <class scalar_type>
void QMMMIntegral<scalar_type>::calc_nuc_energy(double& Ens) {
  Ens = 0.0;
  for (uint i = 0; i < G2G::fortran_vars.atoms; i++) {
    double3 qm_pos = G2G::fortran_vars.atom_positions(i);
    for (uint j = 0; j < integral_vars.clatoms; j++) {
      double3 mm_pos = integral_vars.clatom_positions(j);
      double3 diff = qm_pos - mm_pos;
      double dist = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
      dist = sqrt(dist);

      double E = integral_vars.clatom_charges(j) *
                 integral_vars.atom_Z(i) / dist;
      Ens += E;
    }
  }
}

#if AINT_MP && !FULL_DOUBLE
template class QMMMIntegral<float>;
#else
template class QMMMIntegral<double>;
#endif
}
