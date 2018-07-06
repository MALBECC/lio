#ifndef __AINTEGRAL_H__
#define __AINTEGRAL_H__

#include "../matrix.h"

namespace AINT {

struct FortranVars {
  uint clatoms;

  /* ---DENSITY BASIS--- */
  uint gaussians_dens, s_gaussians_dens, p_gaussians_dens;
  uint s_funcs_dens, p_funcs_dens, d_funcs_dens, spd_funcs_dens, m_dens;
  G2G::FortranMatrix<uint> atom_Z;
  G2G::FortranMatrix<uint> nucleii_dens, contractions_dens;
  G2G::FortranMatrix<double> a_values_dens, c_values_dens;
  G2G::FortranMatrix<double> af_input_ndens1;
  /* ---DENSITY BASIS--- */

  G2G::FortranMatrix<double> clatom_positions_pointer;
  G2G::FortranMatrix<double> clatom_charges_pointer;
  G2G::HostMatrix<double3> clatom_positions;
  G2G::HostMatrix<double> clatom_charges;
  G2G::FortranMatrix<double> rmm_1e_output;
  G2G::FortranMatrix<double> Ginv_input;
  // Arrays used for numerical evaluation of F(m,U) functions in Obara-Saika
  // recursion
  G2G::FortranMatrix<double> str, fac;
  // Gaussian argument cut-off for QM/MM, Coulomb
  double rmax;
};

extern AINT::FortranVars integral_vars;

extern uint NUM_TERM_TYPES;
}

#endif
