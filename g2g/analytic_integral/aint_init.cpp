/* headers */
#include <iostream>
#include <fstream>
#include "../common.h"
#include "../matrix.h"
#include "../init.h"
#include "../cuda_includes.h"
#include "aint_common.h"
#include "aint_init.h"

#include "os_integral.h"
#include "qmmm_integral.h"
#include "coulomb_integral.h"

#ifndef AINT_GPU_LEVEL
#define AINT_GPU_LEVEL AINT_DEFAULT_GPU_LEVEL
#endif

using std::cout;
using std::endl;
using std::boolalpha;
using std::runtime_error;
using std::ifstream;
using std::string;

using namespace AINT;

namespace AINT {
FortranVars integral_vars;

#if !AINT_MP || FULL_DOUBLE
OSIntegral<double> os_integral;
QMMMIntegral<double> qmmm_integral(os_integral);
CoulombIntegral<double> coulomb_integral(os_integral);
#else
OSIntegral<float> os_integral;
QMMMIntegral<float> qmmm_integral(os_integral);
CoulombIntegral<float> coulomb_integral(os_integral);
#endif
}
//==========================================================================================
extern "C" void aint_parameter_init_(const unsigned int& Md,
                                     unsigned int* ncontd,
                                     const unsigned int* nshelld, double* cd,
                                     double* ad, unsigned int* Nucd, double* af,
                                     double* RMM, const unsigned int& M9,
                                     const unsigned int& M11, double* str,
                                     double* fac, double& rmax, uint* atomZ_i) {
  /* DENSITY BASIS SET */
  integral_vars.s_funcs_dens = nshelld[0];
  integral_vars.p_funcs_dens = nshelld[1] / 3;
  integral_vars.d_funcs_dens = nshelld[2] / 6;
  // Md =	# of contractions
  integral_vars.m_dens = Md;

  if (G2G::verbose > 3) {
     cout << "AINT initialisation." << endl;
     cout << "  Density basis - s: " << integral_vars.s_funcs_dens
          << " p: " << integral_vars.p_funcs_dens
          << " d: " << integral_vars.d_funcs_dens
          << " - Total (w/contractions): " << integral_vars.m_dens << endl;
  }
  integral_vars.spd_funcs_dens = integral_vars.s_funcs_dens +
                                 integral_vars.p_funcs_dens +
                                 integral_vars.d_funcs_dens;  /* DENSITY BASIS SET */
  integral_vars.nucleii_dens =
      G2G::FortranMatrix<uint>(Nucd, integral_vars.m_dens, 1, 1);
  integral_vars.contractions_dens =
      G2G::FortranMatrix<uint>(ncontd, integral_vars.m_dens, 1, 1);
  for (uint i = 0; i < integral_vars.m_dens; i++) {
    if ((integral_vars.contractions_dens(i) - 1) > MAX_CONTRACTIONS)
      throw runtime_error("Maximum functions per contraction reached!");
  }
  uint num_dens_gauss = 0, num_s_gauss = 0, num_p_gauss = 0;
  for (uint i = 0; i < integral_vars.s_funcs_dens; i++) {
    num_dens_gauss += integral_vars.contractions_dens(i);
    num_s_gauss += integral_vars.contractions_dens(i);
  }
  for (uint i = integral_vars.s_funcs_dens;
       i < integral_vars.s_funcs_dens + integral_vars.p_funcs_dens * 3; i++) {
    num_dens_gauss += integral_vars.contractions_dens(i);
    num_p_gauss += integral_vars.contractions_dens(i);
  }
  for (uint i = integral_vars.s_funcs_dens + integral_vars.p_funcs_dens * 3;
       i < integral_vars.m_dens; i++) {
    num_dens_gauss += integral_vars.contractions_dens(i);
  }
  integral_vars.gaussians_dens = num_dens_gauss;
  integral_vars.s_gaussians_dens = num_s_gauss;
  integral_vars.p_gaussians_dens = num_p_gauss;
  integral_vars.a_values_dens = G2G::FortranMatrix<double>(
      ad, integral_vars.m_dens, MAX_CONTRACTIONS, num_dens_gauss);
  integral_vars.c_values_dens = G2G::FortranMatrix<double>(
      cd, integral_vars.m_dens, MAX_CONTRACTIONS, num_dens_gauss);
  // 1e Fock matrix
  integral_vars.rmm_1e_output = G2G::FortranMatrix<double>(
      RMM + (M11 - 1), (G2G::fortran_vars.m * (G2G::fortran_vars.m + 1)) / 2);
  // Inverse of G matrix for Coulomb fitting
  integral_vars.Ginv_input = G2G::FortranMatrix<double>(
      RMM + (M9 - 1), (integral_vars.m_dens * (integral_vars.m_dens + 1)) / 2);
  // Fitted density in density basis
  integral_vars.af_input_ndens1 =
      G2G::FortranMatrix<double>(af, integral_vars.m_dens);
  // Arrays needed in the calculation of F(m,U) (sort of incomplete gamma
  // functions) used in the Obara-Saika recursion relations
  // F(m,U) is calculated by one of two methods, branching based on the value of
  // U: Taylor expansion and asymptotic value
  // STR: Table of precalculated values used in the Taylor expansion branch
  // (TODO: the second index is accessed at most by (m+5),and we only go up to
  // m=5, yet we store up to 22...)
  integral_vars.str = G2G::FortranMatrix<double>(str, 880, 22, 880);
  // FAC: Pre-factor for the different m values in the asymptotic branch (TODO:
  // we only ever go up to m=5, yet FAC is length 17...)
  integral_vars.fac = G2G::FortranMatrix<double>(fac, 17, 1, 17);
  // Maximum Gaussian argument used as basis function overlap cut-off (Coulomb
  // and QM/MM)
  integral_vars.rmax = rmax;
  // Atomic number for each atom. It is NOT the atom type, since Ghost Atoms
  // (basis placed but Z=0) may be used.

  integral_vars.atom_Z =
      G2G::FortranMatrix<uint>(atomZ_i, G2G::fortran_vars.atoms, 1, 1);

#if GPU_KERNELS
  os_integral.load_params();
#endif
}
//============================================================================================================
extern "C" void aint_deinit_(void) {
#if GPU_KERNELS
  int previous_device;
  cudaGetDevice(&previous_device);
  if (cudaSetDevice(os_integral.my_device) != cudaSuccess)
    std::cout << "Error: can't set the device " << os_integral.my_device
              << std::endl;

  os_integral.deinit();
  coulomb_integral.clear();

  cudaSetDevice(previous_device);
#endif
}
//==============================================================================================================
extern "C" void aint_new_step_(void) {
  int stat = 0;
#if GPU_KERNELS
  int previous_device;
  cudaGetDevice(&previous_device);
  if (cudaSetDevice(os_integral.my_device) != cudaSuccess)
    std::cout << "Error: can't set the device " << os_integral.my_device
              << std::endl;

  os_integral.new_cutoff();
  if (!os_integral.load_input()) stat = 1;
  if (!os_integral.alloc_output()) stat = 2;

  cudaSetDevice(previous_device);
#endif

  if (stat != 0) {
    cout << "Could not initialize AINT module; probably due to not enough GPU "
            "memory"
         << endl;
    exit(-1);
  }

  integral_vars.clatoms = 0;
}
//==============================================================================================================
extern "C" void aint_qmmm_init_(const unsigned int& nclatom, double* r_all,
                                double* pc) {
  integral_vars.clatoms = nclatom;
  if (integral_vars.clatoms > 0) {
    integral_vars.clatom_positions_pointer = G2G::FortranMatrix<double>(
        r_all + G2G::fortran_vars.atoms, integral_vars.clatoms, 3,
        G2G::fortran_vars.atoms + nclatom);
    integral_vars.clatom_charges_pointer = G2G::FortranMatrix<double>(
        pc + G2G::fortran_vars.atoms, integral_vars.clatoms, 1,
        G2G::fortran_vars.atoms + nclatom);
  }

  // Read in MM atom information ... locations/charges get sent to the device
  // before the QM/MM kernel
  if (integral_vars.clatoms > 0) {
    integral_vars.clatom_positions.resize(
        integral_vars.clatoms);  // cpu version (double3)
    integral_vars.clatom_charges.resize(integral_vars.clatoms);
    for (uint i = 0; i < integral_vars.clatoms; i++) {
      double3 pos = make_double3(integral_vars.clatom_positions_pointer(i, 0),
                                 integral_vars.clatom_positions_pointer(i, 1),
                                 integral_vars.clatom_positions_pointer(i, 2));
      double charge = integral_vars.clatom_charges_pointer(i);
      integral_vars.clatom_positions(i) = pos;
      integral_vars.clatom_charges(i) = charge;
    }
  }
}
//==============================================================================================================
extern "C" void aint_coulomb_init_(void) {
#if GPU_KERNELS
  int previous_device;
  cudaGetDevice(&previous_device);
  if (cudaSetDevice(os_integral.my_device) != cudaSuccess)
    std::cout << "Error: can't set the device " << os_integral.my_device
              << std::endl;

  int stat = 0;
  coulomb_integral.clear();
  if (!coulomb_integral.load_aux_basis()) stat = 1;
  if (!coulomb_integral.load_input()) stat = 2;
  if (!coulomb_integral.alloc_output()) stat = 3;
  if (stat != 0) {
    cout << "Could not initialize COULOMB module; probably due to not enough "
            "GPU memory"
         << endl;
    exit(-1);
  }

  cudaSetDevice(previous_device);
#endif
}
//===============================================================================================================
//                                  QM/MM routines
//===============================================================================================================
extern "C" void aint_qmmm_forces_(double* qm_forces, double* mm_forces) {
#if GPU_KERNELS
  int previous_device;
  cudaGetDevice(&previous_device);
  if (cudaSetDevice(os_integral.my_device) != cudaSuccess)
    std::cout << "Error: can't set the device " << os_integral.my_device
              << std::endl;

  int stat = 0;
  if (integral_vars.clatoms > 0) {
    if (!qmmm_integral.load_clatoms()) stat = 1;
    if (!qmmm_integral.alloc_output()) stat = 2;
  }
  if (stat != 0) {
    cout << "Could not initialize QM/MM module; probably due to not enough GPU "
            "memory"
         << endl;
    exit(-1);
  }

  qmmm_integral.calc_nuc_gradient(qm_forces, mm_forces);
  if (AINT_GPU_LEVEL > 3) {
    if (integral_vars.clatoms > 0) {
      qmmm_integral.calc_gradient(qm_forces, mm_forces, true, true);
    } else {
      qmmm_integral.calc_gradient(qm_forces, mm_forces, false, true);
    }
  } else {
    if (integral_vars.clatoms > 0) {
      qmmm_integral.calc_gradient(qm_forces, mm_forces, true, false);
    } else {
      qmmm_integral.calc_gradient(qm_forces, mm_forces, false,
                                  false);  // This branch does nothing...
    }
  }

  qmmm_integral.clear();

  cudaSetDevice(previous_device);
#endif
}
extern "C" void aint_qmmm_fock_(double& Es, double& Ens) {
  Ens = 0.0;
  Es = 0.0;
#if GPU_KERNELS
  int previous_device;
  cudaGetDevice(&previous_device);
  if (cudaSetDevice(os_integral.my_device) != cudaSuccess)
    std::cout << "Error: can't set the device " << os_integral.my_device
              << std::endl;

  int stat = 0;
  if (integral_vars.clatoms > 0) {
    if (!qmmm_integral.load_clatoms()) stat = 1;
  }
  if (stat != 0) {
    cout << "Could not initialize QM/MM module; probably due to not enough GPU "
            "memory"
         << endl;
    exit(-1);
  }

  qmmm_integral.calc_nuc_energy(Ens);
  if (AINT_GPU_LEVEL > 3) {
    if (integral_vars.clatoms > 0) {
      qmmm_integral.calc_fock(Es, true, true);
    } else {
      qmmm_integral.calc_fock(Es, false, true);
    }
  } else {
    if (integral_vars.clatoms > 0) {
      qmmm_integral.calc_fock(Es, true, false);
    } else {
      qmmm_integral.calc_fock(Es, false, false);
    }
  }

  qmmm_integral.clear();

  cudaSetDevice(previous_device);
#endif
}
//===============================================================================================================
//                                  Coulomb routines
//===============================================================================================================
extern "C" void aint_coulomb_forces_(double* qm_forces) {
#if GPU_KERNELS
  int previous_device;
  cudaGetDevice(&previous_device);
  if (cudaSetDevice(os_integral.my_device) != cudaSuccess)
    std::cout << "Error: can't set the device " << os_integral.my_device
              << std::endl;
  if (AINT_GPU_LEVEL >= 5) {
    coulomb_integral.calc_gradient(qm_forces, false);
  } else {
    coulomb_integral.calc_gradient(qm_forces, true);
  }
  cudaSetDevice(previous_device);
#endif
}
extern "C" void aint_coulomb_fock_(double& Es) {
  Es = 0.0;
#if GPU_KERNELS
  int previous_device;
  cudaGetDevice(&previous_device);
  if (cudaSetDevice(os_integral.my_device) != cudaSuccess)
    std::cout << "Error: can't set the device " << os_integral.my_device
              << std::endl;

  coulomb_integral.fit_aux_density();
  coulomb_integral.calc_fock(Es);

  cudaSetDevice(previous_device);
#endif
}
//===============================================================================================================
extern "C" void aint_query_gpu_level_(int& gpu_level) {
#if !GPU_KERNELS
  gpu_level = 0;
#else

  // 0 - No GPU acceleration
  // 1 - Only XC on GPU (no difference from 0 for us here)
  // 2 - 1 + QM/MM energy and gradient terms on GPU
  // 3 - 1-2 + Coulomb gradient terms on GPU
  // 4 - 1-3 + Nuclear attraction gradient terms on GPU
  // 5 - 1-4 + Coulomb auxiliary basis fitting and energy terms on GPU
  if (AINT_GPU_LEVEL < 0) {
    gpu_level = 0;
  } else if (AINT_GPU_LEVEL > 5) {
    gpu_level = 5;
  } else {
    gpu_level = AINT_GPU_LEVEL;
  }
#endif
}
//===============================================================================================================
