/* headers */
#include <iostream>
#include <fstream>
#include "../common.h"
#include "../matrix.h"
#include "../init.h"
#include "aint_common.h"
#include "aint_init.h"

#include "os_integral.h"
#include "qmmm_integral.h"
#include "coulomb_integral.h"

using std::cout;
using std::endl;
using std::boolalpha;
using std::runtime_error;
using std::ifstream;
using std::string;

using namespace AINT;

namespace AINT
{
	FortranVars integral_vars;

#if !AINT_MP || FULL_DOUBLE
	OSIntegral<double> os_integral;
	QMMMIntegral<double> qmmm_integral( os_integral );
	CoulombIntegral<double> coulomb_integral( os_integral );
#else
	OSIntegral<float> os_integral;
	QMMMIntegral<float> qmmm_integral( os_integral );
	CoulombIntegral<float> coulomb_integral( os_integral );
#endif

	//void gpu_set_os_variables(void);
	//template<class scalar_type> void clean_gamma(void);
	//template<class scalar_type> void gpu_set_gamma_arrays(void);
}
//==========================================================================================
extern "C" void aint_parameter_init_(const unsigned int& Md, unsigned int* ncontd, const unsigned int* nshelld, double* cd, double* ad, unsigned int* Nucd, double* af,
                                     double* RMM, const unsigned int& M9, const unsigned int& M11, double* str, double* fac, double& rmax)
{
        /* DENSITY BASIS SET */
	integral_vars.s_funcs_dens = nshelld[0];
	integral_vars.p_funcs_dens = nshelld[1] / 3;
	integral_vars.d_funcs_dens = nshelld[2] / 6;
        cout << "density basis: s: " << integral_vars.s_funcs_dens  << " p: " << integral_vars.p_funcs_dens << " d: " << integral_vars.d_funcs_dens << endl;
	integral_vars.spd_funcs_dens = integral_vars.s_funcs_dens + integral_vars.p_funcs_dens + integral_vars.d_funcs_dens;
// Md =	# of contractions
	integral_vars.m_dens = Md;
        cout << "density basis: m: " << integral_vars.m_dens << endl;
        /* DENSITY BASIS SET */
	integral_vars.nucleii_dens = G2G::FortranMatrix<uint>(Nucd, integral_vars.m_dens, 1, 1);
	integral_vars.contractions_dens = G2G::FortranMatrix<uint>(ncontd, integral_vars.m_dens, 1, 1);
        for (uint i = 0; i < integral_vars.m_dens; i++) {
        	if ((integral_vars.contractions_dens(i) - 1) > MAX_CONTRACTIONS)  throw runtime_error("Maximum functions per contraction reached!");
 	}
        uint num_dens_gauss = 0,num_s_gauss=0,num_p_gauss=0;
        for (uint i = 0; i < integral_vars.s_funcs_dens; i++) {
          num_dens_gauss += integral_vars.contractions_dens(i);
          num_s_gauss += integral_vars.contractions_dens(i);
        }
        for (uint i = integral_vars.s_funcs_dens; i < integral_vars.s_funcs_dens+integral_vars.p_funcs_dens*3; i++) {
          num_dens_gauss += integral_vars.contractions_dens(i);
          num_p_gauss += integral_vars.contractions_dens(i);
        }
        for (uint i = integral_vars.s_funcs_dens+integral_vars.p_funcs_dens*3; i < integral_vars.m_dens; i++) {
          num_dens_gauss += integral_vars.contractions_dens(i);
        }
        integral_vars.gaussians_dens = num_dens_gauss;
        integral_vars.s_gaussians_dens = num_s_gauss;
        integral_vars.p_gaussians_dens = num_p_gauss;
	integral_vars.a_values_dens = G2G::FortranMatrix<double>(ad, integral_vars.m_dens, MAX_CONTRACTIONS, num_dens_gauss);
	integral_vars.c_values_dens = G2G::FortranMatrix<double>(cd, integral_vars.m_dens, MAX_CONTRACTIONS, num_dens_gauss);
        // 1e Fock matrix
	integral_vars.rmm_1e_output = G2G::FortranMatrix<double>(RMM + (M11 - 1), (G2G::fortran_vars.m * (G2G::fortran_vars.m + 1)) / 2);
        // Inverse of G matrix for Coulomb fitting
        integral_vars.Ginv_input = G2G::FortranMatrix<double>(RMM + (M9 - 1), (integral_vars.m_dens * (integral_vars.m_dens + 1)) / 2);
        // Fitted density in density basis
        integral_vars.af_input_ndens1 = G2G::FortranMatrix<double>(af, integral_vars.m_dens);
        // Arrays needed in the calculation of F(m,U) (sort of incomplete gamma functions) used in the Obara-Saika recursion relations
        // F(m,U) is calculated by one of two methods, branching based on the value of U: Taylor expansion and asymptotic value
        // STR: Table of precalculated values used in the Taylor expansion branch (TODO: the second index is accessed at most by (m+5),and we only go up to m=5, yet we store up to 22...)
        integral_vars.str = G2G::FortranMatrix<double>(str, 880, 22, 880);
        // FAC: Pre-factor for the different m values in the asymptotic branch (TODO: we only ever go up to m=5, yet FAC is length 17...)
        integral_vars.fac = G2G::FortranMatrix<double>(fac, 17, 1, 17);
        // Maximum Gaussian argument used as basis function overlap cut-off (Coulomb and QM/MM)
        integral_vars.rmax = rmax;

        os_integral.load_params();

}
//============================================================================================================
extern "C" void aint_deinit_( void )
{
	os_integral.deinit();
	qmmm_integral.clear();
	coulomb_integral.clear();
}
//==============================================================================================================
extern "C" void aint_new_step_( void )
{
	os_integral.new_cutoff();
	os_integral.load_input();
	os_integral.alloc_output();

        integral_vars.clatoms = 0;
}
//==============================================================================================================
extern "C" void aint_qmmm_init_( const unsigned int& nclatom, double* r_all, double* pc )
{

	integral_vars.clatoms = nclatom;
	cout << "MM point charges: " << integral_vars.clatoms << endl;
        if (integral_vars.clatoms > 0) {
	    integral_vars.clatom_positions_pointer = G2G::FortranMatrix<double>(r_all + G2G::fortran_vars.atoms, integral_vars.clatoms, 3, G2G::fortran_vars.atoms+nclatom);
	    integral_vars.clatom_charges_pointer = G2G::FortranMatrix<double>(pc + G2G::fortran_vars.atoms, integral_vars.clatoms, 1, G2G::fortran_vars.atoms+nclatom);
        }

        // Read in MM atom information ... locations/charges get sent to the device before the QM/MM kernel
        if (integral_vars.clatoms > 0) {
	    integral_vars.clatom_positions.resize(integral_vars.clatoms);	// cpu version (double3)
  	    integral_vars.clatom_charges.resize(integral_vars.clatoms);
	    for (uint i = 0; i < integral_vars.clatoms; i++) {
		double3 pos = make_double3(integral_vars.clatom_positions_pointer(i, 0), integral_vars.clatom_positions_pointer(i, 1), integral_vars.clatom_positions_pointer(i, 2));
                double charge = integral_vars.clatom_charges_pointer(i);
		integral_vars.clatom_positions(i) = pos;
                integral_vars.clatom_charges(i) = charge;
	    }
  	    qmmm_integral.load_clatoms();
	    qmmm_integral.alloc_output();

        }
}
//==============================================================================================================
extern "C" void aint_coulomb_init_( void )
{
        coulomb_integral.load_aux_basis();
	coulomb_integral.load_input();
	coulomb_integral.alloc_output();
}
//===============================================================================================================
//                                  QM/MM routines
//===============================================================================================================
extern "C" void aint_qmmm_forces_(double* qm_forces, double* mm_forces)
{
  qmmm_integral.calc_nuc_gradient(qm_forces, mm_forces);
  qmmm_integral.calc_gradient(qm_forces, mm_forces);
}
extern "C" void aint_qmmm_fock_(double& Es, double& Ens)
{
  Ens = 0.0; Es = 0.0;
  qmmm_integral.calc_nuc_energy(Ens);
  qmmm_integral.calc_fock(Es);
}
//===============================================================================================================
//                                  Coulomb routines
//===============================================================================================================
extern "C" void aint_coulomb_forces_(double* qm_forces)
{
  coulomb_integral.calc_gradient(qm_forces);
}
extern "C" void aint_coulomb_fock_(double& Es)
{
  Es = 0.0;
  coulomb_integral.fit_aux_density( );
  coulomb_integral.calc_fock(Es);
}
//===============================================================================================================
extern "C" void aint_query_gpu_level_(int& gpu_level)
{
#if CPU_KERNELS
    gpu_level = 0;
#else

// 0 - No GPU acceleration
// 1 - Only XC on GPU (no difference from 0 for us here)
// 2 - 1 + QM/MM energy and gradient terms on GPU
// 3 - 1-2 + Coulomb gradient terms on GPU
// 4 - 1-3 + Nuclear attraction gradient terms on GPU
// 5 - 1-4 + Coulomb auxiliary basis fitting and energy terms on GPU
#ifdef AINT_GPU_LEVEL
    if (AINT_GPU_LEVEL < 0)      { gpu_level = 0; }
    else if (AINT_GPU_LEVEL > 5) { gpu_level = 5; }
    else                         { gpu_level = AINT_GPU_LEVEL; }
#else
    gpu_level = 2;
#endif

#endif  
}
//===============================================================================================================
