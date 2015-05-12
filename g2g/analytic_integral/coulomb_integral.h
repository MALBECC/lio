#ifndef __COULOMB_H__
#define __COULOMB_H__

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include "../scalar_vector_types.h"
#include "../timer.h"
#include "os_integral.h"

namespace AINT {

template<class scalar_type>
class CoulombIntegral {
  public:

    //
    // Common data
    //
    OSIntegral<scalar_type>& os_int;

    //
    // On host
    //
    uint s_end;
    uint p_offset, p_end;
    uint d_offset, d_end;
    uint input_size;

    std::vector<uint> input_ind_cpu;

#if !CPU_KERNELS
    //
    // On device
    //
    G2G::CudaMatrix<scalar_type> fit_dens_dev;

    // Coulomb-only input for O-S integral evaluation
    G2G::CudaMatrix<G2G::vec_type<scalar_type, 2> > factor_ac_dens_dev;
    G2G::CudaMatrix<G2G::vec_type<scalar_type, 3> > nuc_dens_dev;
    G2G::CudaMatrix<uint> nuc_ind_dens_dev;
    G2G::CudaMatrix<uint> input_ind_dev;
    G2G::CudaMatrix<double> Ginv_dev;
    G2G::CudaMatrix<double> rc_partial_dev;
#endif

    //
    // Functions
    //
    CoulombIntegral( OSIntegral<scalar_type>& _os_int ) : os_int(_os_int) { }

    // Allocate Coulomb auxiliary basis on GPU
    // (Needs to be done only once)
    void load_aux_basis( void );

    // Allocate Coulomb-only input on GPU
    void load_input( void );
    void alloc_output( void );

    // Main kernel calls
    void fit_aux_density( void );
    void calc_fock( double& Es );
    void calc_gradient( double* qm_forces, bool cpu_fit_dens );

    void clear( void );
    
};

}

#endif
