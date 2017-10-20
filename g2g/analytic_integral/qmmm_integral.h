#ifndef __QMMM_H__
#define __QMMM_H__

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include "../scalar_vector_types.h"
#include "../timer.h"
#include "os_integral.h"

namespace AINT {

template<class scalar_type>
class QMMMIntegral {
  public:

    //
    // Common data
    //
    OSIntegral<scalar_type>& os_int;

#if GPU_KERNELS
    //
    // On device
    //
    // QM/MM-only input for O-S integral evaluation
    G2G::CudaMatrix<G2G::vec_type<scalar_type, 3> > clatom_pos_dev;
    G2G::CudaMatrix<scalar_type> clatom_chg_dev;

    // QM/MM-only output for O-S integral evaluation (gradient)
    G2G::CudaMatrix<G2G::vec_type<scalar_type,3> > partial_mm_forces_dev;
#endif

    //
    // Functions
    //
    QMMMIntegral( OSIntegral<scalar_type>& _os_int ) : os_int(_os_int) { }

    // Nuclear-nuclear interactions
    void calc_nuc_energy( double& Ens );
    void calc_nuc_gradient( double* qm_forces, double* mm_forces );

    bool load_clatoms( void );
    bool alloc_output( void );

    // Main kernel calls
    void calc_fock( double& Es, bool do_cl, bool do_qm );

    void calc_gradient( double* qm_forces, double* mm_forces, bool do_cl, bool do_qm );

    // Get results from GPU and send to final output
    void get_gradient_output( double* mm_forces, uint partial_size );

    void clear( void );
    
};

}

#endif
