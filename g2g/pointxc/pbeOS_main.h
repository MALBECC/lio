//----------------------------------------------------------------------------//
// EASYPBE is a driver for the PBE subroutines, using simple inputs  
// (K. Burke, May 14, 1996.)
//
// inputs:
//       : dens_a =up density
//       : dgrad_a =|grad up|
//       : delgrad_a =(grad up).(grad |grad up|)
//       : rlap_a =grad^2 up=Laplacian of up
//       : dens_b,dgrad_b,delgrad_b,rlap_b =corresponding down quantities
//       : dgrad =|grad rho|
//       : delgrad =(grad rho).(grad |grad rho|)
//       :?? lcor=flag to do correlation(=0=>don't)
//       :?? lpot=flag to do potential(=0=>don't)
//
// outputs:
//       : expbe,vxpbe_a,vxpbe_b,ecpbe,vcpbe_a,vcpbe_b  =PBE quantities
//
//----------------------------------------------------------------------------//
namespace pointxc {

template<class scalar_type> __host__ __device__
void pbeOS_main( scalar_type dens_a,    scalar_type dgrad_a,
                 scalar_type delgrad_a, scalar_type rlap_a,
                 scalar_type dens_b,    scalar_type dgrad_b,
                 scalar_type delgrad_b, scalar_type rlap_b,
                 scalar_type dgrad,     scalar_type delgrad,
                 scalar_type& expbe,    scalar_type& vxpbe_a,
                 scalar_type& vxpbe_b,  scalar_type& ecpbe,
                 scalar_type& corr1,    scalar_type& corr2,
                 scalar_type& vcpbe_a,  scalar_type& vcpbe_b );

}
