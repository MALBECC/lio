
namespace pointxc {

template<class scalar_type>
__host__ __device__
static void pbeCS_v1( scalar_type rho,
                      scalar_type agrad,
                      scalar_type delgrad,
                      scalar_type rlap,
                      scalar_type& expbe,
                      scalar_type& vxpbe,
                      scalar_type& ecpbe,
                      scalar_type& vcpbe );

}
