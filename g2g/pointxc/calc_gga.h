//namespace {
#include "../matrix.h"
//#include "../scalar_vector_types.h"
using G2G::vec_type;
//}

namespace pointxc {

template<class scalar_type> __host__ __device__
void calc_ggaOS( scalar_type dens_a, scalar_type dens_b,
                 const vec_type<scalar_type,4>& grad_a,
                 const vec_type<scalar_type,4>& grad_b,
                 const vec_type<scalar_type,4>& hess1_a,
                 const vec_type<scalar_type,4>& hess1_b,
                 const vec_type<scalar_type,4>& hess2_a,
                 const vec_type<scalar_type,4>& hess2_b,
                 scalar_type& exc_corr, scalar_type& exc,
                 scalar_type& corr,     scalar_type& corr1,
                 scalar_type& corr2,    scalar_type& v_a,
                 scalar_type& v_b );

template<class scalar_type, int iexch> __host__ __device__
void calc_ggaCS( scalar_type dens, const vec_type<scalar_type,3>& grad,
                 const vec_type<scalar_type,3>& hess1,
                 const vec_type<scalar_type,3>& hess2,
                 scalar_type& ex, scalar_type& ec, scalar_type& y2a);

template<class scalar_type>
void calc_ggaCS_in( scalar_type dens, const vec_type<scalar_type,3>& grad,
                    const vec_type<scalar_type,3>& hess1,
                    const vec_type<scalar_type,3>& hess2,
                    scalar_type& ex, scalar_type& ec, scalar_type& y2a,
                    const int iexch);

}
