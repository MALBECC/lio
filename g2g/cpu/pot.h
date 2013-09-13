#ifndef _CPU_POT_H
#define	_CPU_POT_H

namespace G2G {
template<class scalar_type> void cpu_pot(scalar_type dens, scalar_type& ex, scalar_type& ec, scalar_type& y2a);
template<class scalar_type> void cpu_potg(scalar_type dens, const vec_type<scalar_type,3>& grad, const vec_type<scalar_type,3>& hess1,
                                          const vec_type<scalar_type,3>& hess2, scalar_type& ex, scalar_type& ec, scalar_type& y2a);
}


#endif	/* _CPU_POT_H */

