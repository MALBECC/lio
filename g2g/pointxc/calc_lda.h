namespace pointxc {

template<class scalar_type, int iexch> __host__ __device__
void calc_ldaCS( scalar_type dens, scalar_type& ex, 
                   scalar_type& ec,  scalar_type& y2a );

template<class scalar_type> __host__ __device__
void calc_ldaCS_in( scalar_type dens, scalar_type& ex,
                      scalar_type& ec,  scalar_type& y2a,
                      const int iexch);

}
