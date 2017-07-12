
namespace pointxc {

template<class scalar_type> __host__ __device__
static void gcorc1( scalar_type rtrs, scalar_type& gg, scalar_type& grrs );

template<class scalar_type> __host__ __device__
static void gcorc2( scalar_type rtrs, scalar_type& gg, scalar_type& grrs );

template<class scalar_type> __host__ __device__
static void gcorc3( scalar_type rtrs, scalar_type& gg, scalar_type& grrs );

}
