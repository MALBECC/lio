#include "../common.h"

#include "pbeOS_main.h"
#include "../fix_compile.h"
#include "../scalar_vector_types.h"

namespace G2G {

template <class scalar_type, unsigned int width>
__host__ __device__ void calc_ggaOS(scalar_type dens_a, scalar_type dens_b,
                                    const vec_type<scalar_type, width>& grad_a,
                                    const vec_type<scalar_type, width>& grad_b,
                                    const vec_type<scalar_type, width>& hess1_a,
                                    const vec_type<scalar_type, width>& hess1_b,
                                    const vec_type<scalar_type, width>& hess2_a,
                                    const vec_type<scalar_type, width>& hess2_b,
                                    scalar_type& exc_corr, scalar_type& exc,
                                    scalar_type& corr, scalar_type& corr1,
                                    scalar_type& corr2, scalar_type& v_a,
                                    scalar_type& v_b, int Vxc_id) {
  scalar_type expbe = 0.0f;
  scalar_type vxpbe_a = 0.0f;
  scalar_type vxpbe_b = 0.0f;
  scalar_type ecpbe = 0.0f;
  scalar_type vcpbe_a = 0.0f;
  scalar_type vcpbe_b = 0.0f;

  exc_corr = exc = corr = corr1 = corr2 = v_a = v_b = 0.0f;

  scalar_type dgrad_a, delgrad_a, rlap_a, dgrad_b, delgrad_b, rlap_b;
  const scalar_type MINIMUM_DENSITY_VALUE = 1e-13f;

  if ((dens_a + dens_b) > MINIMUM_DENSITY_VALUE) {
    if (Vxc_id == 9) {  // PBE

      vec_type<scalar_type, width> grad;
      grad.x = grad_a.x + grad_b.x;
      grad.y = grad_a.y + grad_b.y;
      grad.z = grad_a.z + grad_b.z;

      vec_type<scalar_type, width> hess1;
      hess1.x = hess1_a.x + hess1_b.x;
      hess1.y = hess1_a.y + hess1_b.y;
      hess1.z = hess1_a.z + hess1_b.z;

      vec_type<scalar_type, width> hess2;
      hess2.x = hess2_a.x + hess2_b.x;
      hess2.y = hess2_a.y + hess2_b.y;
      hess2.z = hess2_a.z + hess2_b.z;

      // Up density
      if (dens_a == ((scalar_type)0.0f)) {
        dgrad_a = 0.0f;
        rlap_a = 0.0f;
        delgrad_a = 0.0f;
      } else {
        scalar_type grad2_a =
            pow(grad_a.x, 2) + pow(grad_a.y, 2) + pow(grad_a.z, 2);
        if (grad2_a == (scalar_type)0.0f) grad2_a = MIN_PRECISION;
        dgrad_a = sqrt(grad2_a);
        rlap_a = hess1_a.x + hess1_a.y + hess1_a.z;  // Laplacian Up

        delgrad_a =
            (pow(grad_a.x, 2) * hess1_a.x +
             (scalar_type)2.0f * grad_a.x * grad_a.y * hess2_a.x +
             (scalar_type)2.0f * grad_a.y * grad_a.z * hess2_a.z +
             (scalar_type)2.0f * grad_a.x * grad_a.z * hess2_a.y +
             pow(grad_a.y, 2) * hess1_a.y + pow(grad_a.z, 2) * hess1_a.z) /
            dgrad_a;
      }

      // Down density
      if (dens_b == ((scalar_type)0.0f)) {
        dgrad_b = 0.0f;
        rlap_b = 0.0f;
        delgrad_b = 0.0f;
      } else {
        scalar_type grad2_b =
            pow(grad_b.x, 2) + pow(grad_b.y, 2) + pow(grad_b.z, 2);
        if (grad2_b == (scalar_type)0.0f) grad2_b = MIN_PRECISION;
        dgrad_b = sqrt(grad2_b);
        rlap_b = hess1_b.x + hess1_b.y + hess1_b.z;

        delgrad_b =
            (pow(grad_b.x, 2) * hess1_b.x +
             (scalar_type)2.0f * grad_b.x * grad_b.y * hess2_b.x +
             (scalar_type)2.0f * grad_b.y * grad_b.z * hess2_b.z +
             (scalar_type)2.0f * grad_b.x * grad_b.z * hess2_b.y +
             pow(grad_b.y, 2) * hess1_b.y + pow(grad_b.z, 2) * hess1_b.z) /
            dgrad_b;
      }

      // Up + Down densities
      scalar_type grad2 = pow(grad.x, 2) + pow(grad.y, 2) + pow(grad.z, 2);
      if (grad2 == (scalar_type)0.0f) grad2 = MIN_PRECISION;
      scalar_type dgrad = sqrt(grad2);
      scalar_type delgrad =
          (pow(grad.x, 2) * hess1.x + pow(grad.y, 2) * hess1.y +
           pow(grad.z, 2) * hess1.z +
           (scalar_type)2.0f * grad.x * grad.y * hess2.x +
           (scalar_type)2.0f * grad.y * grad.z * hess2.z +
           (scalar_type)2.0f * grad.x * grad.z * hess2.y) /
          dgrad;

      pbeOS_main<scalar_type>(dens_a, dgrad_a, delgrad_a, rlap_a, dens_b,
                              dgrad_b, delgrad_b, rlap_b, dgrad, delgrad, expbe,
                              vxpbe_a, vxpbe_b, ecpbe, corr1, corr2, vcpbe_a,
                              vcpbe_b);
    } else {
      // NO HAY IMPLEMENTADO OTRO FUNCIONAL DE XC
    }
  }

  exc = expbe;
  corr = ecpbe;
  exc_corr = expbe + ecpbe;
  v_a = vxpbe_a + vcpbe_a;
  v_b = vxpbe_b + vcpbe_b;
  return;
}
}
