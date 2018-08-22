//----------------------------------------------------------------------------//
// EASYPBE is a driver for the PBE subroutines, using simple inputs
// (K. Burke, May 14, 1996.)
//
// inputs:
//       : dens_a    = up density, dgrad_a = |grad up|
//       : delgrad_a = (grad up).(grad |grad up|)
//       : rlap_a    = Laplacian of up
//       : dens_b, dgrad_b, delgrad_b, rlap_b =corresponding down quantities
//       : dgrad   = |grad rho|
//       : delgrad = (grad rho).(grad |grad rho|)
//       : lcor    = Flag to do correlation(0 = do not)
//       : lpot    = Flag to do potential  (0 = do not)
//
// outputs: (PBE quantities)
//       : expbe, vxpbe_a, vxpbe_b, ecpbe, vcpbe_a, vcpbe_b
//
//----------------------------------------------------------------------------//
#include "pbeOS_corr.h"
#include "pbeOS_exch.h"
#include "../fix_compile.h"
#include "../common.h"

namespace G2G {
#define EASYPBE_PI ((scalar_type)3.141592653589793238462643383279f)
#define EASYPBE_PI32 ((scalar_type)29.608813203268075856503472999628f)

template <class scalar_type>
__host__ __device__ void pbeOS_main(
    scalar_type dens_a, scalar_type dgrad_a, scalar_type delgrad_a,
    scalar_type rlap_a, scalar_type dens_b, scalar_type dgrad_b,
    scalar_type delgrad_b, scalar_type rlap_b, scalar_type dgrad,
    scalar_type delgrad, scalar_type& expbe, scalar_type& vxpbe_a,
    scalar_type& vxpbe_b, scalar_type& ecpbe, scalar_type& corr1,
    scalar_type& corr2, scalar_type& vcpbe_a, scalar_type& vcpbe_b) {
  /*-----------------------------------------------------//
  // PBE exchange
  //-----------------------------------------------------//
  // Use exact spin-scaling:
  //    Ex[up, dn] = 0.5*(Ex[2*up] + Ex[2*dn])
  // Perdew et al. Phys Rev Lett 77,3865. 1996 (Eq. 11)
  //
  // Local Fermi wavevector for 2*up
  // fk = (3 pi^2 (2up))^(1/3)
  //
  // Dimensionless density gradient
  // s = |grad rho|/ (2*fk*rho)    where (rho=2*up)
  //
  // Other magnitudes
  // u = delgrad/(rho^2*(2*fk)**3) where (rho=2*up)
  // v = Laplacian/(rho*(2*fk)**2) where (rho=2*up)
  //-----------------------------------------------------*/
   scalar_type expbe_a, expbe_b;
   scalar_type rho13, fk1, fk, twofk, twofk2, twofk3, s, u, v;
   bool small_dens = false;

   // Output initialization
   expbe = (scalar_type)0.0f;
   ecpbe = (scalar_type)0.0f;
   corr1 = (scalar_type)0.0f;
   corr2 = (scalar_type)0.0f;
   vcpbe_a = (scalar_type)0.0f;
   vcpbe_b = (scalar_type)0.0f;

   // Density Up
   scalar_type twodens = (scalar_type)2.0f * dens_a;
   scalar_type twodens2 = twodens * twodens;
   scalar_type twodens5 = twodens2 * twodens2 * twodens;

   scalar_type flt_minimum = (scalar_type) (100.0 * MIN_PRECISION);
   if (twodens5 < flt_minimum) {
     expbe_a = (scalar_type)0.0f;
     vxpbe_a = (scalar_type)0.0f;
     small_dens = true;
   } else {
     rho13 = cbrt(twodens);
     fk1 = cbrt((scalar_type)EASYPBE_PI32);
     fk = fk1 * rho13;

     twofk = (scalar_type)2.0f * fk;
     twofk2 = twofk * twofk;
     twofk3 = twofk * twofk2;

     s = ((scalar_type)2.0f * dgrad_a)   / (twodens * twofk);
     v = ((scalar_type)2.0f * rlap_a)    / (twodens * twofk2);
     u = ((scalar_type)4.0f * delgrad_a) / (twodens2 * twofk3);
     pbeOS_exch(twodens, s, u, v, expbe_a, vxpbe_a);
   }

   // Density Down
   twodens = (scalar_type)2.0f * dens_b;
   twodens2 = twodens * twodens;
   twodens5 = twodens2 * twodens2 * twodens;

   if (twodens5 < flt_minimum) {
     expbe_a = (scalar_type)0.0f;
     vxpbe_a = (scalar_type)0.0f;
     ecpbe   = (scalar_type)0.0f;
     if (!small_dens) {
        expbe = expbe_a;
     }
     return;
   } else {
     rho13 = cbrt((scalar_type)twodens);
     fk1 = cbrt((scalar_type)EASYPBE_PI32);
     fk = fk1 * rho13;

     twofk = (scalar_type)2.0f * fk;
     twofk2 = twofk * twofk;
     twofk3 = twofk * twofk2;
 
     s = ((scalar_type)2.0f * dgrad_b)   / (twodens * twofk);
     v = ((scalar_type)2.0f * rlap_b)    / (twodens * twofk2);
     u = ((scalar_type)4.0f * delgrad_b) / (twodens2 * twofk3);
     pbeOS_exch(twodens, s, u, v, expbe_b, vxpbe_b);
   }

   if (small_dens) {
      expbe = expbe_b;
      ecpbe = (scalar_type)0.0f;
      return;
   }

   // Construct total density and contribution to Ex
   scalar_type rho = dens_a + dens_b;
   expbe = (expbe_a * dens_a + expbe_b * dens_b) / rho;

  /*-------------------------------------------------------------//
  // PBE correlation
  //-------------------------------------------------------------//
  // zet=(up-dn)/rho || g=phi(zeta)
  //
  // Local Seitz radius (alpha/fk)
  // rs = (3/(4pi*rho))^(1/3)
  //
  // Thomas-Fermi screening wavevector
  // sk = Ks = sqrt(4fk/pi) || twoksg = 2*Ks*phi
  //
  // Correlation dimensionless gradient
  // t = |grad rho| / (2*Ks*phi*rho)
  //
  // rholap=Laplacian
  // uu =delgrad/(rho^2*twoksg^3) || vv = Laplacian/(rho*twoksg^2)
  // ww = (|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
  //
  // ec    = lsd correlation energy
  // vcup  = lsd up correlation potential
  // vcdn  = lsd down correlation potential
  // h     = Gradient correction to correlation energy
  // dvcup = Gradient correction to up correlation potential
  // dvcdn = Gradient correction to down correlation potential
  //--------------------------------------------------------------*/

  scalar_type zet = (dens_a - dens_b) / rho;
  scalar_type g1 = (scalar_type)1.0f + zet;
  scalar_type g1_23 = cbrt(g1 * g1);
  scalar_type g2 = (scalar_type)1.0f - zet;
  scalar_type g2_23 = cbrt(g2 * g2);
  scalar_type g = (g1_23 + g2_23) / (scalar_type)2.0f;

  rho13 = cbrt(rho);
  fk1 = cbrt(EASYPBE_PI32);
  fk = fk1 * rho13;

  scalar_type alpha =
      ((scalar_type)9.0f * (scalar_type)EASYPBE_PI) / (scalar_type)4.0f;
  scalar_type alpha13 = cbrt(alpha);
  scalar_type rs = alpha13 / fk;

  scalar_type sk = sqrt(4.0f * fk / (scalar_type)EASYPBE_PI);
  scalar_type twoksg = (scalar_type)2.0f * sk * g;
  scalar_type t = dgrad / (twoksg * rho);

  scalar_type uu = delgrad / (rho * rho * twoksg * twoksg * twoksg);
  scalar_type rho_lap = rlap_a + rlap_b;
  scalar_type vv = rho_lap / (rho * twoksg * twoksg);

  // Reconstruction of WW
  scalar_type num =
      (dgrad_a * dgrad_a) - (dgrad_b * dgrad_b) - (zet * dgrad * dgrad);
  scalar_type den = (rho * rho) * (twoksg * twoksg);
  scalar_type ww = num / den;

  scalar_type ec, h;
  scalar_type vc_a, vc_b, dvc_a, dvc_b;

  pbeOS_corr<scalar_type>(rho, rs, zet, t, uu, vv, ww, ec, vc_a, vc_b, h, dvc_a,
                          dvc_b);

  ecpbe = ec + h;
  corr1 = ec;
  corr2 = h;
  vcpbe_a = vc_a + dvc_a;
  vcpbe_b = vc_b + dvc_b;
//  if (expbe_a != expbe_a) { printf("NaN in expbe_a \n");};
//  if (expbe_b != expbe_b) { printf("NaN in expbe_b \n");};
//  if (ec != ec) { printf("NaN in ec \n");};
//  if (h != h) { printf("NaN in h \n");};
//  if (vc_a != vc_a) { printf("NaN in vc_a \n");};
//  if (dvc_a != dvc_a) { printf("NaN in dvc_a \n");};
//  if (vc_b != vc_b) { printf("NaN in vc_b \n");};
//  if (dvc_b!= dvc_b) { printf("NaN in dvc_b \n");};
}  // pbeOS_main

#undef EASYPBE_PI
#undef EASYPBE_PI32
}
