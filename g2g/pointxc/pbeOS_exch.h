//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
//  K Burke's modification of PW91 codes, May 14, 1996
//  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  INPUT rho : DENSITY
//  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
//  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
//  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
//   (for U,V, see PW86(24))
//  input lgga:  (=0=>don't put in gradient corrections, just LDA) ?????
//  input lpot:  (=0=>don't get potential and don't need U and V) ??????
//  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
//----------------------------------------------------------------------
// References:
// [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
// [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986); {\bf 40},
// 3399  (1989) (E).
//----------------------------------------------------------------------
// Formulas:
//      e_x[unif]=ax*rho^(4/3)  [LDA]
//
// AX = -0.75*(3/pi)^(1/3)
//
//       e_x[PBE]=e_x[unif]*FxPBE(s)
//       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](14)
//
// UK, UL defined after [a](14)
//----------------------------------------------------------------------
#include "../fix_compile.h"

namespace G2G {
#define EASYPBE_AX ((scalar_type)-0.738558766382022405884230032680836f)
#define EASYPBE_UM ((scalar_type)0.2195149727645171f)
#define EASYPBE_UK ((scalar_type)0.804f)
#define EASYPBE_UL ((scalar_type)0.273028573090195f)
// UL = UM / UK

template <class scalar_type>
__host__ __device__ void pbeOS_exch(scalar_type rho, scalar_type s,
                                    scalar_type u, scalar_type v,
                                    scalar_type& ex, scalar_type& vx) {
  // Construct LDA exchange energy density
  scalar_type rho13 = cbrt(rho);
  scalar_type exunif = (scalar_type)EASYPBE_AX * rho13;

  // Construct PBE enhancement factor and calculate Ex
  scalar_type s2 = s * s;
  scalar_type s3 = s * s2;
  scalar_type p0 = (scalar_type)1.0f + EASYPBE_UL * s2;

  scalar_type fxpbe = (scalar_type)1.0f + EASYPBE_UK - EASYPBE_UK / p0;
  ex = exunif * fxpbe;

  // Find first and second derivatives of Fx w.r.t s.
  // Fs = (1/s)*dFxPBE/ds  || Fss = d Fs/ds
  scalar_type p2 = p0 * p0;
  scalar_type Fs = 2.0f * EASYPBE_UM / p2;

  scalar_type F1 = -4.0f * EASYPBE_UL * s * Fs;
  scalar_type Fss = F1 / p0;

  // Calculate potential from [b](24)
  scalar_type vx2 = (4.0f / 3.0f) * fxpbe;
  scalar_type vx3 = v * Fs;
  scalar_type vx4 = (u - (4.0f / 3.0f) * s3) * Fss;
  vx = exunif * (vx2 - vx4 - vx3);
}
#undef EASYPBE_AX
#undef EASYPBE_UM
#undef EASYPBE_UK
#undef EASYPBE_UL
}
