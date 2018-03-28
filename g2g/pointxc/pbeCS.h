#include <float.h>
#include <math.h>

#include "gcorc.h"
#include "../fix_compile.h"

namespace G2G {

#define CLOSEDPBE_PI32 ((scalar_type)29.608813203268075856503472999628f)
#define CLOSEDPBE_AX ((scalar_type)-0.738558766382022405884230032680836f)
#define CLOSEDPBE_UM ((scalar_type)0.2195149727645171f)
#define CLOSEDPBE_UK ((scalar_type)0.804f)
#define CLOSEDPBE_UL ((scalar_type)0.273028573090195f)
#define CLOSEDPBE_GAMMA ((scalar_type)0.03109069086965489503494086371273f)
#define CLOSEDPBE_GAMMAINV ((scalar_type)32.1639684429148f)
#define CLOSEDPBE_BETA ((scalar_type)0.06672455060314922f)
#define CLOSEDPBE_DELTA ((scalar_type)2.14612633996736f)

/*---------------------------------------------------------//
// pbeCS inputs:
//   rho:     density.
//   agrad:   module of the gradient.
//   delgrad: grad . grad(agrad)
//   rlap:    laplacian.
// pbeCS outputs:
//   expbe, ecpbe: energies of exchange and correlation.
//   vxpbe, vcpbe: potentials of exchange and correlation.
//---------------------------------------------------------*/

template <class scalar_type>
__host__ __device__ void pbeCS(scalar_type rho, scalar_type agrad,
                               scalar_type delgrad, scalar_type rlap,
                               scalar_type& expbe, scalar_type& vxpbe,
                               scalar_type& ecpbe, scalar_type& vcpbe) {

  scalar_type rho2 = rho * rho;
  scalar_type rho5 = rho2*rho2*rho;
  scalar_type flt_minimum = 100.0 * (scalar_type)FLT_MIN;
  if (rho5 < flt_minimum) {
    expbe = (scalar_type)0.0f;
    vxpbe = (scalar_type)0.0f;
    ecpbe = (scalar_type)0.0f;
    vcpbe = (scalar_type)0.0f;
    return;
  }
  /*----------------------------------//
  // PBE exchange
  //----------------------------------//
  // Local Fermi wavevector for 2*up
  // fk = (3 pi^2 (rho))^(1/3)
  //
  // Dimensionless density gradient
  // s = |grad rho| / (2*fk*rho)
  //
  // Other magnitudes
  // u = delgrad /(rho^2 * (2*fk)^3)
  // v = rlap    /(rho   * (2*fk)^2)
  //---------------------------------*/

  scalar_type rho13 = cbrt(rho);
  scalar_type fk1 = cbrt((scalar_type)CLOSEDPBE_PI32);
  scalar_type fk = fk1 * rho13;

  scalar_type twofk = 2.0f * fk;
  scalar_type twofk2 = twofk * twofk;
  scalar_type twofk3 = twofk * twofk2;

  scalar_type s = agrad / (twofk * rho);
  scalar_type s2 = s * s;
  scalar_type s3 = s * s2;

  // LDA exchange contribution:
  scalar_type exlda = CLOSEDPBE_AX * rho13;

  // PBE enhancement function and exchange energy:
  scalar_type p0 = 1.0f + CLOSEDPBE_UL * s2;
  scalar_type fxpbe = 1.0f + CLOSEDPBE_UK - CLOSEDPBE_UK / p0;

  expbe = exlda * fxpbe;

  // Now the potentials:
  scalar_type v = rlap / (twofk2 * rho);
  scalar_type u = (delgrad == 0.0f ? 0.0f : delgrad / (twofk3 * rho2));

  // First and second derivatives.
  scalar_type P2 = p0 * p0;
  scalar_type Fs = 2.0f * CLOSEDPBE_UM / P2;
  scalar_type F1 = -4.0f * CLOSEDPBE_UL * s * Fs;
  scalar_type Fss = F1 / p0;

  // Potential Vx
  scalar_type vx2 = (4.0f / 3.0f) * fxpbe;
  scalar_type vx3 = v * Fs;
  scalar_type vx4 = (u - (4.0f / 3.0f) * s3) * Fss;
  vxpbe = exlda * (vx2 - vx4 - vx3);

  /*-----------------------------------------------//
  // PBE correlation
  //-----------------------------------------------//
  // Local Seitz radius (alpha/fk)
  // rs = (3/(4pi*rho))^(1/3)
  //
  // Thomas-Fermi screening wavevector
  // sk = Ks = sqrt(4fk/pi) || twoksg = 2*Ks*phi
  //
  // Correlation dimensionless gradient
  // t = |grad rho| / (2*Ks*phi*rho)
  //
  // uu = delgrad / (rho^2 * twoksg^3)
  // vv = rlap    / (rho   * twoksg^2)
  //
  // ec    = LSD correlation energy
  // vclda = LSD correlation potential
  // h     = Gradient correlation energy
  // COMM  = Gradient correlation potential
  //-----------------------------------------------*/

  // LSD contribution to correlation energy.
  scalar_type pirho = 4.0f * (scalar_type)M_PI * rho;
  scalar_type rs = cbrt(3.0f / pirho);
  scalar_type rtrs = sqrt(rs);
  scalar_type sk = sqrt(4.0f * fk / (scalar_type)M_PI);

  scalar_type twoks = 2.0f * sk;
  scalar_type twoks2 = twoks * twoks;
  scalar_type twoks3 = twoks2 * twoks;

  scalar_type t = agrad / (twoks * rho);
  scalar_type t2 = t * t;

  scalar_type UU = (delgrad == 0.0f ? 0.0f : delgrad / (rho2 * twoks3));
  scalar_type VV = rlap / (rho * twoks2);

  scalar_type ec, eurs;
  gcorc1(rtrs, ec, eurs);
  if (ec == (scalar_type)0.0f) ec = (scalar_type)FLT_MIN;
  scalar_type eclda = ec;
  scalar_type ecrs = eurs;
  scalar_type vclda = eclda - rs * (1.0f / 3.0f) * ecrs;

  // H function to evaluate the GGA contribution to the correlation energy.
  scalar_type PON = -ec * CLOSEDPBE_GAMMAINV;
  scalar_type B = CLOSEDPBE_DELTA / (exp(PON) - 1.0f);
  scalar_type B2 = B * B;
  scalar_type T4 = t2 * t2;
  scalar_type Q4 = 1.0f + B * t2;
  scalar_type Q5 = 1.0f + B * t2 + B2 * T4;
  scalar_type H = (CLOSEDPBE_BETA / CLOSEDPBE_DELTA) *
                  log(1.0f + CLOSEDPBE_DELTA * Q4 * t2 / Q5);

  ecpbe = eclda + H;

  // GGA Contribution to the potential.
  scalar_type T6 = T4 * t2;
  scalar_type RSTHRD = rs / 3.0f;
  scalar_type FAC = CLOSEDPBE_DELTA / B + 1.0f;
  scalar_type BEC = B2 * FAC / CLOSEDPBE_BETA;
  scalar_type Q8 = Q5 * Q5 + CLOSEDPBE_DELTA * Q4 * Q5 * t2;
  scalar_type Q9 = 1.0f + 2.0f * B * t2;
  scalar_type hB = -CLOSEDPBE_BETA * B * T6 * (2.0f + B * t2) / Q8;
  scalar_type hRS = -RSTHRD * hB * BEC * ecrs;
  scalar_type FACT0 = 2.0f * CLOSEDPBE_DELTA - 6.0f * B;
  scalar_type FACT1 = Q5 * Q9 + Q4 * Q9 * Q9;
  scalar_type hBT = 2.0f * CLOSEDPBE_BETA * T4 *
                    ((Q4 * Q5 * FACT0 - CLOSEDPBE_DELTA * FACT1) / Q8) / Q8;
  scalar_type hRST = RSTHRD * t2 * hBT * BEC * ecrs;
  scalar_type hT = 2.0f * CLOSEDPBE_BETA * Q9 / Q8;
  scalar_type FACT2 = Q4 * Q5 + B * t2 * (Q4 * Q9 + Q5);
  scalar_type FACT3 = 2.0f * B * Q5 * Q9 + CLOSEDPBE_DELTA * FACT2;
  scalar_type hTT =
      4.0f * CLOSEDPBE_BETA * t * (2.0f * B / Q8 - (Q9 * FACT3 / Q8) / Q8);
  scalar_type COMM =
      H + hRS + hRST + t2 * hT / 6.0f + 7.0f * t2 * t * hTT / 6.0f;

  COMM = COMM - UU * hTT - VV * hT;
  vcpbe = vclda + COMM;
}  // pbeCS
}
