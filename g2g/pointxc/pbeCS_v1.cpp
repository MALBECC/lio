#include "pbeCS_v1.h"
//#include <limits>
//#include "../common.h"
//#include "../init.h"
//#include "../cuda/cuda_extra.h"
//#include "../matrix.h"
//#include <float.h>
//#include "gcorc.h"


namespace pointxc {

#define CLOSEDPBE_PI32     ((scalar_type) 29.608813203268075856503472999628)
#define CLOSEDPBE_AX       ((scalar_type)-0.738558766382022405884230032680836)
#define CLOSEDPBE_UM       ((scalar_type) 0.2195149727645171)
#define CLOSEDPBE_UK       ((scalar_type) 0.804)
#define CLOSEDPBE_UL       ((scalar_type) 0.273028573090195) // um / uk
#define CLOSEDPBE_GAMMA    ((scalar_type) 0.03109069086965489503494086371273)
#define CLOSEDPBE_GAMMAINV ((scalar_type) 32.1639684429148) // 1 / gamma
#define CLOSEDPBE_BETA     ((scalar_type) 0.06672455060314922)
#define CLOSEDPBE_DELTA    ((scalar_type) 2.14612633996736) // beta/gamma

template<class scalar_type>
__host__ __device__ 
static void pbeCS_v1( scalar_type rho,
                      scalar_type agrad,
                      scalar_type delgrad,
                      scalar_type rlap,
                      scalar_type& expbe,
                      scalar_type& vxpbe,
                      scalar_type& ecpbe,
                      scalar_type& vcpbe ){


   if (rho < 2e-18) {
      expbe = vxpbe = ecpbe = vcpbe = 0;
      return;
   }

   scalar_type rho2  = rho * rho;
   scalar_type rho13 = cbrt( (scalar_type) rho );
   scalar_type fk1   = cbrt( (scalar_type) CLOSEDPBE_PI32 );
   scalar_type fk    = fk1 * rho13;

   scalar_type twofk  = 2.0f * fk;
   scalar_type twofk2 = twofk * twofk;
   scalar_type twofk3 = twofk * twofk2;

   // S = |grad(rho)|/(2*fk*rho)
   scalar_type s  = agrad / (twofk * rho);
   scalar_type s2 = s * s;
   scalar_type s3 = s * s2;

   // LDA exchange contribution:
   // ex*rho ==> energy, we will calculate ex ==> energy density
   // ex*rho = -(3/4Pi)*(e^2)*(3pi)^2/3*rho^1/3*rho
   // ex*rho = -0.75*(3/Pi)^1/3*rho^4/3
   // ex*rho = ax*rho^4/3
   scalar_type exlda = CLOSEDPBE_AX * rho13;

   // In order to calculate the PBE contribution
   // to exchange energy, we have to calculate the
   // enhancement function Fx:
   // Fx = 1+uk -(uk/(1+(um*s^2)/uk)
   // um/uk = ul
   // P0 = 1 + (um*s^2)/uk
   scalar_type p0    = 1.0f + CLOSEDPBE_UL * s2;
   scalar_type fxpbe = 1.0f + CLOSEDPBE_UK - CLOSEDPBE_UK/p0;

   // exchange pbe energy
   expbe = exlda * fxpbe;

   // Now the potential:
   scalar_type v = rlap / (twofk2 * rho);
   scalar_type u = (delgrad == 0 ? 0 : delgrad / (twofk3 * rho2));

   // Calculation of first and second derivatives
   scalar_type P2  = p0 * p0;
   scalar_type Fs  =  2.0f * CLOSEDPBE_UM / P2;
   scalar_type F1  = -4.0f * CLOSEDPBE_UL * s * Fs;
   scalar_type Fss = F1/p0;

   // Now we calculate the potential Vx
   scalar_type vx2 = (4.0f / 3.0f) * fxpbe;
   scalar_type vx3 = v * Fs;
   scalar_type vx4 = (u - (4.0f / 3.0f) * s3) * Fss;
   vxpbe = exlda * (vx2 - vx4 - vx3);

   // Now we need to calculate the Correlation contribution
   // to the energy
   // ecpbe = eclsd*rho + h*rho
   // first we calculate the lsd contribution to the correlation energy
   // we will use the subroutine GCOR.
   // We need only the  rs (seitz radius) rs = (3/4pi*rho)^1/3
   scalar_type pirho = 4.0f * M_PI * rho;
   scalar_type rs    = cbrt( 3.0f / pirho );
   scalar_type rtrs  = sqrt(rs);

   scalar_type sk    = sqrt(4.0f * fk / M_PI);
   scalar_type twoks = 2.0f * sk;

   scalar_type t  = agrad / (twoks * rho);
   scalar_type t2 = t * t;

   scalar_type twoks2 = twoks * twoks;
   scalar_type twoks3 = twoks * twoks2;

   scalar_type UU = (delgrad == 0 ? 0 : delgrad / (rho2 * twoks3));
   scalar_type VV = rlap / (rho * twoks2);

   scalar_type ec, eurs;
   gcorc1(rtrs, ec, eurs);

   if (ec == 0) ec = numeric_limits<scalar_type>::min();
   scalar_type eclda = ec;
   scalar_type ecrs  = eurs;
   scalar_type vclda = eclda - rs * (1.0f / 3.0f) * ecrs;

   // Now we have to calculate the H function in order to evaluate
   // the GGA contribution to the correlation energy
   scalar_type PON = -ec * CLOSEDPBE_GAMMAINV;
   scalar_type B   = CLOSEDPBE_DELTA / (exp(PON) - 1.0f);
   scalar_type B2  = B * B;
   scalar_type T4  = t2 * t2;

   scalar_type Q4 = 1.0f + B * t2;
   scalar_type Q5 = 1.0f + B * t2 + B2 * T4;

   scalar_type H = (CLOSEDPBE_BETA/CLOSEDPBE_DELTA) * log(1.0f + CLOSEDPBE_DELTA * Q4 * t2/Q5);

   // So the correlation energy for pbe is:
   ecpbe = eclda + H;
   //cout << expl(PON) << " " << t2 << endl;

   // Now we have to calculate the potential contribution of GGA
   scalar_type T6 = T4 * t2;
   scalar_type RSTHRD = rs / 3.0f;
   scalar_type FAC = CLOSEDPBE_DELTA / B + 1.0f;
   scalar_type BEC = B2 * FAC / CLOSEDPBE_BETA;
   scalar_type Q8 = Q5 * Q5 + CLOSEDPBE_DELTA * Q4 * Q5 * t2;
   scalar_type Q9 = 1.0f + 2.0f * B * t2;
   scalar_type hB = -CLOSEDPBE_BETA * B * T6 * (2.0f + B * t2)/Q8;
   scalar_type hRS = -RSTHRD * hB * BEC * ecrs;
   scalar_type FACT0 = 2.0f * CLOSEDPBE_DELTA - 6.0f * B;
   scalar_type FACT1 = Q5 * Q9 + Q4 * Q9 * Q9;
   scalar_type hBT = 2.0f * CLOSEDPBE_BETA * T4 * ((Q4 * Q5 * FACT0 - CLOSEDPBE_DELTA * FACT1)/Q8)/Q8;
   scalar_type hRST = RSTHRD * t2 * hBT * BEC * ecrs;
   scalar_type hT = 2.0f * CLOSEDPBE_BETA * Q9/Q8;
   scalar_type FACT2 = Q4 * Q5 + B * t2 * (Q4 * Q9 + Q5);
   scalar_type FACT3 = 2.0f * B * Q5 * Q9 + CLOSEDPBE_DELTA * FACT2;
   scalar_type hTT = 4.0f * CLOSEDPBE_BETA * t * (2.0f * B/Q8 -(Q9 * FACT3 / Q8)/Q8);
   scalar_type COMM = H + hRS + hRST + t2 * hT/6.0f + 7.0f * t2 * t * hTT/6.0f;

   COMM = COMM - UU * hTT - VV * hT;

   // Then, the potential for PBE is:
   vcpbe = vclda + COMM;
}

}
