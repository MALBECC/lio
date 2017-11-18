//----------------------------------------------------------------------------//
// SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
//
// slimmed down version of GCOR used in PW91 routines, to interpolate
// LSD correlation energy, as given by (10) of
// J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
// K. Burke, May 11, 1996.
//----------------------------------------------------------------------------//
#include <float.h>
#include <math.h>
#include "../fix_compile.h"

namespace G2G {

#define GCORC_A0_1 ((scalar_type)0.0310907f)
#define GCORC_A1_1 ((scalar_type)0.21370f)
#define GCORC_B1_1 ((scalar_type)7.5957f)
#define GCORC_B2_1 ((scalar_type)3.5876f)
#define GCORC_B3_1 ((scalar_type)1.6382f)
#define GCORC_B4_1 ((scalar_type)0.49294f)
template <class scalar_type>
__host__ __device__ void gcorc1(scalar_type rtrs, scalar_type& gg,
                                scalar_type& grrs) {
  scalar_type Q0 = -2.0f * GCORC_A0_1 * (1.0f + GCORC_A1_1 * rtrs * rtrs);
  scalar_type Q1 =
      2.0f * GCORC_A0_1 * rtrs *
      (GCORC_B1_1 +
       rtrs * (GCORC_B2_1 + rtrs * (GCORC_B3_1 + GCORC_B4_1 * rtrs)));
  scalar_type Q2 = logf(1.0f + 1.0f / Q1);
  gg = Q0 * Q2;
  scalar_type Q3 =
      GCORC_A0_1 * (GCORC_B1_1 / rtrs + 2.0f * GCORC_B2_1 +
                    rtrs * (3.0f * GCORC_B3_1 + 4.0f * GCORC_B4_1 * rtrs));
  grrs = -2.0f * GCORC_A0_1 * GCORC_A1_1 * Q2 - Q0 * Q3 / (Q1 * (1.0f + Q1));
}  // gcorc1

#define GCORC_A0_2 ((scalar_type)0.01554535f)
#define GCORC_A1_2 ((scalar_type)0.20548f)
#define GCORC_B1_2 ((scalar_type)14.1189f)
#define GCORC_B2_2 ((scalar_type)6.1977f)
#define GCORC_B3_2 ((scalar_type)3.3662f)
#define GCORC_B4_2 ((scalar_type)0.62517f)
template <class scalar_type>
__host__ __device__ void gcorc2(scalar_type rtrs, scalar_type& gg,
                                scalar_type& grrs) {
  scalar_type Q0 = -2.0f * GCORC_A0_2 * (1.0f + GCORC_A1_2 * rtrs * rtrs);
  scalar_type Q1 =
      2.0f * GCORC_A0_2 * rtrs *
      (GCORC_B1_2 +
       rtrs * (GCORC_B2_2 + rtrs * (GCORC_B3_2 + GCORC_B4_2 * rtrs)));
  scalar_type Q2 = logf(1.0f + 1.0f / Q1);
  gg = Q0 * Q2;
  scalar_type Q3 =
      GCORC_A0_2 * (GCORC_B1_2 / rtrs + 2.0f * GCORC_B2_2 +
                    rtrs * (3.0f * GCORC_B3_2 + 4.0f * GCORC_B4_2 * rtrs));
  grrs = -2.0f * GCORC_A0_2 * GCORC_A1_2 * Q2 - Q0 * Q3 / (Q1 * (1.0f + Q1));
}  // gcorc2

#define GCORC_A0_3 ((scalar_type)0.0168869f)
#define GCORC_A1_3 ((scalar_type)0.11125f)
#define GCORC_B1_3 ((scalar_type)10.357f)
#define GCORC_B2_3 ((scalar_type)3.6231f)
#define GCORC_B3_3 ((scalar_type)0.88026f)
#define GCORC_B4_3 ((scalar_type)0.49671f)
template <class scalar_type>
__host__ __device__ void gcorc3(scalar_type rtrs, scalar_type& gg,
                                scalar_type& grrs) {
  scalar_type Q0 = -2.0f * GCORC_A0_3 * (1.0f + GCORC_A1_3 * rtrs * rtrs);
  scalar_type Q1 =
      2.0f * GCORC_A0_3 * rtrs *
      (GCORC_B1_3 +
       rtrs * (GCORC_B2_3 + rtrs * (GCORC_B3_3 + GCORC_B4_3 * rtrs)));
  scalar_type Q2 = logf(1.0f + 1.0f / Q1);
  gg = Q0 * Q2;
  scalar_type Q3 =
      GCORC_A0_3 * (GCORC_B1_3 / rtrs + 2.0f * GCORC_B2_3 +
                    rtrs * (3.0f * GCORC_B3_3 + 4.0f * GCORC_B4_3 * rtrs));
  grrs = -2.0f * GCORC_A0_3 * GCORC_A1_3 * Q2 - Q0 * Q3 / (Q1 * (1.0f + Q1));
}  // gcorc3
}
