#include <cassert>
#include <stdexcept>
#include "../fix_compile.h"

namespace G2G {

// POT_ALPHA = -(3/PI)^(1/3)
#define POT_ALPHA ((scalar_type)-0.738558766382022447)
#define POT_GL ((scalar_type)0.620350490899400087)
#define POT_VOSKO_A1 ((scalar_type)0.03109205)
#define POT_VOSKO_B1 ((scalar_type)3.72744)
#define POT_VOSKO_C1 ((scalar_type)12.9352)
#define POT_VOSKO_X0 ((scalar_type)-0.10498)
#define POT_VOSKO_Q ((scalar_type)6.15199066246304849)
#define POT_VOSKO_A16 ((scalar_type)0.005182008333)

// POT_T4          = POT_VOSKO_B1 * POT_VOSKO_X0 / POT_XX0
// POT_VOSKO_2C1   = 2 * POT_VOSKO_C1
// POT_VOSKO_2B1Q  = 2 * POT_VOSKO_B1 / POT_VOSKO_Q
// POT_VOSKO_B2X0Q = 2 * (POT_VOSKO_B1 + 2 * POT_VOSKO_X0) / POT_VOSKO_Q
// POT_VOSKO_4B1   = 4.0 * POT_VOSKO_B1
// POT_VOSKO_QSQ   = POT_VOSKO_Q * POT_VOSKO_Q
// POT_VOSKO_B1X0  = (1.0f - t6 * (POT_VOSKO_B1 - 2.0f * POT_VOSKO_X0))
#define POT_T4 ((scalar_type)-0.0311676086789438)
#define POT_VOSKO_2C1 ((scalar_type)25.8704)
#define POT_VOSKO_2B1Q ((scalar_type)1.21178337371132)
#define POT_VOSKO_B2X0Q ((scalar_type)1.14352579286644)
#define POT_VOSKO_4B1 ((scalar_type)14.90976)
#define POT_VOSKO_QSQ ((scalar_type)37.8469891110325)
#define POT_VOSKO_B1X0 ((scalar_type)1.0329232240928)

template <class scalar_type, int iexch>
__host__ __device__ void calc_ldaCS(scalar_type dens, scalar_type& ex,
                                    scalar_type& ec, scalar_type& y2a) {
  // data X alpha

  if (dens == 0) {
    ex = 0;
    ec = 0;
    y2a = 0;
    return;
  }

  scalar_type y = cbrt(dens);  // rho^(1/3)
  scalar_type v0 =
      (scalar_type)-0.984745021842697 * y;  // -4/3 * (3/PI)^(1/3) * rho^(1/3)
  ex = POT_ALPHA * y;                       // -(3/PI)^(1/3) * rho^(1/3)

  switch (iexch) {
    case 1: {
      ec = 0;
      y2a = v0;
    } break;

    case 2: {
      scalar_type rs = POT_GL / y;
      scalar_type x1 = rs / (scalar_type)11.4;
      scalar_type vc;

      if (x1 > 1.0) {
        ec = (scalar_type)-0.0333 *
             ((scalar_type)0.5 * x1 - (scalar_type)0.33333333333333);
        vc = (scalar_type)0.0111 * x1 * (scalar_type)0.5;
      } else {
        scalar_type t1 = ((scalar_type)1.0 + x1 * x1 * x1);
        scalar_type t2 = logf((scalar_type)1.0 + (scalar_type)1.0 / x1);
        scalar_type t3 = x1 * x1;
        ec = (scalar_type)-0.0333 * (t1 * t2 - t3 + (scalar_type)0.5 * x1 -
                                     (scalar_type)0.33333333333333);
        vc = (scalar_type)0.0111 * x1 *
             ((scalar_type)3.0 * t3 * t2 - t1 / (x1 * (x1 + (scalar_type)1.0)) -
              (scalar_type)2.0 * x1 + (scalar_type)0.5);
      }
      y2a = v0 + ec + vc;
    } break;

    case 3: {
      scalar_type rs = POT_GL / y;
      scalar_type x1 = sqrt(rs);
      scalar_type Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
      scalar_type t1 = (scalar_type)2.0 * x1 + POT_VOSKO_B1;
      scalar_type t2 = log(Xx);
      scalar_type t3 = atan(POT_VOSKO_Q / t1);
      scalar_type t5 = (POT_VOSKO_B1 * x1 + POT_VOSKO_2C1) / x1;

      ec = ((scalar_type)2.0 * logf(x1) - t2 + POT_VOSKO_2B1Q * t3 -
            POT_T4 * ((scalar_type)2.0 * log(x1 - POT_VOSKO_X0) - t2 +
                      POT_VOSKO_B2X0Q * t3)) *
           POT_VOSKO_A1;

      scalar_type vc;
      vc = ec -
           POT_VOSKO_A16 * x1 *
               (t5 / Xx -
                POT_VOSKO_4B1 / (t1 * t1 + POT_VOSKO_QSQ) * POT_VOSKO_B1X0 -
                POT_T4 * ((scalar_type)2.0 / (x1 - POT_VOSKO_X0) - t1 / Xx));
      y2a = v0 + vc;
    } break;
    default:
      throw std::runtime_error("Error, invalid iexch");
  }
}

template <class scalar_type>
__host__ __device__ void calc_ldaCS_in(scalar_type dens, scalar_type& ex,
                                       scalar_type& ec, scalar_type& y2a,
                                       int iexch) {
  switch (iexch) {
    case 0:
      return calc_ldaCS<scalar_type, 0>(dens, ex, ec, y2a);
    case 1:
      return calc_ldaCS<scalar_type, 1>(dens, ex, ec, y2a);
    case 2:
      return calc_ldaCS<scalar_type, 2>(dens, ex, ec, y2a);
    case 3:
      return calc_ldaCS<scalar_type, 3>(dens, ex, ec, y2a);
    case 4:
      return calc_ldaCS<scalar_type, 4>(dens, ex, ec, y2a);
    case 5:
      return calc_ldaCS<scalar_type, 5>(dens, ex, ec, y2a);
    case 6:
      return calc_ldaCS<scalar_type, 6>(dens, ex, ec, y2a);
    case 7:
      return calc_ldaCS<scalar_type, 7>(dens, ex, ec, y2a);
    case 8:
      return calc_ldaCS<scalar_type, 8>(dens, ex, ec, y2a);
    case 9:
      return calc_ldaCS<scalar_type, 9>(dens, ex, ec, y2a);
    default:
      assert(false);
  }
}

#undef POT_ALPHA
#undef POT_GL
#undef POT_VOSKO_A1
#undef POT_VOSKO_B1
#undef POT_VOSKO_C1
#undef POT_VOSKO_X0
#undef POT_VOSKO_Q
#undef POT_VOSKO_A16
#undef POT_T4
#undef POT_VOSKO_2C1
#undef POT_VOSKO_2B1Q
#undef POT_VOSKO_B2X0Q
#undef POT_VOSKO_4B1
#undef POT_VOSKO_QSQ
#undef POT_VOSKO_B1X0
}
