//----------------------------------------------------------------------
//      SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
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
// [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986); {\bf 40},  3399  (1989) (E).
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
namespace pointxc {

template<class scalar_type> __host__ __device__
void pbeOS_exch( scalar_type rho, scalar_type s, scalar_type u, scalar_type v,
                 scalar_type& ex, scalar_type& vx);

}
