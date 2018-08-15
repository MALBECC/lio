//----------------------------------------------------------------------
//      SUBROUTINE
//      CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,H,DVCUP,DVCDN,rho)
//----------------------------------------------------------------------
//  Official PBE correlation code. K. Burke, May 14, 1996.
//  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
//       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
//       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
//       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
//       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
//       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
//       :  UU,VV,WW, only needed for PBE potential
//       : lgga=flag to do gga (0=>LSD only)  ???
//       : lpot=flag to do potential (0=>energy only) ??
//  output: ec=lsd correlation energy from [a]
//        : vcup=lsd up correlation potential
//        : vcdn=lsd dn correlation potential
//        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
//        : dvcup=nonlocal correction to vcup
//        : dvcdn=nonlocal correction to vcdn
//----------------------------------------------------------------------
// References:
// [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
//     {\sl Generalized gradient approximation made simple}, sub.
//     to Phys. Rev.Lett. May 1996.
// [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
//     construction of a generalized gradient approximation:  The PW91
//     density functional}, submitted to Phys. Rev. B, Feb. 1996.
// [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
//----------------------------------------------------------------------
// thrd*=various multiples of 1/3
// numbers for use in LSD energy spin-interpolation formula, [c](9).
//      GAM= 2^(4/3)-2
//      FZZ=f''(0)= 8/(9*GAM)
// numbers for construction of PBE
//      gama=(1-log(2))/pi^2
//      bet=coefficient in gradient expansion for correlation, [a](4).
//      eta=small number to stop d phi/ dzeta from blowing up at
//          |zeta|=1.
//
//      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
//      parameter(sixthm=thrdm/2.d0)
//      parameter(thrd4=4.d0*thrd)
//      parameter(GAM=0.5198420997897463295344212145565d0)
//      parameter(fzz=8.d0/(9.d0*GAM))
//      parameter(gama=0.03109069086965489503494086371273d0)
//      parameter(bet=0.06672455060314922d0,delt=bet/gama)
//      parameter(eta=1.d-12)
//----------------------------------------------------------------------
#include <float.h>
#include "../fix_compile.h"

namespace G2G {

#define EASYPBE_GAMMA \
  ((scalar_type)0.03109069086965489503494086371273f)      // THIS
#define EASYPBE_BETA ((scalar_type)0.06672455060314922f)  // THIS
#define EASYPBE_DELTA \
  ((scalar_type)(EASYPBE_BETA / EASYPBE_GAMMA))     // beta/gamma // THIS
#define EASYPBE_ETA ((scalar_type)0.000000000001f)  // THIS

template <class scalar_type>
__host__ __device__ void pbeOS_corr(scalar_type rho, scalar_type rs,
                                    scalar_type zet, scalar_type t,
                                    scalar_type uu, scalar_type vv,
                                    scalar_type ww, scalar_type& ec,
                                    scalar_type& vc_a, scalar_type& vc_b,
                                    scalar_type& h, scalar_type& dvc_a,
                                    scalar_type& dvc_b) {
  // find LSD energy contributions, using [c](10) and Table I[c].
  // EU=unpolarized LSD correlation energy
  // EURS=dEU/drs
  // EP=fully polarized LSD correlation energy
  // EPRS=dEP/drs
  // ALFM=-spin stiffness, [c](3).
  // ALFRSM=-dalpha/drs
  // F=spin-scaling factor from [c](9).
  // construct ec, using [c](8)

  scalar_type rtrs = sqrt(rs);
  scalar_type EU, EURS, EP, EPRS, ALFM, ALFRSM;
  gcorc1(rtrs, EU, EURS);
  gcorc2(rtrs, EP, EPRS);
  gcorc3(rtrs, ALFM, ALFRSM);

  //  scalar_type ALFC = -ALFM;
  scalar_type Z2 = zet * zet;
  scalar_type Z4 = Z2 * Z2;
  scalar_type ZET1 = cbrt(pow((scalar_type)(1.0f + zet), 4));
  scalar_type ZET2 = cbrt(pow((scalar_type)(1.0f - zet), 4));
  scalar_type GAM =
      pow((scalar_type)(2.0f), (scalar_type)4.0f / 3.0f) - (scalar_type)2.0f;

  // F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
  // EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
  scalar_type F = (ZET1 + ZET2 - (scalar_type)2.0f) / GAM;
  scalar_type FZZ = (scalar_type)8.0f / (9.0f * GAM);
  ec = EU * (1.0f - F * Z4) + EP * F * Z4 - ALFM * F * (1.0f - Z4) / FZZ;

  //----------------------------------------------------------------------
  // LSD potential from [c](A1)
  // ECRS = dEc/drs [c](A2)
  // ECZET=dEc/dzeta [c](A3)
  // FZ = dF/dzeta [c](A4)
  //
  // ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
  // FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM

  scalar_type ECRS =
      EURS * (1.0f - F * Z4) + EPRS * F * Z4 - ALFRSM * F * (1.0f - Z4) / FZZ;
  scalar_type ZET3 = cbrt((scalar_type)(1.0f + zet));
  scalar_type ZET4 = cbrt((scalar_type)(1.0f - zet));
  scalar_type FZ = (scalar_type)(4.0f / 3.0f) * (ZET3 - ZET4) / GAM;

  // ECZET =
  // 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU-(1.D0-Z4)*ALFM/FZZ)
  // COMM = EC -RS*ECRS/3.D0-ZET*ECZET
  scalar_type ECZET = (scalar_type)4.0f * pow((scalar_type)zet, 3) * F *
                          (EP - EU + ALFM / FZZ) +
                      FZ * (Z4 * EP - Z4 * EU - (1.0f - Z4) * ALFM / FZZ);
  scalar_type COMM = ec - rs * ECRS / (scalar_type)3.0f - zet * ECZET;
  vc_a = COMM + ECZET;
  vc_b = COMM - ECZET;

  //----------------------------------------------------------------------
  // PBE correlation energy
  //
  //     if(lgga.eq.0)return  // WHAT IS THIS?
  //----------------------------------------------------------------------
  // G=phi(zeta), given after [a](3)
  // DELT=bet/gamma
  // B=A of [a](8)
  //
  // G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
  // G3 = G**3
  // PON=-EC/(G3*gama)
  scalar_type ZET5 = cbrt(pow((scalar_type)(1.0f + zet), 2));
  scalar_type ZET6 = cbrt(pow((scalar_type)(1.0f - zet), 2));
  scalar_type G = (ZET5 + ZET6) / (scalar_type)2.0f;
  scalar_type G3 = pow(G, 3);
  scalar_type PON = -ec / (G3 * EASYPBE_GAMMA);

  // B = DELT/(DEXP(PON)-1.D0)
  // B2 = B*B
  // T2 = T*T
  // T4 = T2*T2

  scalar_type B1 = exp(PON);
  scalar_type B = EASYPBE_DELTA / (B1 - (scalar_type)1.0f);
  scalar_type B2 = B * B;
  scalar_type T2 = t * t;
  scalar_type T4 = T2 * T2;

  // RS2 = RS*RS
  // RS3 = RS2*RS
  // Q4 = 1.D0+B*T2
  // Q5 = 1.D0+B*T2+B2*T4
  // H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
  scalar_type RS2 = pow(rs, 2);
  scalar_type Q4 = (scalar_type)1.0f + B * T2;
  scalar_type Q5 = (scalar_type)1.0f + B * T2 + B2 * T4;
  h = G3 * (scalar_type)(EASYPBE_BETA / EASYPBE_DELTA) *
      log((scalar_type)1.0f + EASYPBE_DELTA * Q4 * T2 / Q5);

  // ENERGY DONE !!!!!!!
  //=========================================================
  // NOW THE POTENTIAL, using appendix E of [b].
  //
  //  if(lpot.eq.0)return  // WHAT IS THIS???
  //
  // G4 = G3*G
  // T6 = T4*T2
  // RSTHRD = RS/3.D0
  scalar_type G4 = G3 * G;
  scalar_type T6 = T4 * T2;
  scalar_type rsthrd = rs / (scalar_type)3.0f;

  // GZ=(((1.d0+zet)**2+eta)**sixthm-((1.d0-zet)**2+eta)**sixthm)/3.d0
  scalar_type ZET7 = pow((scalar_type)(1.0f + zet), 2);
  scalar_type ZET8 = pow((scalar_type)(1.0f - zet), 2);
  scalar_type GZ1 = pow((ZET7 + EASYPBE_ETA), (scalar_type)(-1.0f / 6.0f));
  scalar_type GZ2 = pow((ZET8 + EASYPBE_ETA), (scalar_type)(-1.0f / 6.0f));
  scalar_type GZ = (GZ1 - GZ2) / (scalar_type)3.0f;

  //  FAC = DELT/B+1.D0
  //  BG = -3.D0*B2*EC*FAC/(BET*G4)
  //  BEC = B2*FAC/(BET*G3)
  //  Q8 = Q5*Q5+DELT*Q4*Q5*T2
  //  Q9 = 1.D0+2.D0*B*T2
  scalar_type FAC = EASYPBE_DELTA / B + (scalar_type)1.0f;
  scalar_type BG = -(scalar_type)3.0f * B2 * ec * FAC / (EASYPBE_BETA * G4);
  scalar_type BEC = B2 * FAC / (EASYPBE_BETA * G3);
  scalar_type Q8 = pow(Q5, 2) + EASYPBE_DELTA * Q4 * Q5 * T2;
  scalar_type Q9 = (scalar_type)1.0f + (scalar_type)2.0f * B * T2;

  //      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
  //      hRS = -RSTHRD*hB*BEC*ECRS
  //      FACT0 = 2.D0*DELT-6.D0*B
  //      FACT1 = Q5*Q9+Q4*Q9*Q9

  scalar_type hB =
      -EASYPBE_BETA * G3 * B * (T6 / Q8) * ((scalar_type)2.0f + B * T2);
  scalar_type hRS = -rsthrd * hB * BEC * ECRS;
  scalar_type FACT0 = (scalar_type)2.0f * EASYPBE_DELTA - (scalar_type)6.0f * B;
  scalar_type FACT1 = Q5 * Q9 + Q4 * Q9 * Q9;

  // hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
  // hRST = RSTHRD*T2*hBT*BEC*ECRS
  // hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
  // hT = 2.d0*BET*G3*Q9/Q8
  // hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
  scalar_type hBT = (scalar_type)2.0f * EASYPBE_BETA * G3 * T4 *
                    ((Q4 * Q5 * FACT0 - EASYPBE_DELTA * FACT1) / Q8) / Q8;
  scalar_type hRST = rsthrd * T2 * hBT * BEC * ECRS;
  scalar_type hz =
      (scalar_type)3.0f * GZ * h / G + hB * (BG * GZ + BEC * ECZET);
  scalar_type hT = (scalar_type)2.0f * EASYPBE_BETA * G3 * Q9 / Q8;
  scalar_type hZT =
      (scalar_type)3.0f * GZ * hT / G + hBT * (BG * GZ + BEC * ECZET);

  // FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
  // FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
  // hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
  // COMM = H+ HRS + HRST + T2*HT/6.D0 + 7.D0*T2*T*HTT/6.D0
  // PREF = HZ-GZ*T2*HT/G
  // FACT5 = GZ*(2.D0*HT+T*HTT)/G
  // COMM = COMM - PREF*ZET - UU*HTT - VV*HT - WW*(HZT-FACT5)
  scalar_type FACT2 = Q4 * Q5 + B * T2 * (Q4 * Q9 + Q5);
  scalar_type FACT3 = (scalar_type)2.0f * B * Q5 * Q9 + EASYPBE_DELTA * FACT2;
  scalar_type hTT = (scalar_type)4.0f * EASYPBE_BETA * G3 * t *
                    (2.0f * B / Q8 - ((Q9 / Q8) * FACT3) / Q8);
  COMM = h + hRS + hRST + T2 * hT / (scalar_type)6.0f +
         (scalar_type)7.0f * T2 * t * hTT / (scalar_type)6.0f;
  scalar_type PREF = hz - GZ * T2 * hT / G;
  scalar_type FACT5 = GZ * ((scalar_type)2.0f * hT + t * hTT) / G;
  COMM = COMM - PREF * zet - uu * hTT - vv * hT - ww * (hZT - FACT5);
//  if ( COMM != COMM ) { printf("NaN in COMM2 - hRS %E hRST %E T2 %E hT %E hTT %E \n",
// hRS, hRST, T2, hT, hTT); };
  dvc_a = COMM + PREF;
  dvc_b = COMM - PREF;
  // PBE POTENTIAL DONE !!!!!
  //--------------------------------------------------
}
}
