#define GCORC_A_1 ((scalar_type)0.0310907f)
#define GCORC_A1_1 ((scalar_type)0.21370f)
#define GCORC_B1_1 ((scalar_type)7.5957f)
#define GCORC_B2_1 ((scalar_type)3.5876f)
#define GCORC_B3_1 ((scalar_type)1.6382f)
#define GCORC_B4_1 ((scalar_type)0.49294f)

#define GCORC_A_2 ((scalar_type)0.01554535f)
#define GCORC_A1_2 ((scalar_type)0.20548f)
#define GCORC_B1_2 ((scalar_type)14.1189f)
#define GCORC_B2_2 ((scalar_type)6.1977f)
#define GCORC_B3_2 ((scalar_type)3.3662f)
#define GCORC_B4_2 ((scalar_type)0.62517f)

#define GCORC_A_3 ((scalar_type)0.0168869f)
#define GCORC_A1_3 ((scalar_type)0.11125f)
#define GCORC_B1_3 ((scalar_type)10.357f)
#define GCORC_B2_3 ((scalar_type)3.6231f)
#define GCORC_B3_3 ((scalar_type)0.88026f)
#define GCORC_B4_3 ((scalar_type)0.49671f)


template<class scalar_type>
__device__ void gcor2(scalar_type A, scalar_type A1, scalar_type B1, scalar_type B2, scalar_type B3, scalar_type B4,scalar_type rtrs,scalar_type& gg, scalar_type& ggrs){
// CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,0.49294D0,rtrs,EU,EURS)
//----------------------------------------------------------------------
//      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
// slimmed down version of GCOR used in PW91 routines, to interpolate
// LSD correlation energy, as given by (10) of
// J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
// K. Burke, May 11, 1996.

  scalar_type Q0 = -(scalar_type)2.0f*A*((scalar_type)1.0f+A1*pow(rtrs,(scalar_type)2.0f));
  scalar_type Q1 = (scalar_type)2.0f*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)));
  scalar_type Q2 = log((scalar_type)1.0f+(scalar_type)1.0f/Q1);
  gg = Q0*Q2;

  scalar_type Q3 = A*(B1/rtrs+(scalar_type)2.0f*B2+rtrs*((scalar_type)3.0f*B3+(scalar_type)4.0f*B4*rtrs));

  ggrs = -(scalar_type)2.0f*A*A1*Q2 - Q0*Q3/(Q1*((scalar_type)1.0f+Q1));

  return;

}

#define EASYPBE_PI   ((scalar_type)3.141592653589793238462643383279f)
#define EASYPBE_PI32 ((scalar_type)29.608813203268075856503472999628f)
#define EASYPBE_AX   ((scalar_type)-0.738558766382022405884230032680836f)
#define EASYPBE_UM   ((scalar_type)0.2195149727645171f)
#define EASYPBE_UK   ((scalar_type)0.804f)
#define EASYPBE_UL   ((scalar_type)0.273028573090195f) // um / uk

template<class scalar_type>
__device__ void exchpbe(scalar_type rho, scalar_type s,scalar_type u, scalar_type v, scalar_type& ex, scalar_type& vx){
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

// construct LDA exchange energy density
  scalar_type rho13 = pow(rho, (scalar_type)(1.0f / 3.0f));
  scalar_type exunif = (scalar_type)EASYPBE_AX*rho13;
//----------------------------------------------------------------------
// construct PBE enhancement factor
  scalar_type s2 = pow(s, (scalar_type)2.0f);
  scalar_type s3 = pow(s, (scalar_type)3.0f);
  scalar_type p0 = 1.0f + EASYPBE_UL * s2;

  scalar_type fxpbe = 1.0f + EASYPBE_UK - EASYPBE_UK/p0;
//
  ex = exunif * fxpbe;
// EXCHANGE PBE ENERGY DONE !!!!!!!!!!

//===============================================================
//  NOW THE POTENTIAL:
//      if(lpot.eq.0)return
//  find first and second derivatives of Fx w.r.t s.
//  Fs=(1/s)*d FxPBE/ ds
//  Fss=d Fs/ds
  scalar_type p2 = p0 * p0;
  scalar_type Fs = 2.0f * EASYPBE_UM / p2;
  scalar_type F1 = -4.0f * EASYPBE_UL * s * Fs;
  scalar_type Fss = F1/p0;
//-----------------------------------------------------------
// calculate potential from [b](24)
  scalar_type vx2 = (4.0f / 3.0f) * fxpbe;
  scalar_type vx3 = v * Fs;
  scalar_type vx4 = (u - (4.0f / 3.0f) * s3) * Fss;
//
  vx = exunif * (vx2 - vx4 - vx3);
// EXCHANGE PBE POTENTIAL DONE !!!!!!!

  return;

}
#define EASYPBE_GAMMA    ((scalar_type)0.03109069086965489503494086371273f)
#define EASYPBE_GAMMAINV ((scalar_type)32.1639684429148f) // 1 / gamma
#define EASYPBE_BETA     ((scalar_type)0.06672455060314922f)
#define EASYPBE_DELTA    ((scalar_type)(EASYPBE_BETA/EASYPBE_GAMMA)) // beta/gamma
#define EASYPBE_ETA      ((scalar_type)0.000000000001f) //
template<class scalar_type>
__device__ void corpbe(scalar_type rho,scalar_type rs,scalar_type zet,scalar_type t,
    scalar_type uu,scalar_type vv,scalar_type ww,scalar_type& ec,scalar_type& vc_a,
    scalar_type& vc_b,scalar_type& h,scalar_type& dvc_a,scalar_type& dvc_b){

//printf("rho,rs,zet,t,uu,vv,ww, %14.8f %16.8f %10.4e %16.8f %20.8f %16.8f %16.10e\n",rho,rs,zet,t,uu,vv,ww);
//----------------------------------------------------------------------
//      SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,H,DVCUP,DVCDN,rho)
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

  gcor2(GCORC_A_1, GCORC_A1_1, GCORC_B1_1, GCORC_B2_1, GCORC_B3_1, GCORC_B4_1, rtrs, EU, EURS);
        gcor2(GCORC_A_2, GCORC_A1_2, GCORC_B1_2, GCORC_B2_2, GCORC_B3_2, GCORC_B4_2, rtrs, EP, EPRS);
        gcor2(GCORC_A_3, GCORC_A1_3, GCORC_B1_3, GCORC_B2_3, GCORC_B3_3, GCORC_B4_3, rtrs, ALFM, ALFRSM);

  scalar_type Z4 = pow((scalar_type)zet, 4);
  scalar_type ZET1 = cbrt(pow((scalar_type)(1.0f+zet), 4));
  scalar_type ZET2 = cbrt(pow((scalar_type)(1.0f-zet), 4));
  scalar_type GAM  = pow ((scalar_type)(2.0f), (scalar_type)(4.0f/3.0f))-(scalar_type)2.0f;

  scalar_type F =( ZET1 + ZET2 - (scalar_type)2.0f )/GAM;

  scalar_type FZZ = (scalar_type)8.0f/(9.0f*GAM);

  ec = EU*(1.0f-F*Z4) + EP*F*Z4 - ALFM*F*(1.0f-Z4)/FZZ;

//----------------------------------------------------------------------
// LSD potential from [c](A1)
// ECRS = dEc/drs [c](A2)
// ECZET=dEc/dzeta [c](A3)
// FZ = dF/dzeta [c](A4)
  scalar_type ECRS = EURS*(1.0f-F*Z4) + EPRS*F*Z4 - ALFRSM*F*(1.0f-Z4)/FZZ;
  scalar_type ZET3 = cbrt((scalar_type)(1.0f+zet));
  scalar_type ZET4 = cbrt((scalar_type)(1.0f-zet));
  scalar_type FZ = (scalar_type)(4.0f/3.0f) * ( ZET3 - ZET4 )/GAM;

  scalar_type ECZET = (scalar_type)4.0f*pow((scalar_type)zet,3)*F*(EP-EU+ALFM/FZZ) + FZ*(Z4*EP-Z4*EU - (1.0f-Z4)*ALFM/FZZ);

  scalar_type COMM = ec - rs*ECRS/(scalar_type)3.0f - zet*ECZET;
//----------------------------------------------------------------------
  vc_a = COMM + ECZET;
  vc_b = COMM - ECZET;
//----------------------------------------------------------------------
// PBE correlation energy
//----------------------------------------------------------------------
// G=phi(zeta), given after [a](3)
// DELT=bet/gamma
// B=A of [a](8)
//
  scalar_type ZET5 = cbrt(pow((scalar_type)(1.0f+zet), 2));
  scalar_type ZET6 = cbrt(pow((scalar_type)(1.0f-zet), 2));
  scalar_type G = (ZET5 + ZET6)/(scalar_type)2.0f;
  scalar_type G3 = pow(G, 3);
  scalar_type PON = -ec/(G3*EASYPBE_GAMMA);
  scalar_type B1 = exp(PON);
  scalar_type B = EASYPBE_DELTA/(B1 - (scalar_type)1.0f);

  scalar_type B2 = pow(B, 2);
  scalar_type T2 = pow(t, 2);
  scalar_type T4 = pow(T2, 2);
  scalar_type RS2 = pow(rs, 2);
  scalar_type Q4 = (scalar_type)1.0f + B * T2;
  scalar_type Q5 = (scalar_type)1.0f + B*T2 + B2*T4;

  h = G3*(scalar_type)(EASYPBE_BETA/EASYPBE_DELTA)*log((scalar_type)1.0f+EASYPBE_DELTA*Q4*T2/Q5);
// ENERGY DONE !!!!!!!
//=========================================================
//NOW THE POTENTIAL, using appendix E of [b].
//
  scalar_type G4=G3*G;
  scalar_type T6 = T4*T2;
  scalar_type rsthrd = rs/(scalar_type)3.0f;
  scalar_type ZET7 = pow((scalar_type)(1.0f+zet), 2);
  scalar_type ZET8 = pow((scalar_type)(1.0f-zet), 2);
  scalar_type GZ1 = pow((ZET7+EASYPBE_ETA), (scalar_type)(-1.0f/6.0f));
  scalar_type GZ2 = pow((ZET8+EASYPBE_ETA), (scalar_type)(-1.0f/6.0f));
  scalar_type GZ = (GZ1-GZ2)/(scalar_type)3.0f;

  scalar_type FAC = EASYPBE_DELTA/B+(scalar_type)1.0f;

  scalar_type BG = -(scalar_type)3.0f*B2*ec*FAC/(EASYPBE_BETA*G4);
  scalar_type BEC = B2*FAC/(EASYPBE_BETA*G3);

  scalar_type Q8 = pow(Q5, 2)+EASYPBE_DELTA*Q4*Q5*T2;
  scalar_type Q9 = (scalar_type)1.0f+(scalar_type)2.0f*B*T2;

  scalar_type hB = -EASYPBE_BETA*G3*B*T6*((scalar_type)2.0f+B*T2)/Q8;
  scalar_type hRS= -rsthrd*hB*BEC*ECRS;

  scalar_type FACT0 = (scalar_type)2.0f*EASYPBE_DELTA-(scalar_type)6.0f*B;
  scalar_type FACT1 = Q5*Q9+Q4*Q9*Q9;

  scalar_type hBT = (scalar_type)2.0f*EASYPBE_BETA*G3*T4*((Q4*Q5*FACT0-EASYPBE_DELTA*FACT1)/Q8)/Q8;
  scalar_type hRST = rsthrd*T2*hBT*BEC*ECRS;
  scalar_type hz = (scalar_type)3.0f*GZ*h/G + hB*(BG*GZ+BEC*ECZET);
  scalar_type hT = (scalar_type)2.0f*EASYPBE_BETA*G3*Q9/Q8;
  scalar_type hZT = (scalar_type)3.0f*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET);

  scalar_type FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5);
  scalar_type FACT3 = (scalar_type)2.0f*B*Q5*Q9+EASYPBE_DELTA*FACT2;

  scalar_type hTT = (scalar_type)4.0f*EASYPBE_BETA*G3*t*(2.0f*B/Q8-(Q9*FACT3/Q8)/Q8);

  COMM = h+ hRS + hRST + T2*hT/(scalar_type)6.0f + (scalar_type)7.0f*T2*t*hTT/(scalar_type)6.0f;

  scalar_type PREF = hz-GZ*T2*hT/G;

  scalar_type FACT5 = GZ*((scalar_type)2.0f*hT+t*hTT)/G;

  COMM = COMM - PREF*zet - uu*hTT - vv*hT - ww*(hZT-FACT5);
//--------------------------------------------------
  dvc_a = COMM + PREF;
  dvc_b = COMM - PREF;
//
// PBE POTENTIAL DONE !!!!!
//--------------------------------------------------
  return;
}

template<class scalar_type>
__device__ void easypbe(scalar_type dens_a, scalar_type dgrad_a, scalar_type delgrad_a,
    scalar_type rlap_a, scalar_type dens_b, scalar_type dgrad_b, scalar_type delgrad_b,
    scalar_type rlap_b, scalar_type dgrad, scalar_type delgrad, scalar_type& expbe,
    scalar_type& vxpbe_a, scalar_type& vxpbe_b,scalar_type& ecpbe,scalar_type& corr1,
    scalar_type& corr2, scalar_type& vcpbe_a, scalar_type& vcpbe_b){
//----------------------------------------------------------------------
// EASYPBE is a driver for the PBE subroutines, using simple inputs  (K. Burke, May 14, 1996.)
// inputs:
//       : dens_a =up density
//       : dgrad_a =|grad up|
//       : delgrad_a =(grad up).(grad |grad up|)
//       : rlap_a =grad^2 up=Laplacian of up
//       : dens_b,dgrad_b,delgrad_b,rlap_b =corresponding down quantities
//       : dgrad =|grad rho|
//       : delgrad =(grad rho).(grad |grad rho|)
//       :?? lcor=flag to do correlation(=0=>don't)
//       :?? lpot=flag to do potential(=0=>don't)
// outputs:
//       : expbe,vxpbe_a,vxpbe_b,ecpbe,vcpbe_a,vcpbe_b  =PBE quantities

//----------------------------------------------------------------------
// PBE exchange
//----------------------------------------------------------------------
// use  Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) (i.e., exact spin-scaling)
// Perdew et al. Phys Rev Lett 77,3865. 1996 (Eq. 11)
//
// do up exchange
// fk=local Fermi wavevector for 2*up
//   =(3 pi^2 (2up))^(1/3)
//
// s=dimensionless density gradient
//  =|grad rho|/ (2*fk*rho)  where (rho=2*up)
//
// u= ??????
//  =delgrad/(rho^2*(2*fk)**3)   where (rho=2*up)
//
// v=??????
//  =Laplacian/(rho*(2*fk)**2)  where  (rho=2*up)

  scalar_type expbe_a,expbe_b;

// Density Up
  scalar_type twodens = (scalar_type)2.0f * dens_a;
  scalar_type twodens2 = pow(twodens, 2);

  if (twodens > ((scalar_type)1e-18f)) {
    scalar_type rho13 = cbrt(twodens);
    scalar_type fk1 = cbrt((scalar_type)EASYPBE_PI32);
    scalar_type fk = fk1 * rho13;
    scalar_type twofk = (scalar_type)2.0f * fk;
    scalar_type twofk2 = pow(twofk, 2);
    scalar_type twofk3 = pow(twofk, 3);

    scalar_type s = ((scalar_type)2.0f * dgrad_a) / (twofk * twodens);
    scalar_type u = ((scalar_type)4.0f * delgrad_a ) / (twodens2 * twofk3);
    scalar_type v = ((scalar_type)2.0f * rlap_a) / (twodens * twofk2);
    exchpbe(twodens,s,u,v,expbe_a,vxpbe_a);
  }
  else {
    expbe_a=(scalar_type)0.0f;
    vxpbe_a=(scalar_type)0.0f;
  }

// Density Down
  twodens = (scalar_type)2.0f * dens_b;
  twodens2 = pow((scalar_type)twodens, 2);

  if (twodens > ((scalar_type)1e-18f)) {
    scalar_type rho13 = cbrt((scalar_type)twodens);
    scalar_type fk1 = cbrt((scalar_type)EASYPBE_PI32);
    scalar_type fk = fk1 * rho13;
    scalar_type twofk = (scalar_type)2.0f * fk;
    scalar_type twofk2 = pow((scalar_type)twofk, 2);
    scalar_type twofk3 = pow((scalar_type)twofk, 3);

    scalar_type s = ((scalar_type)2.0f * dgrad_b) / (twofk * twodens);
    scalar_type u = ((scalar_type)4.0f * delgrad_b ) / (twodens2 * twofk3);
    scalar_type v = ((scalar_type)2.0f * rlap_b) / (twodens * twofk2);

    exchpbe(twodens,s,u,v,expbe_b,vxpbe_b);
  }
  else {
    expbe_b=(scalar_type)0.0f;
    vxpbe_b=(scalar_type)0.0f;
  }

// construct total density and contribution to ex
  scalar_type rho = dens_a + dens_b;
  expbe=(expbe_a * dens_a + expbe_b * dens_b)/rho;

//----------------------------------------------------------------------
// Now do correlation
//----------------------------------------------------------------------
// zet=(up-dn)/rho
// g=phi(zeta)
// rs=local Seitz radius=alpha/fk
//   =(3/(4pi*rho))^(1/3)
// sk=Ks=Thomas-Fermi screening wavevector
//      =sqrt(4fk/pi)
// twoksg=2*Ks*phi
//
// t=correlation dimensionless gradient
//  =|grad rho|/(2*Ks*phi*rho)
// uu=delgrad/(rho^2*twoksg^3)
// rholap=Laplacian
// vv=Laplacian/(rho*twoksg^2)
// ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
//
// ec=lsd correlation energy
// vcup=lsd up correlation potential
// vcdn=lsd down correlation potential
// h=gradient correction to correlation energy
// dvcup=gradient correction to up correlation potential
// dvcdn=gradient correction to down correlation potential

  if(rho < ((scalar_type)1e-18f)){
    ecpbe = (scalar_type)0.0f; //ec+h
    vcpbe_a = (scalar_type)0.0f; // vc_a+dvc_a
    vcpbe_b = (scalar_type)0.0f; //vc_b+dvc_b
    return;
  }
  scalar_type zet=(dens_a-dens_b)/rho;
  scalar_type g1= (scalar_type)1.0f + zet;
  scalar_type g1_23=cbrt(pow(g1, 2));
  scalar_type g2= (scalar_type)1.0f - zet;
  scalar_type g2_23=cbrt(pow(g2, 2));
  scalar_type g=(g1_23+g2_23) / (scalar_type)2.0f;

  scalar_type rho13 = cbrt(rho);
  scalar_type fk1 = cbrt(EASYPBE_PI32);
  scalar_type fk = fk1 * rho13;

  scalar_type alpha = ((scalar_type)9.0f * (scalar_type)EASYPBE_PI) / (scalar_type)4.0f;
  scalar_type alpha13 = cbrt(alpha);
  scalar_type rs = alpha13/fk;

  scalar_type sk = sqrt(4.0f*fk/(scalar_type)EASYPBE_PI);
  scalar_type twoksg = (scalar_type)2.0f*sk*g;
  scalar_type t = dgrad/(twoksg*rho);

  scalar_type uu=delgrad/(pow(rho,2)*pow(twoksg,3));
  scalar_type rho_lap = rlap_a + rlap_b;
  scalar_type vv = rho_lap/(rho*pow(twoksg,2));

// C----reconstruction of ww-----------------------
  scalar_type ww =
    (pow(dgrad_a,2)-pow(dgrad_b, 2)-zet*pow(dgrad,2)) /
    (pow(rho,2)*pow(twoksg,2));

  scalar_type ec,h;
  scalar_type vc_a,vc_b,dvc_a,dvc_b;

  corpbe<scalar_type>(rho,rs,zet,t,uu,vv,ww,ec,vc_a,vc_b,h,dvc_a,dvc_b);

  ecpbe = ec+h;
  corr1 = ec;
  corr2 = h;
  vcpbe_a = vc_a+dvc_a;
  vcpbe_b = vc_b+dvc_b;

  return;
}


template<class scalar_type>
__device__ void gpu_pot_open(scalar_type dens_a, scalar_type dens_b, const vec_type<scalar_type,4>& grad_a,
    const vec_type<scalar_type,4>& grad_b, const vec_type<scalar_type,4>& hess1_a,const vec_type<scalar_type,4>& hess1_b,
    const vec_type<scalar_type,4>& hess2_a,const vec_type<scalar_type,4>& hess2_b, scalar_type& exc_corr,
    scalar_type& exc,scalar_type& corr,scalar_type& corr1,scalar_type& corr2, scalar_type& v_a, scalar_type& v_b)
{
// Default values
// v_a, v_b are y2a, y2b

  scalar_type expbe,vxpbe_a,vxpbe_b,ecpbe,vcpbe_a,vcpbe_b;
  expbe = vxpbe_a = vxpbe_b = ecpbe = vcpbe_a = vcpbe_b = exc_corr = exc = corr = corr1 = corr2 = v_a = v_b = 0.0f;
  scalar_type dgrad_a,delgrad_a,rlap_a,dgrad_b,delgrad_b,rlap_b;
//===============================================================
// THRESHOLD VALUE for density ...
  if ( (dens_a + dens_b) > ((scalar_type)1e-13f) ) {
//===============================================================
// PBE (gpu_Iexch == 9)
    if (gpu_Iexch == 9) {
      vec_type<scalar_type,4> grad;
      grad.x = grad_a.x + grad_b.x;
      grad.y = grad_a.y + grad_b.y;
      grad.z = grad_a.z + grad_b.z;
      vec_type<scalar_type,4> hess1;
      hess1.x = hess1_a.x + hess1_b.x;
      hess1.y = hess1_a.y + hess1_b.y;
      hess1.z = hess1_a.z + hess1_b.z;
      vec_type<scalar_type,4> hess2;
      hess2.x = hess2_a.x + hess2_b.x;
      hess2.y = hess2_a.y + hess2_b.y;
      hess2.z = hess2_a.z + hess2_b.z;

//   Up density
      if (dens_a == ((scalar_type)0.0f)) {
        dgrad_a = 0.0f;
        rlap_a = 0.0f;
        delgrad_a = 0.0f;
      }
      else{
        scalar_type grad2_a = pow(grad_a.x, 2) + pow(grad_a.y, 2) + pow(grad_a.z, 2);
        if (grad2_a == (scalar_type)0.0f) grad2_a = FLT_MIN;
        dgrad_a = sqrt(grad2_a);
//   Laplacian Up
        rlap_a = hess1_a.x + hess1_a.y + hess1_a.z;
        delgrad_a = (pow(grad_a.x, 2)*hess1_a.x + (scalar_type)2.0f*grad_a.x*grad_a.y*hess2_a.x +
            (scalar_type)2.0f*grad_a.y*grad_a.z*hess2_a.z + (scalar_type)2.0f*grad_a.x*grad_a.z*hess2_a.y +
            pow(grad_a.y, 2)*hess1_a.y + pow(grad_a.z, 2)*hess1_a.z) / dgrad_a;
      }

//   Down density
      if (dens_b == ((scalar_type)0.0f)) {
        dgrad_b = 0.0f;
        rlap_b = 0.0f;
        delgrad_b = 0.0f;
      }
      else {
        scalar_type grad2_b = pow(grad_b.x, 2)+ pow(grad_b.y, 2)+ pow(grad_b.z, 2);
        if (grad2_b == (scalar_type)0.0f) grad2_b = FLT_MIN;
        dgrad_b = sqrt(grad2_b);
//   Laplacian Down
        rlap_b = hess1_b.x + hess1_b.y + hess1_b.z;
        delgrad_b = (pow(grad_b.x, 2)*hess1_b.x + (scalar_type)2.0f*grad_b.x*grad_b.y*hess2_b.x +
            (scalar_type)2.0f*grad_b.y*grad_b.z*hess2_b.z + (scalar_type)2.0f*grad_b.x*grad_b.z*hess2_b.y +
            pow(grad_b.y, 2)*hess1_b.y + pow(grad_b.z, 2)*hess1_b.z) / dgrad_b;
      }
//   Up + Down densities
      scalar_type grad2 =  pow(grad.x, 2) + pow(grad.y, 2) + pow(grad.z, 2);
      if (grad2 == (scalar_type)0.0f) grad2 = FLT_MIN;
      scalar_type dgrad = sqrt(grad2);
      scalar_type delgrad = (pow(grad.x, 2)*hess1.x + pow(grad.y, 2)*hess1.y +
          pow(grad.z, 2)*hess1.z + (scalar_type)2.0f*grad.x*grad.y*hess2.x +
          (scalar_type)2.0f*grad.y*grad.z*hess2.z + (scalar_type)2.0f*grad.x*grad.z*hess2.y)/dgrad;
      easypbe<scalar_type>(dens_a, dgrad_a, delgrad_a, rlap_a, dens_b, dgrad_b, delgrad_b,
                           rlap_b, dgrad, delgrad, expbe, vxpbe_a, vxpbe_b, ecpbe, corr1, corr2,
                           vcpbe_a, vcpbe_b);
    }
    else {
//       NO OTHER XC FUNCIONAL
    }
  }
  exc=expbe;
  corr=ecpbe;
  exc_corr = expbe + ecpbe;

  v_a = vxpbe_a + vcpbe_a;
  v_b = vxpbe_b + vcpbe_b;

  return;
}
