#define GCORC_A_1 ((scalar_type)0.0310907f)
#define GCORC_A1_1 ((scalar_type)0.21370f)
#define GCORC_B1_1 ((scalar_type)7.5957f)
#define GCORC_B2_1 ((scalar_type)3.5876f)
#define GCORC_B3_1 ((scalar_type)1.6382f)
#define GCORC_B4_1 ((scalar_type)0.49294f)

#define GCORC_A_2 ((scalar_type)0.0310907f)
#define GCORC_A1_2 ((scalar_type)0.21370f)
#define GCORC_B1_2 ((scalar_type)7.5957f)
#define GCORC_B2_2 ((scalar_type)3.5876f)
#define GCORC_B3_2 ((scalar_type)1.6382f)
#define GCORC_B4_2 ((scalar_type)0.49294f)

#define GCORC_A_3 ((scalar_type)0.0310907f)
#define GCORC_A1_3 ((scalar_type)0.21370f)
#define GCORC_B1_3 ((scalar_type)7.5957f)
#define GCORC_B2_3 ((scalar_type)3.5876f)
#define GCORC_B3_3 ((scalar_type)1.6382f)
#define GCORC_B4_3 ((scalar_type)0.49294f)


template<class scalar_type>
__device__ void gcor2(scalar_type A, scalar_type A1, scalar_type B1, scalar_type B2, scalar_type B3, scalar_type B4,scalar_type rtrs,scalar_type& gg, scalar_type& ggrs){
// CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,0.49294D0,rtrs,EU,EURS)
//----------------------------------------------------------------------
//      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
// slimmed down version of GCOR used in PW91 routines, to interpolate
// LSD correlation energy, as given by (10) of
// J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
// K. Burke, May 11, 1996.
    
	scalar_type Q0 = -2.0f*A*(1.0f+A1+rtrs*rtrs);
	scalar_type Q1 = 2.0f*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)));
	scalar_type Q2 = log(1.0f+1.0f/Q1);
	gg = Q0*Q2;

	scalar_type Q3 = A*(B1/rtrs+2.0f*B2+rtrs*(3.0f*B3+4.0f*B4*rtrs));
 
	ggrs = -2.0f*A*A1*Q2 -Q0*Q3/(Q1*(1.0f+Q1));

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
//      if(lgga.eq.0)then
//        ex=exunif
//        vx=ex*thrd4
//        return
//      endif
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
#define EASYPBE_DELTA    ((scalar_type)2.14612633996736f) // beta/gamma
#define EASYPBE_ETA      ((scalar_type)1.0f-12) // 
template<class scalar_type>
__device__ void corpbe(scalar_type rho,scalar_type rs,scalar_type zet,scalar_type t,scalar_type uu,scalar_type vv,scalar_type ww,scalar_type& ec,scalar_type& vc_a,scalar_type& vc_b,scalar_type& h,scalar_type& dvc_a,scalar_type& dvc_b){
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
	
	//gcor2(0.0310907f ,0.21370f, 7.5957f,3.5876f,1.6382f ,0.49294f,rtrs,EU,EURS);
	//gcor2(0.01554535f,0.20548f,14.1189f,6.1977f,3.3662f ,0.62517f,rtrs,EP,EPRS);
	//gcor2(0.0168869f ,0.11125f,10.357f ,3.6231f,0.88026f,0.49671f,rtrs,ALFM,ALFRSM);

	gcor2(GCORC_A_1, GCORC_A1_1, GCORC_B1_1, GCORC_B2_1, GCORC_B3_1, GCORC_B4_1, rtrs, EU, EURS);
        gcor2(GCORC_A_2, GCORC_A1_2, GCORC_B1_2, GCORC_B2_2, GCORC_B3_2, GCORC_B4_2, rtrs, EP, EPRS);
        gcor2(GCORC_A_3, GCORC_A1_3, GCORC_B1_3, GCORC_B2_3, GCORC_B3_3, GCORC_B4_3, rtrs, ALFM, ALFRSM);

//	scalar_type ALFC = -ALFM;
	scalar_type Z4 = pow(zet, (scalar_type)4.0);

	scalar_type ZET1 = pow((scalar_type)(1.0f+zet), (scalar_type)(4.0f/3.0f));
	scalar_type ZET2 = pow((scalar_type)(1.0f-zet), (scalar_type)(4.0f/3.0f));
	scalar_type GAM = pow ((scalar_type)(2.0f), (scalar_type)((4.0f/3.0f)-2.0f));

	scalar_type F =( ZET1 + ZET2 - 2.0f )/GAM;
//	F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM

	scalar_type FZZ = 8.0f/(9.0f*GAM);

	ec = EU*(1.0f-F*Z4) + EP*F*Z4 - ALFM*F*(1.0f-Z4)/FZZ;
//      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ

//----------------------------------------------------------------------
// LSD potential from [c](A1)
// ECRS = dEc/drs [c](A2)
// ECZET=dEc/dzeta [c](A3)
// FZ = dF/dzeta [c](A4)
	scalar_type ECRS = EURS*(1.0f-F*Z4) + EPRS*F*Z4 - ALFRSM*F*(1.0f-Z4)/FZZ;
//      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ

	scalar_type ZET3 = pow((scalar_type)(1.0f+zet), (scalar_type)(1.0f/3.0f));
        scalar_type ZET4 = pow((scalar_type)(1.0f-zet), (scalar_type)(1.0f/3.0f));
	scalar_type FZ = (4.0f/3.0f) * ( ZET3 - ZET4 )/GAM;
//      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM

	scalar_type ECZET = 4.0f*(zet*zet*zet)*F*(EP-EU+ALFM/FZZ) + FZ*(Z4*EP-Z4*EU - (1.0f-Z4)*ALFM/FZZ);
//      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU-(1.D0-Z4)*ALFM/FZZ)

	scalar_type COMM = ec - rs*ECRS/3.0f - zet*ECZET;
//      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
//----------------------------------------------------------------------
	vc_a = COMM + ECZET;
	vc_b = COMM - ECZET;
//----------------------------------------------------------------------
// PBE correlation energy
//
//     if(lgga.eq.0)return
//----------------------------------------------------------------------
// G=phi(zeta), given after [a](3)
// DELT=bet/gamma
// B=A of [a](8)
//
//      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
	scalar_type ZET5 = pow((scalar_type)(1.0f+zet), (scalar_type)(2.0f/3.0f));
	scalar_type ZET6 = pow((scalar_type)(1.0f-zet), (scalar_type)(2.0f/3.0f));
	scalar_type G = (ZET5 + ZET6)/2.0f;
//      G3 = G**3
	scalar_type G3 = pow(G, (scalar_type)3.0f);
//      PON=-EC/(G3*gama)
	scalar_type PON = -ec/(G3*EASYPBE_GAMMA);
//      B = DELT/(DEXP(PON)-1.D0)
	scalar_type B1 = exp(PON);
	scalar_type B = EASYPBE_DELTA/(B1 - 1.0f);
//      B2 = B*B
	scalar_type B2 = pow(B, (scalar_type)2.0f);
//      T2 = T*T
	scalar_type T2 = pow(t, (scalar_type)2.0f);
//      T4 = T2*T2
	scalar_type T4 = pow(T2, (scalar_type)2.0f);
//      RS2 = RS*RS
	scalar_type RS2 = pow(rs, (scalar_type)2.0f);
//      RS3 = RS2*RS
//	scalar_type RS3 = RS2*rs;
//      Q4 = 1.D0+B*T2
	scalar_type Q4 = 1.0f + B * T2;
//      Q5 = 1.D0+B*T2+B2*T4
	scalar_type Q5 = 1.0f + B*T2 + B2*T4;

//      H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
	h = G3*(EASYPBE_BETA/EASYPBE_DELTA)*log(1.0f+EASYPBE_DELTA*Q4*T2/Q5);

// ENERGY DONE !!!!!!! 
//=========================================================
//NOW THE POTENTIAL, using appendix E of [b].
//
//	if(lpot.eq.0)return
//	G4 = G3*G
	scalar_type G4=G3*G;
//	T6 = T4*T2
	scalar_type T6 = T4*T2;
//	RSTHRD = RS/3.D0
	scalar_type rsthrd = rs/3.0f;
//	GZ=(((1.d0+zet)**2+eta)**sixthm-((1.d0-zet)**2+eta)**sixthm)/3.d0
	scalar_type ZET7 = pow((scalar_type)(1.0f+zet), (scalar_type)2.0f);
	scalar_type ZET8 = pow((scalar_type)(1.0f-zet), (scalar_type)2.0f);
	scalar_type GZ1 = pow((ZET7+EASYPBE_ETA), (scalar_type)(-1.0f/6.0f));
	scalar_type GZ2 = pow((ZET8+EASYPBE_ETA), (scalar_type)(-1.0f/6.0f));
	scalar_type GZ = (GZ1-GZ2)/3.0f;

//	FAC = DELT/B+1.D0
	scalar_type FAC = EASYPBE_DELTA/(B+1.0f);

//	BG = -3.D0*B2*EC*FAC/(BET*G4)
	scalar_type BG = -3.0f*B2*ec*FAC/(EASYPBE_BETA*G4);
//	BEC = B2*FAC/(BET*G3)
	scalar_type BEC = B2*FAC/(EASYPBE_BETA*G3);

//	Q8 = Q5*Q5+DELT*Q4*Q5*T2
	scalar_type Q8 = Q5*Q5+EASYPBE_DELTA*Q4+Q5*T2;

//	Q9 = 1.D0+2.D0*B*T2
	scalar_type Q9 = 1.0f+2.0f*B*T2;

//      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
	scalar_type hB = -EASYPBE_BETA*G3*B*T6*(2.0f+B*T2)/Q8;

//      hRS = -RSTHRD*hB*BEC*ECRS
	scalar_type hRS= -rsthrd*hB*BEC*ECRS;

//      FACT0 = 2.D0*DELT-6.D0*B
	scalar_type FACT0 = 2.0f*EASYPBE_DELTA-6.0f*B;

//      FACT1 = Q5*Q9+Q4*Q9*Q9
	scalar_type FACT1 = Q5*Q9+Q4*Q9*Q9;

//      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
	scalar_type hBT = 2.0f*EASYPBE_BETA*G3*T4*((Q4*Q5*FACT0-EASYPBE_DELTA*FACT1)/Q8)/Q8;

//      hRST = RSTHRD*T2*hBT*BEC*ECRS
	scalar_type hRST = rsthrd*T2*hBT*BEC*ECRS;

//      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
	scalar_type hz = 3.0f*GZ*h/G + hB*(BG*GZ+BEC*ECZET);

//      hT = 2.d0*BET*G3*Q9/Q8
	scalar_type hT = 2.0f*EASYPBE_BETA*G3*Q9/Q8;

//      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
	scalar_type hZT = 3.0f*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET);

//      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
	scalar_type FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5);

//      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
	scalar_type FACT3 = 2.0f*B*Q5*Q9+EASYPBE_DELTA*FACT2;

//      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
	scalar_type hTT = 4.0f*EASYPBE_BETA*G3*t*(2.0f*B/Q8-(Q9*FACT3/Q8)/Q8);

//      COMM = H+ HRS + HRST + T2*HT/6.D0 + 7.D0*T2*T*HTT/6.D0
	COMM = h+ hRS + hRST + T2*hT/6.0f + 7.0f*T2*t*hTT/6.0f;

//      PREF = HZ-GZ*T2*HT/G
	scalar_type PREF = hz-GZ*T2*hT/G;

//      FACT5 = GZ*(2.D0*HT+T*HTT)/G
	scalar_type FACT5 = GZ*(2.0f*hT+t*hTT)/G;

//      COMM = COMM - PREF*ZET - UU*HTT - VV*HT - WW*(HZT-FACT5)
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
__device__ void easypbe(scalar_type dens_a, scalar_type dgrad_a, scalar_type delgrad_a, scalar_type rlap_a, scalar_type dens_b, scalar_type dgrad_b, scalar_type delgrad_b, scalar_type rlap_b, scalar_type dgrad, scalar_type delgrad, scalar_type& expbe, scalar_type& vxpbe_a, scalar_type& vxpbe_b,scalar_type& ecpbe, scalar_type& vcpbe_a, scalar_type& vcpbe_b){
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
	scalar_type twodens = 2.0f * dens_a;
	scalar_type twodens2 = pow(twodens, (scalar_type)2.0f);	

	if (twodens > ((scalar_type)1e-18f)) {
//   		fk=(pi32*rho2)**thrd
		scalar_type rho13 = pow(twodens, (scalar_type)(1.0f / 3.0f));
  		scalar_type fk1 = pow((scalar_type)EASYPBE_PI32, (scalar_type)(1.0f / 3.0f));
  		scalar_type fk = fk1 * rho13;
		scalar_type twofk = 2.0f * fk;
		scalar_type twofk2 = pow(twofk, (scalar_type)2.0f);
		scalar_type twofk3 = pow(twofk, (scalar_type)3.0f);


//		s=2.d0*agrup/(2.d0*fk*rho2)
		scalar_type s = (2.0f * dgrad_a) / (twofk * twodens);

//		u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
		scalar_type u = (4.0f * delgrad_a ) / (twodens2 * twofk3);

//		v=2.d0*uplap/(rho2*(2.d0*fk)**2)
		scalar_type v = (2.0f * rlap_a) / (twodens * twofk2);
	
		exchpbe(twodens,s,u,v,expbe_a,vxpbe_a);
	}	
	else{
        	expbe_a=0.0f;
        	vxpbe_a=0.0f;
	}

// Density Down
	twodens = 2.0f * dens_b;
        twodens2 = pow(twodens, (scalar_type)2.0f);

        if (twodens > ((scalar_type)1e-18f)) {
//              fk=(pi32*rho2)**thrd
                scalar_type rho13 = pow(twodens, (scalar_type)(1.0f / 3.0f));
                scalar_type fk1 = pow((scalar_type)EASYPBE_PI32, (scalar_type)(1.0f / 3.0f));
                scalar_type fk = fk1 * rho13;
                scalar_type twofk = 2.0f * fk;
                scalar_type twofk2 = pow(twofk, (scalar_type)2.0f);
                scalar_type twofk3 = pow(twofk, (scalar_type)3.0f);

//              s=2.d0*agrup/(2.d0*fk*rho2)
                scalar_type s = (2.0f * dgrad_b) / (twofk * twodens);

//              u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
                scalar_type u = (4.0f * delgrad_b ) / (twodens2 * twofk3);

//              v=2.d0*uplap/(rho2*(2.d0*fk)**2)
                scalar_type v = (2.0f * rlap_b) / (twodens * twofk2);

	        exchpbe(twodens,s,u,v,expbe_b,vxpbe_b);
        }
        else{
                expbe_b=0.0f;
                vxpbe_b=0.0f;
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

//	if(rho < ((scalar_type)1e-18f)){
//		ecpbe = 0 //ec+h
//        	vcpbe_a = 0 // vc_a+dvc_a
//        	vcpbe_b = 0 //vc_b+dvc_b
//	return;
//	}

	scalar_type zet=(dens_a-dens_b)/rho;
	scalar_type g1= 1.0f + zet;
	scalar_type g1_23=pow(g1, (scalar_type)(2.0f / 3.0f));
	scalar_type g2= 1.0f - zet;
	scalar_type g2_23=pow(g2, (scalar_type)(2.0f / 3.0f));
	scalar_type g=(g1_23+g2_23) / 2.0f;
	
//	fk=(pi32*rho)**thrd
	scalar_type rho13 = pow(rho, (scalar_type)(1.0f / 3.0f));
        scalar_type fk1 = pow(EASYPBE_PI32, (scalar_type)(1.0f / 3.0f));
        scalar_type fk = fk1 * rho13;
//        twofk = 2.0f * fk;
//        twofk2 = pow(twofk, (scalar_type)2.0f);
//        twofk3 = pow(twofk, (scalar_type)3.0f);	

// alpha=(9pi/4)**(1/3)
	scalar_type alpha = (9.0f * (scalar_type)EASYPBE_PI) / 4.0f;
	scalar_type alpha13 = pow(alpha, (scalar_type)(1.0f / 3.0f));
	scalar_type rs = alpha13/fk;
	scalar_type sk = sqrt(4.0f*fk/(scalar_type)EASYPBE_PI);
//
	scalar_type twoksg = 2.0f*sk*g;
	scalar_type t = dgrad/(twoksg*rho);

	scalar_type uu=delgrad/(rho*rho*twoksg*twoksg*twoksg);

	scalar_type rho_lap = rlap_a + rlap_b;
	scalar_type vv = rho_lap/(rho*twoksg*twoksg);

// C----reconstruction of ww-----------------------

	scalar_type ww = ((dgrad_a*dgrad_a)-(dgrad_b*dgrad_b)-(zet*dgrad*dgrad)) / (rho*rho*twoksg*twoksg);

//      call CORPBE(RS,ZET,T,UU,VV,WW,1,lpot,ec,vcup,vcdn,H,DVCUP,DVCDN,rho)
	
	scalar_type ec,h;
	scalar_type vc_a,vc_b,dvc_a,dvc_b;

	corpbe(rho,rs,zet,t,uu,vv,ww,ec,vc_a,vc_b,h,dvc_a,dvc_b);

	ecpbe = ec+h;
	vcpbe_a = vc_a+dvc_a;
	vcpbe_b = vc_b+dvc_b;

	return;
}


template<class scalar_type>
__device__ void gpu_potop(scalar_type dens_a, scalar_type dens_b, const vec_type<scalar_type,4>& grad_a, const vec_type<scalar_type,4>& grad_b, const vec_type<scalar_type,4>& hess1_a,const vec_type<scalar_type,4>& hess1_b, const vec_type<scalar_type,4>& hess2_a,const vec_type<scalar_type,4>& hess2_b, scalar_type& exc_corr, scalar_type& v_a, scalar_type& v_b)
{
// Default values
// v_a, v_b son y2a, y2b
//	ex = ec = v_a = v_b = 0.0f;
	exc_corr = v_a = v_b = 0.0f;
	scalar_type dgrad_a,delgrad_a,rlap_a,dgrad_b,delgrad_b,rlap_b;

// PBE (gpu_Iexch == 9)
	if (gpu_Iexch == 9){

		scalar_type dens = dens_a + dens_b;
//===============================================================
// VALOR UMBRAL DEL DENSIDAD  para continuar ...
		if (dens < ((scalar_type)1e-18f)) return;
//===============================================================
//		scalar_type dens2 = dens*dens;
	
		vec_type<scalar_type,4> grad;
		grad.x = grad_a.x + grad_b.x;
		grad.y = grad_a.y + grad_b.y;
		grad.z = grad_a.z + grad_b.z;

// // hess1: xx, yy, zz
// // hess2: xy, xz, yz
		vec_type<scalar_type,4> hess1;
		hess1.x = hess1_a.x + hess1_b.x;
		hess1.y = hess1_a.y + hess1_b.y;
		hess1.z = hess1_a.z + hess1_b.z;
		vec_type<scalar_type,4> hess2;
        	hess2.x = hess2_a.x + hess2_b.x;
       		hess2.y = hess2_a.y + hess2_b.y;
        	hess2.z = hess2_a.z + hess2_b.z;

//	scalar_type y = pow(dens, (1.0f/3.0f));  // rho^(1/3)

// Up density
		if (dens_a == ((scalar_type)0.0f)){ 
			dgrad_a = 0.0f;
			rlap_a = 0.0f;
			delgrad_a = 0.0f;
		}
		else{
			scalar_type grad2_a = grad_a.x * grad_a.x + grad_a.y * grad_a.y + grad_a.z * grad_a.z;
			if (grad2_a == (scalar_type)0.0f) grad2_a = FLT_MIN;
			dgrad_a = sqrt(grad2_a);
// Laplacian Up
			rlap_a = hess1_a.x + hess1_a.y + hess1_a.z;
// misterious thing !!! Up
			delgrad_a = ((grad_a.x * grad_a.x) * hess1_a.x + 2.0f * grad_a.x * grad_a.y * hess2_a.x + 2.0f * grad_a.y * grad_a.z * hess2_a.z + 2.0f * grad_a.x * grad_a.z * hess2_a.y + (grad_a.y * grad_a.y) * hess1_a.y + (grad_a.z * grad_a.z) * hess1_a.z) / dgrad_a;
		}

// Down density
        	if (dens_b == ((scalar_type)0.0f)){
        		dgrad_b = 0.0f;
        		rlap_b = 0.0f;
        		delgrad_b = 0.0f;
        	}
        	else{
        		scalar_type grad2_b = grad_b.x * grad_b.x + grad_b.y * grad_b.y + grad_b.z * grad_b.z;
        		if (grad2_b == (scalar_type)0.0f) grad2_b = FLT_MIN;
        		dgrad_b = sqrt(grad2_b);
// Laplacian Down
        		rlap_b = hess1_b.x + hess1_b.y + hess1_b.z;
// misterious thing !!! Down
        		delgrad_b = ((grad_b.x * grad_b.x) * hess1_b.x + 2.0f * grad_b.x * grad_b.y * hess2_b.x + 2.0f * grad_b.y * grad_b.z * hess2_b.z + 2.0f * grad_b.x * grad_b.z * hess2_b.y + (grad_b.y * grad_b.y) * hess1_b.y + (grad_b.z * grad_b.z) * hess1_b.z) / dgrad_b;
        	}

// Up + Down densities
//
		scalar_type grad2 = grad.x * grad.x + grad.y * grad.y + grad.z * grad.z;
		if (grad2 == (scalar_type)0.0f) grad2 = FLT_MIN;
		scalar_type dgrad = sqrt(grad2);
//
//		scalar_type rlap = hess1.x + hess1.y + hess1.z;
		scalar_type delgrad = ((grad.x * grad.x) * hess1.x + 2.0f * grad.x * grad.y * hess2.x + 2.0f * grad.y * grad.z * hess2.z + 2.0f * grad.x * grad.z * hess2.y + (grad.y * grad.y) * hess1.y + (grad.z * grad.z) * hess1.z) / dgrad;

// 

		scalar_type expbe,vxpbe_a,vxpbe_b,ecpbe,vcpbe_a,vcpbe_b;
   
		easypbe(dens_a,dgrad_a,delgrad_a,rlap_a,dens_b,dgrad_b,delgrad_b,rlap_b,dgrad,delgrad,expbe,vxpbe_a,vxpbe_b,ecpbe,vcpbe_a,vcpbe_b);

      		//ex = expbe;
		//ec = ecpbe;
        	exc_corr = expbe + ecpbe;
     		v_a = vxpbe_a + vcpbe_a;
		v_b = vxpbe_b + vcpbe_b;

		return;
	}
	else{
// NO HAY IMPLEMENTADO OTRO FUNCIONAL DE XC	
	}
}
