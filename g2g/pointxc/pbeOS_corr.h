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
namespace pointxc {
template<class scalar_type> __host__ __device__
void pbeOS_corr( scalar_type rho,   scalar_type rs,  scalar_type zet,
                 scalar_type t,     scalar_type uu,  scalar_type vv,
                 scalar_type ww,    scalar_type& ec, scalar_type& vc_a,
                 scalar_type& vc_b, scalar_type& h,  scalar_type& dvc_a,
                 scalar_type& dvc_b );

}
