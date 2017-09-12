/* -*- mode: c -*- */
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <limits>
#include "../common.h"
#include "../init.h"
#include "../cuda/cuda_extra.h"
#include "../matrix.h"
using namespace std;

namespace G2G {

#define POT_ALPHA 		((scalar_type)-0.738558766382022447) // -(3/PI)^(1/3)
#define POT_GL 				((scalar_type)0.620350490899400087)

#define POT_VOSKO_A1 	((scalar_type)0.03109205)
#define POT_VOSKO_B1 	((scalar_type)3.72744)
#define POT_VOSKO_C1 	((scalar_type)12.9352)
#define POT_VOSKO_X0 	((scalar_type)-0.10498)

#define POT_VOSKO_Q 	((scalar_type)6.15199066246304849)
#define POT_VOSKO_A16 ((scalar_type)0.005182008333)
#define POT_VOSKO_A2 	((scalar_type)0.015546025)
#define POT_VOSKO_B2 	((scalar_type)7.06042)
#define POT_VOSKO_C2 	((scalar_type)18.0578)
#define POT_VOSKO_X02 ((scalar_type)-0.32500)
#define POT_VOSKO_Q2 	((scalar_type)4.7309269)
#define POT_VOSKO_A26 ((scalar_type)0.0025910042)

#define POT_XX0 ((scalar_type)12.5549141492) // POT_VOSKO_X0 * POT_VOSKO_X0 + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1
#define POT_T6 ((scalar_type)-0.00836166609762834) // POT_VOSKO_X0 / POT_XX0
#define POT_T4 ((scalar_type)-0.0311676086789438) // POT_VOSKO_B1 * POT_VOSKO_X0 / POT_XX0
#define POT_VOSKO_2C1 ((scalar_type)25.8704) // 2 * POT_VOSKO_C1

#define POT_VOSKO_2B1Q ((scalar_type)1.21178337371132) // 2 * POT_VOSKO_B1 / POT_VOSKO_Q
#define POT_VOSKO_B2X0Q ((scalar_type)1.14352579286644) // 2 * (POT_VOSKO_B1 + 2 * POT_VOSKO_X0) / POT_VOSKO_Q
#define POT_VOSKO_4B1 ((scalar_type)14.90976) // 4.0 * POT_VOSKO_B1
#define POT_VOSKO_QSQ ((scalar_type)37.8469891110325) // POT_VOSKO_Q * POT_VOSKO_Q
#define POT_VOSKO_B1X0 ((scalar_type)1.0329232240928) // (1.0f - t6 * (POT_VOSKO_B1 - 2.0f * POT_VOSKO_X0))

template<class scalar_type, int iexch>
void cpu_pot_tmpl(scalar_type dens, scalar_type& ex, scalar_type& ec, scalar_type& y2a)
{
	// data X alpha

	if (dens == 0) {
		ex = 0; ec = 0;
		y2a = 0;
		return;
	}

	scalar_type y = pow(dens, (scalar_type)0.333333333333333333);  // rho^(1/3)
	scalar_type v0 = (scalar_type)-0.984745021842697 * y; // -4/3 * (3/PI)^(1/3) * rho^(1/3)

	ex = POT_ALPHA * y; // -(3/PI)^(1/3) * rho^(1/3)

	switch(iexch) {
		case 1:
		{
			ec = 0;
			y2a = v0;
		}
		break;
		case 2:
		{
			scalar_type rs = POT_GL / y;
			scalar_type x1 = rs / (scalar_type)11.4;
			scalar_type vc;

			if (x1 > 1.0) {
				ec = (scalar_type)-0.0333 * ((scalar_type)0.5 * x1 - (scalar_type)0.33333333333333);
				vc = (scalar_type)0.0111 * x1 * (scalar_type)0.5;
			}
			else {
				scalar_type t1 = ((scalar_type)1.0 + x1 * x1 * x1);
				scalar_type t2 = logf((scalar_type)1.0 + (scalar_type)1.0 / x1);
				scalar_type t3 = x1 * x1;
                                ec = (scalar_type)-0.0333 * (t1 * t2 - t3 + (scalar_type)0.5 * x1 - (scalar_type)0.33333333333333);
                                vc = (scalar_type)0.0111 * x1 * ((scalar_type)3.0 * t3 * t2 - t1 / (x1 * (x1 + (scalar_type)1.0)) - (scalar_type)2.0 * x1 + (scalar_type)0.5);
			}
			y2a = v0 + ec + vc;
		}
		break;
		case 3:
		{
			scalar_type rs = POT_GL / y;
			scalar_type x1 = sqrt(rs);
			scalar_type Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
			scalar_type t1 = (scalar_type)2.0 * x1 + POT_VOSKO_B1;
			scalar_type t2 = log(Xx);
			scalar_type t3 = atan(POT_VOSKO_Q/t1);
                        scalar_type t5 = (POT_VOSKO_B1 * x1 + POT_VOSKO_2C1) / x1;

                        ec = POT_VOSKO_A1 * ((scalar_type)2.0 * logf(x1) - t2 + POT_VOSKO_2B1Q * t3 - POT_T4 * ((scalar_type)2.0 * log(x1 - POT_VOSKO_X0) - t2 + POT_VOSKO_B2X0Q * t3));

			scalar_type vc;
                        vc = ec - POT_VOSKO_A16 * x1 * (t5 / Xx - POT_VOSKO_4B1 / (t1 * t1 + POT_VOSKO_QSQ) * POT_VOSKO_B1X0 - POT_T4 * ((scalar_type)2.0 / (x1 - POT_VOSKO_X0) - t1 / Xx));
			y2a = v0 + vc;
		}
		break;
    default:
     throw std::runtime_error("Error, invalid iexch");
    break;
	}
}

template<class scalar_type>
void cpu_pot(scalar_type dens, scalar_type& ex, scalar_type& ec, scalar_type& y2a, const int iexch)
{
    switch(iexch) {
    case 0: return cpu_pot_tmpl<scalar_type,0>(dens, ex, ec, y2a);
    case 1: return cpu_pot_tmpl<scalar_type,1>(dens, ex, ec, y2a);
    case 2: return cpu_pot_tmpl<scalar_type,2>(dens, ex, ec, y2a);
    case 3: return cpu_pot_tmpl<scalar_type,3>(dens, ex, ec, y2a);
    case 4: return cpu_pot_tmpl<scalar_type,4>(dens, ex, ec, y2a);
    case 5: return cpu_pot_tmpl<scalar_type,5>(dens, ex, ec, y2a);
    case 6: return cpu_pot_tmpl<scalar_type,6>(dens, ex, ec, y2a);
    case 7: return cpu_pot_tmpl<scalar_type,7>(dens, ex, ec, y2a);
    case 8: return cpu_pot_tmpl<scalar_type,8>(dens, ex, ec, y2a);
    case 9: return cpu_pot_tmpl<scalar_type,9>(dens, ex, ec, y2a);
    default: assert(false);
    }
}

#define POT_ALYP ((scalar_type)0.04918)
#define POT_BLYP ((scalar_type)0.132)
#define POT_CLYP ((scalar_type)0.2533)
#define POT_CLYP3 ((scalar_type)0.0844333333)
#define POT_DLYP ((scalar_type)0.349)
#define POT_DLYP3 ((scalar_type)0.116333333)
#define POT_CF ((scalar_type)2.87123400018819)
#define POT_BETA ((scalar_type)0.0042)


#define POT_ALF ((scalar_type)0.023266)
#define POT_BET ((scalar_type)7.389)
#define POT_GAM ((scalar_type)8.723)
#define POT_DEL ((scalar_type)0.472)

#include <float.h>

template<class scalar_type>
static void closedpbe(scalar_type rho, scalar_type agrad, scalar_type delgrad, scalar_type rlap, scalar_type& expbe, scalar_type& vxpbe, scalar_type& ecpbe, scalar_type& vcpbe);
template<class scalar_type>
static void gcorc(scalar_type rtrs, scalar_type& gg, scalar_type& grrs);

/*
 subroutine for evaluating exchange correlation density and
 potential, for non local density functionals
 4 : Exchange given by Perdew : Phys. Rev B 33 8800 (1986)
     Correlation :     Perdew : Phys. Rev B 33 8822 (1986)

 5 : Exchange given by Becke  : Phys. Rev A 38 3098 (1988)
     Correlation :     Perdew : Phys. Rev B 33 8822 (1986)
 6 : only  exchange given by Becke
 7 : BLYP, Exchange: given by Becke
           Correlation: given by LYP: PRB 37 785 (1988)
 8: Exchange given by Perdew : Phys. Rev B 33 8800 (1986)
       Correlation: given by LYP: PRB 37 785 (1988)
 to be used in concert with dnsg subroutine, the one that gives
 gradient and second derivatives of the electronic density
 19-1-93, last version: 25/11/97.
 9: PBE
*/

template<class scalar_type, int iexch>
void cpu_potg_tmpl(scalar_type dens, const vec_type<scalar_type,3>& grad, const vec_type<scalar_type,3>& hess1, const vec_type<scalar_type,3>& hess2, scalar_type& ex, scalar_type& ec, scalar_type& y2a)
{
  // hess1: xx, yy, zz
  // hess2: xy, xz, yz
  //if (dens < 1e-12) { ex = ec = 0; return; }
  if (dens < 1e-13) { ex = ec = 0; return; }

  scalar_type y = pow((scalar_type)dens, (scalar_type)0.333333333333333333);  // rho^(1/3)

  scalar_type grad2 = grad.x * grad.x + grad.y * grad.y + grad.z * grad.z;
  if (grad2 == 0) grad2 = numeric_limits<scalar_type>::min();
  scalar_type dgrad = sqrt(grad2);

  scalar_type d0 = hess1.x + hess1.y + hess1.z;
  scalar_type u0 = ((grad.x * grad.x) * hess1.x 
                 + 2.0 * grad.x * grad.y * hess2.x 
                 + 2.0 * grad.y * grad.z * hess2.z 
                 + 2.0 * grad.x * grad.z * hess2.y 
                 + (grad.y * grad.y) * hess1.y 
                 + (grad.z * grad.z) * hess1.z) / dgrad;

  y2a = 0;

  /** Exchange **/
  if (iexch == 4 || iexch == 8) {   // Perdew : Phys. Rev B 33 8800 (1986)
    scalar_type dens2 = (dens * dens);
    scalar_type ckf = (scalar_type)3.0936677 * y;
    scalar_type s = dgrad / ((scalar_type)2.0 * ckf * dens);

    scalar_type fx = (1.0 / 15.0);
    scalar_type s2 = (s * s);
    scalar_type s3 = (s * s * s);
    scalar_type g0 = 1.0 + 1.296 * s2 + 14.0 * pow(s, 4) + 0.2 * pow(s, 6);
    scalar_type F = pow(g0, fx);
    scalar_type e = POT_ALPHA * F * y;
    ex = e;

    scalar_type t = d0 / (dens * 4.0 * (ckf * ckf));
    scalar_type u = u0 / (pow(2.0 * ckf, 3) * dens2);
    //cout << t << " " << u << endl;

    scalar_type g2 = 2.592 * s + 56.0 * s3 + 1.2 * pow(s, 5);
    scalar_type g3 = 2.592 + 56.0 * s2 + 1.2 * pow(s, 4);
    scalar_type g4 = 112.0 * s + 4.8 * s3;
    scalar_type dF = fx * F/g0 * g2;
    scalar_type dsF = fx * F/g0 * (-14.0 * fx * g3 * g2/g0 + g4);
    y2a = POT_ALPHA * (1.33333333333 * F - t/s * dF - (u-1.3333333333 * s3) * dsF) * y;
  }
  else if (iexch >= 5 && iexch <= 7) { // Becke  : Phys. Rev A 38 3098 (1988)
    scalar_type e0 = POT_ALPHA * y;
    scalar_type y2 = dens / 2.0;
    scalar_type r13 = pow(y2, (scalar_type)(1.0 / 3.0));
    scalar_type r43 = pow(y2, (scalar_type)(4.0 / 3.0));
    scalar_type Xs = dgrad / (2.0 * r43);
    scalar_type siper = asinh(Xs);
    scalar_type DN = 1.0 + 6.0 * POT_BETA * Xs * siper;
    scalar_type ect = -2.0 * POT_BETA * r43 * Xs * Xs/(DN * dens);
    scalar_type e = e0 + ect;
    ex = e;

    // potential
    scalar_type v0 = 1.33333333333333 * e0;
    scalar_type Fb = 1.0 / DN;
    scalar_type XA1 = Xs / sqrt(1.0 + Xs * Xs);
    scalar_type DN1 = 1.0 + Fb * (1.0 - 6.0 * POT_BETA * Xs * XA1);
    scalar_type DN2 = 1.0 / (1.0 + Xs * Xs) + 2.0 * Fb * (2.0 - 6.0 * POT_BETA * Xs * XA1);
    scalar_type DN3 = siper * (1.0 + 2.0 * Fb) + XA1 * DN2;
    scalar_type D02 = d0 / 2.0;
    scalar_type de1 = 1.33333333333333 / (pow((scalar_type)dens,(scalar_type)2.33333333333333));
    scalar_type DGRADx = (grad.x * hess1.x + grad.y * hess2.x + grad.z * hess2.y) / dgrad;
    scalar_type GRADXx = pow(2.0, 0.33333333333333) * (1.0 / (dens * y) * DGRADx - de1 * grad.x * dgrad);
    scalar_type DGRADy = (grad.x * hess2.x + grad.y * hess1.y + grad.z * hess2.z) / dgrad;
    scalar_type GRADXy = pow(2.0, 0.33333333333333) * (1.0 / (dens * y) * DGRADy - de1 * grad.y * dgrad);
    scalar_type DGRADz = (grad.x * hess2.y + grad.y * hess2.z + grad.z * hess1.z) / dgrad;
    scalar_type GRADXz = pow(2.0, 0.33333333333333) * (1.0 / (dens * y) * DGRADz - de1 * grad.z * dgrad);
    scalar_type T1 = grad.x / 2.0 * GRADXx;
    scalar_type T2 = grad.y / 2.0 * GRADXy;
    scalar_type T3 = grad.z / 2.0 * GRADXz;
    scalar_type DN4 = 6.0 * POT_BETA * Fb * (T1 + T2 + T3);
    scalar_type DN5 = 1.33333333333333 * r43 * r13 * Xs * Xs;
    scalar_type TOT2 = DN5 - D02 * DN1 + DN4 * DN3;
    scalar_type vxc = -POT_BETA * Fb/r43 * TOT2;
    //cout << "vars: " << hess1 << endl;
    y2a = v0 + vxc;
  }
  else { // PBE (Iexch 9 complete)
    scalar_type dgrad1 = grad.x * grad.x * hess1.x;
    scalar_type dgrad2 = grad.y * grad.y * hess1.y;
    scalar_type dgrad3 = grad.z * grad.z * hess1.z;
    scalar_type dgrad4 = grad.x * grad.y * hess2.x;
    scalar_type dgrad5 = grad.x * grad.z * hess2.y;
    scalar_type dgrad6 = grad.y * grad.z * hess2.z;
    scalar_type delgrad = (dgrad1 + dgrad2 + dgrad3 + 2 * (dgrad4 + dgrad5 + dgrad6)) / dgrad;
    scalar_type rlap = hess1.x + hess1.y + hess1.z;

    scalar_type expbe, vxpbe, ecpbe, vcpbe;
    closedpbe(dens, dgrad, delgrad, rlap, expbe, vxpbe, ecpbe, vcpbe);
    ex = expbe;
    ec = ecpbe;
    y2a = vxpbe + vcpbe;
    return;
  }

  /** Correlation **/
  if (iexch >= 4 && iexch <= 6) { // Perdew : Phys. Rev B 33 8822 (1986)
    // TODO: hay algun problema con 4 y 5, probablemente este aca
    scalar_type dens2 = (dens * dens);
    scalar_type rs = POT_GL / y;
    scalar_type x1 = sqrt(rs);
    scalar_type Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
    scalar_type Xxo = (POT_VOSKO_X0 * POT_VOSKO_X0) + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1;
    scalar_type t1 = 2.0 * x1 + POT_VOSKO_B1;
    scalar_type t2 = log(Xx);
    scalar_type t3 = atan(POT_VOSKO_Q/t1);
    scalar_type t4 = POT_VOSKO_B1 * POT_VOSKO_X0/Xxo;

    ec = POT_VOSKO_A1 * (2.0 * log(x1) - t2 + 2.0 * POT_VOSKO_B1/POT_VOSKO_Q * t3 - t4 *(2.0 * log(x1 - POT_VOSKO_X0) - t2 + 2.0 * (POT_VOSKO_B1 + 2.0 * POT_VOSKO_X0)/POT_VOSKO_Q * t3));
    //cout << "vars: " << 2.0 * log(x1) << " " << t2 << " " << t3 << " " << t4 << endl;

    scalar_type t5 = (POT_VOSKO_B1 * x1 + 2.0 * POT_VOSKO_C1)/x1;
    scalar_type t6 = POT_VOSKO_X0/Xxo;
    scalar_type vc = ec - POT_VOSKO_A16 * x1 * (t5/Xx - 4.0 * POT_VOSKO_B1 / ((t1 * t1)+(POT_VOSKO_Q * POT_VOSKO_Q2)) * (1.0 - t6 * (POT_VOSKO_B1 - 2.0 * POT_VOSKO_X0)) - t4 * (2.0 / (x1 - POT_VOSKO_X0) - t1/Xx));

    if (iexch == 6) {
      y2a = y2a + vc;
    }
    else { // ?? citation??
      scalar_type rs2 = (rs * rs);
      scalar_type Cx1 = 0.002568 + POT_ALF * rs + POT_BET * rs2;
      scalar_type Cx2 = 1.0 + POT_GAM * rs + POT_DEL * rs2 + 1.0e4 * POT_BET * (rs * rs * rs);
      scalar_type C = 0.001667 + Cx1/Cx2;
      scalar_type Cx3 = POT_ALF + 2.0 * POT_BET * rs;
      scalar_type Cx4 = POT_GAM + 2.0 * POT_DEL * rs + 3.0e4 * POT_BET * rs2;
      scalar_type dC = Cx3/Cx2 - Cx1/(Cx2 * Cx2) * Cx4;
      dC = -0.333333333333333 * dC * POT_GL / (y * dens);
      scalar_type phi = 0.0008129082/C * dgrad/pow((scalar_type)dens, (scalar_type)(7.0 / 6.0));
      scalar_type expo = exp(-phi);
      scalar_type ex0 = expo * C;
      ec = ec + ex0 * grad2 / (y * dens2);

      scalar_type D1 = (2.0f - phi) * d0/dens;
      scalar_type phi2 = (phi * phi);
      scalar_type D2 = 1.33333333333333333 - 3.666666666666666666 * phi + 1.166666666666666 * phi2;
      D2 = D2 * grad2/dens2;
      scalar_type D3 = phi * (phi - 3.0) * u0/(dens * dgrad);
      scalar_type D4 = expo * grad2 / (y * dens) * (phi2 - phi - 1.0) * dC;
      vc = vc - 1.0 * (ex0 / y * (D1 - D2 + D3) - D4);
      y2a = y2a + vc;
      //cout << ec << " " << y2a << " " << D1 << " " << D2 << " " << D3 << endl;
    }
  }
  else if (iexch == 7 || iexch == 8) { // Correlation: given by LYP: PRB 37 785 (1988)
    scalar_type rom13 = pow(dens, -0.3333333333333f);
    scalar_type rom53 = pow(dens, 1.666666666666f);
    scalar_type ecro = expf(-POT_CLYP * rom13);
    scalar_type f1 = 1.0f / (1.0f + POT_DLYP * rom13);
    scalar_type tw = 1.0f / 8.0f * (grad2/dens - d0);
    scalar_type term = (tw / 9.0f + d0 / 18.0f) - 2.0f * tw + POT_CF * rom53;
    term = dens + POT_BLYP * (rom13 * rom13) * ecro * term;
    ec = -POT_ALYP * f1 * term/dens;

    // y2a
    scalar_type h1 = ecro/rom53;
    scalar_type g1 = f1 * h1;
    scalar_type tm1 = POT_DLYP3 * (rom13/dens);
    scalar_type fp1 = tm1 * (f1 * f1);
    scalar_type tm2 = -1.666666666f + POT_CLYP3 * rom13;
    scalar_type hp1 = h1 * tm2/dens;
    scalar_type gp1 = fp1 * h1 + hp1 * f1;
    scalar_type fp2 = tm1 * 2.0f * f1 * (fp1 - 0.6666666666f * f1/dens);
    scalar_type tm3 = 1.6666666666f - POT_CLYP3 * 1.3333333333f * rom13;
    scalar_type hp2 = hp1 * tm2/dens + h1 * tm3/(dens * dens);
    scalar_type gp2 = fp2 * h1 + 2.0f * fp1 * hp1 + hp2 * f1;

    scalar_type term3 = -POT_ALYP * (fp1 * dens + f1) - POT_ALYP * POT_BLYP * POT_CF * (gp1 * dens + 8.0f/3.0f * g1) * rom53;
    scalar_type term4 = (gp2 * dens * grad2 + gp1 * (3.0f * grad2 + 2.0f * dens * d0) + 4.0f * g1 * d0) * POT_ALYP * POT_BLYP/4.0f;
    scalar_type term5 = (3.0f * gp2 * dens * grad2 + gp1 * (5.0f * grad2 + 6.0f * dens * d0) + 4.0f * g1 * d0) * POT_ALYP * POT_BLYP/72.0f;
    //cout <<  term3 << " " << term4 << " " << term5 << endl;

    y2a = y2a + (term3 - term4 - term5);
  }
}

template<class scalar_type>
void cpu_potg(scalar_type dens, const vec_type<scalar_type,3>& grad, const vec_type<scalar_type,3>& hess1, const vec_type<scalar_type,3>& hess2, scalar_type& ex, scalar_type& ec, scalar_type& y2a, const int iexch)
{
    switch(iexch) {
    case 0: return cpu_potg_tmpl<scalar_type, 0>(dens, grad, hess1, hess2, ex, ec, y2a);
    case 1: return cpu_potg_tmpl<scalar_type, 1>(dens, grad, hess1, hess2, ex, ec, y2a);
    case 2: return cpu_potg_tmpl<scalar_type, 2>(dens, grad, hess1, hess2, ex, ec, y2a);
    case 3: return cpu_potg_tmpl<scalar_type, 3>(dens, grad, hess1, hess2, ex, ec, y2a);
    case 4: return cpu_potg_tmpl<scalar_type, 4>(dens, grad, hess1, hess2, ex, ec, y2a);
    case 5: return cpu_potg_tmpl<scalar_type, 5>(dens, grad, hess1, hess2, ex, ec, y2a);
    case 6: return cpu_potg_tmpl<scalar_type, 6>(dens, grad, hess1, hess2, ex, ec, y2a);
    case 7: return cpu_potg_tmpl<scalar_type, 7>(dens, grad, hess1, hess2, ex, ec, y2a);
    case 8: return cpu_potg_tmpl<scalar_type, 8>(dens, grad, hess1, hess2, ex, ec, y2a);
    case 9: return cpu_potg_tmpl<scalar_type, 9>(dens, grad, hess1, hess2, ex, ec, y2a);
    default: assert(false);
    }
}

#define CLOSEDPBE_PI32 ((scalar_type)29.608813203268075856503472999628)
#define CLOSEDPBE_AX ((scalar_type)-0.738558766382022405884230032680836)
#define CLOSEDPBE_UM ((scalar_type)0.2195149727645171)
#define CLOSEDPBE_UK ((scalar_type)0.804)
#define CLOSEDPBE_UL ((scalar_type)0.273028573090195) // um / uk
#define CLOSEDPBE_GAMMA ((scalar_type)0.03109069086965489503494086371273)
#define CLOSEDPBE_GAMMAINV ((scalar_type)32.1639684429148) // 1 / gamma
#define CLOSEDPBE_BETA ((scalar_type)0.06672455060314922)
#define CLOSEDPBE_DELTA ((scalar_type)2.14612633996736) // beta/gamma

template<class scalar_type>
static void closedpbe(scalar_type rho, scalar_type agrad, scalar_type delgrad, scalar_type rlap, scalar_type& expbe, scalar_type& vxpbe, scalar_type& ecpbe, scalar_type& vcpbe)
{
  if (rho < 2e-18) {
    expbe = vxpbe = ecpbe = vcpbe = 0;
    return;
  }

  scalar_type rho2 = rho * rho;
  scalar_type rho13 = pow((scalar_type)rho, (scalar_type)(1.0 / 3.0));
  scalar_type fk1 = pow(CLOSEDPBE_PI32, 1.0 / 3.0);
  scalar_type fk = fk1 * rho13;

  scalar_type twofk = 2.0 * fk;
  scalar_type twofk2 = pow(twofk, 2);
  scalar_type twofk3 = pow(twofk, 3);

  // S = |grad(rho)|/(2*fk*rho)
  scalar_type s = agrad / (twofk * rho);
  scalar_type s2 = pow(s, 2);
  scalar_type s3 = pow(s, 3);

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
  scalar_type p0 = 1.0 + CLOSEDPBE_UL * s2;
  scalar_type fxpbe = 1.0 + CLOSEDPBE_UK - CLOSEDPBE_UK/p0;

  // exchange pbe energy
  expbe = exlda * fxpbe;

  // Now the potential:
  scalar_type v = rlap / (twofk2 * rho);
  scalar_type u = (delgrad == 0 ? 0 : delgrad / (twofk3 * rho2));

  // Calculation of first and second derivatives
  scalar_type P2 = p0 * p0;
  scalar_type Fs = 2.0 * CLOSEDPBE_UM / P2;

  scalar_type F1 = -4.0 * CLOSEDPBE_UL * s * Fs;
  scalar_type Fss = F1/p0;

  // Now we calculate the potential Vx
  scalar_type vx2 = (4.0 / 3.0) * fxpbe;
  scalar_type vx3 = v * Fs;
  scalar_type vx4 = (u - (4.0 / 3.0) * s3) * Fss;

  vxpbe = exlda * (vx2 - vx4 - vx3);

  // Now we need to calculate the Correlation contribution
  // to the energy
  // ecpbe = eclsd*rho + h*rho
  // first we calculate the lsd contribution to the correlation energy
  // we will use the subroutine GCOR.
  // We need only the  rs (seitz radius) rs = (3/4pi*rho)^1/3
  scalar_type pirho = 4.0 * M_PI * rho;
  scalar_type rs = pow(3.0 / pirho, 1.0 / 3.0);
  scalar_type rtrs = sqrt(rs);

  scalar_type sk = sqrt(4.0 * fk / M_PI);
  scalar_type twoks = 2.0 * sk;

  scalar_type t = agrad / (twoks * rho);
  scalar_type t2 = t * t;

  scalar_type twoks2 = twoks * twoks;
  scalar_type twoks3 = twoks2 * twoks;

  scalar_type UU = (delgrad == 0 ? 0 : delgrad / (rho2 * twoks3));
  scalar_type VV = rlap / (rho * twoks2);

  scalar_type ec, eurs;
  gcorc(rtrs, ec, eurs);
	if (ec == 0) ec = numeric_limits<scalar_type>::min();

  scalar_type eclda = ec;
  scalar_type ecrs = eurs;
  scalar_type vclda = eclda - rs * (1.0 / 3.0) * ecrs;

  // Now we have to calculate the H function in order to evaluate
  // the GGA contribution to the correlation energy
  scalar_type PON = -ec * CLOSEDPBE_GAMMAINV;
  scalar_type B = CLOSEDPBE_DELTA / (exp(PON) - 1.0);
  scalar_type B2 = B * B;
  scalar_type T4 = t2 * t2;

  scalar_type Q4 = 1.0 + B * t2;
  scalar_type Q5 = 1.0 + B * t2 + B2 * T4;

  scalar_type H = (CLOSEDPBE_BETA/CLOSEDPBE_DELTA) * log(1.0 + CLOSEDPBE_DELTA * Q4 * t2/Q5);

  // So the correlation energy for pbe is:
  ecpbe = eclda + H;
  //cout << expl(PON) << " " << t2 << endl;

  // Now we have to calculate the potential contribution of GGA
  scalar_type T6 = T4 * t2;
  scalar_type RSTHRD = rs / 3.0;
  scalar_type FAC = CLOSEDPBE_DELTA / B + 1.0;
  scalar_type BEC = B2 * FAC / CLOSEDPBE_BETA;
  scalar_type Q8 = Q5 * Q5 + CLOSEDPBE_DELTA * Q4 * Q5 * t2;
  scalar_type Q9 = 1.0 + 2.0 * B * t2;
  scalar_type hB = -CLOSEDPBE_BETA * B * T6 * (2.0 + B * t2)/Q8;
  scalar_type hRS = -RSTHRD * hB * BEC * ecrs;
  scalar_type FACT0 = 2.0 * CLOSEDPBE_DELTA - 6.0 * B;
  scalar_type FACT1 = Q5 * Q9 + Q4 * Q9 * Q9;
  scalar_type hBT = 2.0 * CLOSEDPBE_BETA * T4 * ((Q4 * Q5 * FACT0 - CLOSEDPBE_DELTA * FACT1)/Q8)/Q8;
  scalar_type hRST = RSTHRD * t2 * hBT * BEC * ecrs;
  scalar_type hT = 2.0 * CLOSEDPBE_BETA * Q9/Q8;
  scalar_type FACT2 = Q4 * Q5 + B * t2 * (Q4 * Q9 + Q5);
  scalar_type FACT3 = 2.0 * B * Q5 * Q9 + CLOSEDPBE_DELTA * FACT2;
  scalar_type hTT = 4.0 * CLOSEDPBE_BETA * t * (2.0 * B/Q8 -(Q9 * FACT3 / Q8)/Q8);
  scalar_type COMM = H + hRS + hRST + t2 * hT/6.0 + 7.0 * t2 * t * hTT/6.0;

  COMM = COMM - UU * hTT - VV * hT;

  // Then, the potential for PBE is:
  vcpbe = vclda + COMM;

  //cout << expbe << " " << Q4 << " " << H << " " << eclda << " " << (double)eclda + H << " " << endl;

	//cout << rho << " " << delgrad << " " << rlap << " ret: " << expbe << " " << vxpbe << " " << ecpbe << " " << vcpbe << endl;
}

#define GCORC_A ((scalar_type)0.0310907)
#define GCORC_A1 ((scalar_type)0.21370)
#define GCORC_B1 ((scalar_type)7.5957)
#define GCORC_B2 ((scalar_type)3.5876)
#define GCORC_B3 ((scalar_type)1.6382)
#define GCORC_B4 ((scalar_type)0.49294)

template<class scalar_type>
static void gcorc(scalar_type rtrs, scalar_type& gg, scalar_type& grrs)
{
  scalar_type Q0 = -2.0 * GCORC_A * (1.0 + GCORC_A1 * rtrs * rtrs);
  scalar_type Q1 = 2.0 * GCORC_A * rtrs * (GCORC_B1 + rtrs * (GCORC_B2 + rtrs * (GCORC_B3 + GCORC_B4 * rtrs)));
  scalar_type Q2 = logf(1.0 + 1.0 / Q1);
  gg = Q0 * Q2;
  scalar_type Q3 = GCORC_A * (GCORC_B1/rtrs + 2.0 * GCORC_B2 + rtrs * (3.0 * GCORC_B3 + 4.0 * GCORC_B4 * rtrs));
  grrs = -2.0 * GCORC_A * GCORC_A1 * Q2 - Q0 * Q3/(Q1 * (1.0 + Q1));
}

template void cpu_pot(float dens, float& ex, float& ec, float& y2a, const int);
template void cpu_potg(float dens, const vec_type<float,3>& grad, const vec_type<float,3>& hess1,
                                          const vec_type<float,3>& hess2, float& ex, float& ec, float& y2a, const int);
template void cpu_pot(double dens, double& ex, double& ec, double& y2a, const int);
template void cpu_potg(double dens, const vec_type<double,3>& grad, const vec_type<double,3>& hess1,
                                          const vec_type<double,3>& hess2, double& ex, double& ec, double& y2a, const int);
}
