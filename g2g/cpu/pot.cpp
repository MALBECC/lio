/* -*- mode: c -*- */
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include "../common.h"
#include "../init.h"
#include "../cuda/cuda_extra.h"
#include "../matrix.h"
using namespace std;

#define POT_ALPHA 		-0.738558766382022447 // -(3/PI)^(1/3)
#define POT_GL 				0.620350490899400087

#define POT_VOSKO_A1 	0.03109205
#define POT_VOSKO_B1 	3.72744
#define POT_VOSKO_C1 	12.9352
#define POT_VOSKO_X0 	-0.10498

#define POT_VOSKO_Q 	6.15199066246304849
#define POT_VOSKO_A16 0.005182008333
#define POT_VOSKO_A2 	0.015546025
#define POT_VOSKO_B2 	7.06042
#define POT_VOSKO_C2 	18.0578
#define POT_VOSKO_X02 -0.32500
#define POT_VOSKO_Q2 	4.7309269
#define POT_VOSKO_A26 0.0025910042

#define POT_XX0 12.5549141492 // POT_VOSKO_X0 * POT_VOSKO_X0 + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1
#define POT_T6 -0.00836166609762834 // POT_VOSKO_X0 / POT_XX0
#define POT_T4 -0.0311676086789438 // POT_VOSKO_B1 * POT_VOSKO_X0 / POT_XX0
#define POT_VOSKO_2C1 25.8704 // 2 * POT_VOSKO_C1

#define POT_VOSKO_2B1Q 1.21178337371132 // 2 * POT_VOSKO_B1 / POT_VOSKO_Q
#define POT_VOSKO_B2X0Q 1.14352579286644 // 2 * (POT_VOSKO_B1 + 2 * POT_VOSKO_X0) / POT_VOSKO_Q
#define POT_VOSKO_4B1 14.90976 // 4.0 * POT_VOSKO_B1
#define POT_VOSKO_QSQ 37.8469891110325 // POT_VOSKO_Q * POT_VOSKO_Q
#define POT_VOSKO_B1X0 1.0329232240928 // (1.0f - t6 * (POT_VOSKO_B1 - 2.0f * POT_VOSKO_X0))

using namespace G2G;

void cpu_pot(float dens, float& ex, float& ec, float& y2a)
{
	// data X alpha

	if (dens == 0) {
		ex = 0; ec = 0;
		y2a = 0;
		return;
	}

	float y = powf(dens, 0.333333333333333333f);  // rho^(1/3)
	float v0 = -0.984745021842697f * y; // -4/3 * (3/PI)^(1/3) * rho^(1/3)

	ex = POT_ALPHA * y; // -(3/PI)^(1/3) * rho^(1/3)

	switch(fortran_vars.iexch) {
		case 1:
		{
			ec = 0;
			y2a = v0;
		}
		break;
		case 2:
		{
			float rs = POT_GL / y;
			float x1 = rs / 11.4f;
			float vc;

			if (x1 > 1.0f) {
				ec = -0.0333f * (0.5f * x1 - 0.33333333333333f);
				vc = 0.0111f * x1 * 0.5f;
			}
			else {
				float t1 = (1.0f + x1 * x1 * x1);
				float t2 = logf(1.0f + 1.0f / x1);
				float t3 = x1 * x1;
        ec = -0.0333f * (t1 * t2 - t3 + 0.5f * x1 - 0.33333333333333f);
        vc = 0.0111f * x1 * (3.0f * t3 * t2 - t1 / (x1 * (x1 + 1.0f)) - 2.0f * x1 + 0.5f);
			}
			y2a = v0 + ec + vc;
		}
		break;
		case 3:
		{
			float rs = POT_GL / y;
			float x1 = sqrtf(rs);
			float Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
			float t1 = 2.0f * x1 + POT_VOSKO_B1;
			float t2 = logf(Xx);
			float t3 = atanf(POT_VOSKO_Q/t1);
      float t5 = (POT_VOSKO_B1 * x1 + POT_VOSKO_2C1) / x1;

      ec = POT_VOSKO_A1 * (2.0f * logf(x1) - t2 + POT_VOSKO_2B1Q * t3 - POT_T4 * (2.0f * logf(x1 - POT_VOSKO_X0) - t2 + POT_VOSKO_B2X0Q * t3));

			float vc;
      vc = ec - POT_VOSKO_A16 * x1 * (t5 / Xx - POT_VOSKO_4B1 / (t1 * t1 + POT_VOSKO_QSQ) * POT_VOSKO_B1X0 - POT_T4 * (2.0f / (x1 - POT_VOSKO_X0) - t1 / Xx));
			y2a = v0 + vc;
		}
		break;
	}
}

#define POT_ALYP 0.04918f
#define POT_BLYP 0.132f
#define POT_CLYP 0.2533f
#define POT_CLYP3 0.0844333333f
#define POT_DLYP 0.349f
#define POT_DLYP3 0.116333333f
#define POT_CF 2.87123400018819f
#define POT_BETA 0.0042f


#define POT_ALF 0.023266f
#define POT_BET 7.389f
#define POT_GAM 8.723f
#define POT_DEL 0.472f

#include <float.h>

static void closedpbe(float rho, double agrad, double delgrad, double rlap, double& expbe, double& vxpbe, double& ecpbe, double& vcpbe);
static void gcorc(double rtrs, double& gg, double& grrs);

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

void cpu_potg(float dens, const double3& grad, const double3& hess1, const double3& hess2, float& ex, float& ec, float& y2a)
{
  // hess1: xx, yy, zz
  // hess2: xy, xz, yz
  //if (dens < 1e-12) { ex = ec = 0; return; }
  if (dens < 1e-13) { ex = ec = 0; return; }

  double y = pow((double)dens, 0.333333333333333333);  // rho^(1/3)

  double grad2 = grad.x * grad.x + grad.y * grad.y + grad.z * grad.z;
  if (grad2 == 0) grad2 = DBL_MIN;
  double dgrad = sqrt(grad2);

  double d0 = hess1.x + hess1.y + hess1.z;
  double u0 = ((grad.x * grad.x) * hess1.x + 2.0 * grad.x * grad.y * hess2.x + 2.0 * grad.y * grad.z * hess2.z + 2.0 * grad.x * grad.z * hess2.y +
      (grad.y * grad.y) * hess1.y + (grad.z * grad.z) * hess1.z) / dgrad; // esto ya difiere
  
  y2a = 0;

  /** Exchange **/
  if (fortran_vars.iexch == 4 || fortran_vars.iexch == 8) {   // Perdew : Phys. Rev B 33 8800 (1986)
    double dens2 = (dens * dens);
    double ckf = 3.0936677 * y;
    double s = dgrad / (2.0 * ckf * dens);

    double fx = 1.0 / 15.0;
    double s2 = (s * s);
    double s3 = (s * s * s);
    double g0 = 1.0 + 1.296 * s2 + 14.0 * pow(s, 4) + 0.2 * pow(s, 6);
    double F = pow(g0, fx);
    double e = POT_ALPHA * F * y;
    ex = e;

    double t = d0 / (dens * 4.0 * (ckf * ckf));
    double u = u0 / (pow(2.0 * ckf, 3) * dens2);
    //cout << t << " " << u << endl;

    double g2 = 2.592 * s + 56.0 * s3 + 1.2 * pow(s, 5);
    double g3 = 2.592 + 56.0 * s2 + 1.2 * pow(s, 4);
    double g4 = 112.0 * s + 4.8 * s3;
    double dF = fx * F/g0 * g2;
    double dsF = fx * F/g0 * (-14.0 * fx * g3 * g2/g0 + g4);
    y2a = POT_ALPHA * (1.33333333333 * F - t/s * dF - (u-1.3333333333 * s3) * dsF) * y;
  }
  else if (fortran_vars.iexch >= 5 && fortran_vars.iexch <= 7) { // Becke  : Phys. Rev A 38 3098 (1988)
    double e0 = POT_ALPHA * y;
    double y2 = dens / 2.0;
    double r13 = pow(y2, (1.0 / 3.0));
    double r43 = pow(y2, (4.0 / 3.0));
    double Xs = dgrad / (2.0 * r43);
    double siper = asinh(Xs);
    double DN = 1.0 + 6.0 * POT_BETA * Xs * siper;
    double ect = -2.0 * POT_BETA * r43 * Xs * Xs/(DN * dens);
    double e = e0 + ect;
    ex = e;

    // potential
    double v0 = 1.33333333333333 * e0;
    double Fb = 1.0 / DN;
    double XA1 = Xs / sqrt(1.0 + Xs * Xs);
    double DN1 = 1.0 + Fb * (1.0 - 6.0 * POT_BETA * Xs * XA1);
    double DN2 = 1.0 / (1.0 + Xs * Xs) + 2.0 * Fb * (2.0 - 6.0 * POT_BETA * Xs * XA1);
    double DN3 = siper * (1.0 + 2.0 * Fb) + XA1 * DN2;
    double D02 = d0 / 2.0;
    double de1 = 1.33333333333333 / (pow((double)dens,2.33333333333333));
    double DGRADx = (grad.x * hess1.x + grad.y * hess2.x + grad.z * hess2.y) / dgrad;
    double GRADXx = pow(2.0, 0.33333333333333) * (1.0 / (dens * y) * DGRADx - de1 * grad.x * dgrad);
    double DGRADy = (grad.x * hess2.x + grad.y * hess1.y + grad.z * hess2.z) / dgrad;
    double GRADXy = pow(2.0, 0.33333333333333) * (1.0 / (dens * y) * DGRADy - de1 * grad.y * dgrad);
    double DGRADz = (grad.x * hess2.y + grad.y * hess2.z + grad.z * hess1.z) / dgrad;
    double GRADXz = pow(2.0, 0.33333333333333) * (1.0 / (dens * y) * DGRADz - de1 * grad.z * dgrad);
    double T1 = grad.x / 2.0 * GRADXx;
    double T2 = grad.y / 2.0 * GRADXy;
    double T3 = grad.z / 2.0 * GRADXz;
    double DN4 = 6.0 * POT_BETA * Fb * (T1 + T2 + T3);
    double DN5 = 1.33333333333333 * r43 * r13 * Xs * Xs;
    double TOT2 = DN5 - D02 * DN1 + DN4 * DN3;
    double vxc = -POT_BETA * Fb/r43 * TOT2;
    //cout << "vars: " << hess1 << endl;
    y2a = v0 + vxc;
  }
  else { // PBE (Iexch 9 complete)
    double dgrad1 = grad.x * grad.x * hess1.x;
    double dgrad2 = grad.y * grad.y * hess1.y;
    double dgrad3 = grad.z * grad.z * hess1.z;
    double dgrad4 = grad.x * grad.y * hess2.x;
    double dgrad5 = grad.x * grad.z * hess2.y;
    double dgrad6 = grad.y * grad.z * hess2.z;
    double delgrad = (dgrad1 + dgrad2 + dgrad3 + 2 * (dgrad4 + dgrad5 + dgrad6)) / dgrad;
    double rlap = hess1.x + hess1.y + hess1.z;

    double expbe, vxpbe, ecpbe, vcpbe;
    closedpbe(dens, dgrad, delgrad, rlap, expbe, vxpbe, ecpbe, vcpbe);
    ex = expbe;
    ec = ecpbe;
    y2a = vxpbe + vcpbe;
    return;
  }

  /** Correlation **/
  if (fortran_vars.iexch >= 4 && fortran_vars.iexch <= 6) { // Perdew : Phys. Rev B 33 8822 (1986)
    double dens2 = (dens * dens);
    double rs = POT_GL / y;
    double x1 = sqrt(rs);
    double Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
    double Xxo = (POT_VOSKO_X0 * POT_VOSKO_X0) + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1;
    double t1 = 2.0 * x1 + POT_VOSKO_B1;
    double t2 = log(Xx);
    double t3 = atan(POT_VOSKO_Q/t1);
    double t4 = POT_VOSKO_B1 * POT_VOSKO_X0/Xxo;

    ec = POT_VOSKO_A1 * (2.0 * log(x1) - t2 + 2.0 * POT_VOSKO_B1/POT_VOSKO_Q * t3 - t4 *(2.0 * log(x1 - POT_VOSKO_X0) - t2 + 2.0 * (POT_VOSKO_B1 + 2.0 * POT_VOSKO_X0)/POT_VOSKO_Q * t3));
    //cout << "vars: " << 2.0 * log(x1) << " " << t2 << " " << t3 << " " << t4 << endl;

    double t5 = (POT_VOSKO_B1 * x1 + 2.0 * POT_VOSKO_C1)/x1;
    double t6 = POT_VOSKO_X0/Xxo;
    double vc = ec - POT_VOSKO_A16 * x1 * (t5/Xx - 4.0 * POT_VOSKO_B1 / ((t1 * t1)+(POT_VOSKO_Q * POT_VOSKO_Q2)) * (1.0 - t6 * (POT_VOSKO_B1 - 2.0 * POT_VOSKO_X0)) - t4 * (2.0 / (x1 - POT_VOSKO_X0) - t1/Xx));

    if (fortran_vars.iexch == 6) {
      y2a = y2a + vc;
    }
    else { // ?? citation??
      double rs2 = (rs * rs);
      double Cx1 = 0.002568 + POT_ALF * rs + POT_BET * rs2;
      double Cx2 = 1.0 + POT_GAM * rs + POT_DEL * rs2 + 1.0e4 * POT_BET * (rs * rs * rs);
      double C = 0.001667 + Cx1/Cx2;
      double Cx3 = POT_ALF + 2.0 * POT_BET * rs;
      double Cx4 = POT_GAM + 2.0 * POT_DEL * rs + 3.0e4 * POT_BET * rs2;
      double dC = Cx3/Cx2 - Cx1/(Cx2 * Cx2) * Cx4;
      dC = -0.333333333333333 * dC * POT_GL / (y * dens);
      double phi = 0.0008129082/C * dgrad/pow((double)dens, 7.0 / 6.0);
      double expo = exp(-phi);
      double ex0 = expo * C;
      ec = ec + ex0 * grad2 / (y * dens2);

      double D1 = (2.0f - phi) * d0/dens;
      double phi2 = (phi * phi);
      double D2 = 1.33333333333333333 - 3.666666666666666666 * phi + 1.166666666666666 * phi2;
      D2 = D2 * grad2/dens2;
      double D3 = phi * (phi - 3.0) * u0/(dens * dgrad);
      double D4 = expo * grad2 / (y * dens) * (phi2 - phi - 1.0) * dC;
      vc = vc - 1.0 * (ex0 / y * (D1 - D2 + D3) - D4);
      y2a = y2a + vc;
      //cout << ec << " " << y2a << " " << D1 << " " << D2 << " " << D3 << endl;
    }
  }
  else if (fortran_vars.iexch == 7 || fortran_vars.iexch == 8) { // Correlation: given by LYP: PRB 37 785 (1988)
    double rom13 = pow(dens, -0.3333333333333f);
    double rom53 = pow(dens, 1.666666666666f);
    double ecro = expf(-POT_CLYP * rom13);
    double f1 = 1.0f / (1.0f + POT_DLYP * rom13);
    double tw = 1.0f / 8.0f * (grad2/dens - d0);
    double term = (tw / 9.0f + d0 / 18.0f) - 2.0f * tw + POT_CF * rom53;
    term = dens + POT_BLYP * (rom13 * rom13) * ecro * term;
    ec = -POT_ALYP * f1 * term/dens;

    // y2a
    double h1 = ecro/rom53;
    double g1 = f1 * h1;
    double tm1 = POT_DLYP3 * (rom13/dens);
    double fp1 = tm1 * (f1 * f1);
    double tm2 = -1.666666666f + POT_CLYP3 * rom13;
    double hp1 = h1 * tm2/dens;
    double gp1 = fp1 * h1 + hp1 * f1;
    double fp2 = tm1 * 2.0f * f1 * (fp1 - 0.6666666666f * f1/dens);
    double tm3 = 1.6666666666f - POT_CLYP3 * 1.3333333333f * rom13;
    double hp2 = hp1 * tm2/dens + h1 * tm3/(dens * dens);
    double gp2 = fp2 * h1 + 2.0f * fp1 * hp1 + hp2 * f1;

    double term3 = -POT_ALYP * (fp1 * dens + f1) - POT_ALYP * POT_BLYP * POT_CF * (gp1 * dens + 8.0f/3.0f * g1) * rom53;
    double term4 = (gp2 * dens * grad2 + gp1 * (3.0f * grad2 + 2.0f * dens * d0) + 4.0f * g1 * d0) * POT_ALYP * POT_BLYP/4.0f;
    double term5 = (3.0f * gp2 * dens * grad2 + gp1 * (5.0f * grad2 + 6.0f * dens * d0) + 4.0f * g1 * d0) * POT_ALYP * POT_BLYP/72.0f;
    //cout <<  term3 << " " << term4 << " " << term5 << endl;

    y2a = y2a + (term3 - term4 - term5);
  }
}

#define CLOSEDPBE_PI32 29.608813203268075856503472999628
#define CLOSEDPBE_AX -0.738558766382022405884230032680836
#define CLOSEDPBE_UM 0.2195149727645171
#define CLOSEDPBE_UK 0.804
#define CLOSEDPBE_UL 0.273028573090195 // um / uk
#define CLOSEDPBE_GAMMA 0.03109069086965489503494086371273
#define CLOSEDPBE_GAMMAINV 32.1639684429148 // 1 / gamma
#define CLOSEDPBE_BETA 0.06672455060314922
#define CLOSEDPBE_DELTA 2.14612633996736 // beta/gamma

static void closedpbe(float rho, double agrad, double delgrad, double rlap, double& expbe, double& vxpbe, double& ecpbe, double& vcpbe)
{
  if (rho < 2e-18) {
    expbe = vxpbe = ecpbe = vcpbe = 0;
    return;
  }
	
  float rho2 = rho * rho;
  double rho13 = pow((double)rho, 1.0 / 3.0);
  double fk1 = pow(CLOSEDPBE_PI32, 1.0 / 3.0);
  double fk = fk1 * rho13;

  double twofk = 2.0 * fk;
  double twofk2 = pow(twofk, 2);
  double twofk3 = pow(twofk, 3);

  // S = |grad(rho)|/(2*fk*rho)
  double s = agrad / (twofk * rho);
  double s2 = pow(s, 2);
  double s3 = pow(s, 3);

  // LDA exchange contribution:
  // ex*rho ==> energy, we will calculate ex ==> energy density
  // ex*rho = -(3/4Pi)*(e^2)*(3pi)^2/3*rho^1/3*rho
  // ex*rho = -0.75*(3/Pi)^1/3*rho^4/3
  // ex*rho = ax*rho^4/3
  double exlda = CLOSEDPBE_AX * rho13;

  // In order to calculate the PBE contribution
  // to exchange energy, we have to calculate the
  // enhancement function Fx:
  // Fx = 1+uk -(uk/(1+(um*s^2)/uk)
  // um/uk = ul
  // P0 = 1 + (um*s^2)/uk
  double p0 = 1.0 + CLOSEDPBE_UL * s2;
  double fxpbe = 1.0 + CLOSEDPBE_UK - CLOSEDPBE_UK/p0;

  // exchange pbe energy
  expbe = exlda * fxpbe;

  // Now the potential:
  double v = rlap / (twofk2 * rho);
  double u = (delgrad == 0 ? 0 : delgrad / (twofk3 * rho2));

  // Calculation of first and second derivatives
  double P2 = p0 * p0;
  double Fs = 2.0 * CLOSEDPBE_UM / P2;

  double F1 = -4.0 * CLOSEDPBE_UL * s * Fs;
  double Fss = F1/p0;

  // Now we calculate the potential Vx
  double vx2 = (4.0 / 3.0) * fxpbe;
  double vx3 = v * Fs;
  double vx4 = (u - (4.0 / 3.0) * s3) * Fss;

  vxpbe = exlda * (vx2 - vx4 - vx3);

  // Now we need to calculate the Correlation contribution
  // to the energy
  // ecpbe = eclsd*rho + h*rho
  // first we calculate the lsd contribution to the correlation energy
  // we will use the subroutine GCOR.
  // We need only the  rs (seitz radius) rs = (3/4pi*rho)^1/3
  double pirho = 4.0 * M_PI * rho;
  double rs = pow(3.0 / pirho, 1.0 / 3.0);
  double rtrs = sqrt(rs);

  double sk = sqrt(4.0 * fk / M_PI);
  double twoks = 2.0 * sk;

  double t = agrad / (twoks * rho);
  double t2 = t * t;

  double twoks2 = twoks * twoks;
  double twoks3 = twoks2 * twoks;

  double UU = (delgrad == 0 ? 0 : delgrad / (rho2 * twoks3));
  double VV = rlap / (rho * twoks2);

  double ec, eurs;
  gcorc(rtrs, ec, eurs);
	if (ec == 0) ec = DBL_MIN;
	
  double eclda = ec;
  double ecrs = eurs;
  double vclda = eclda - rs * (1.0 / 3.0) * ecrs;

  // Now we have to calculate the H function in order to evaluate
  // the GGA contribution to the correlation energy
  double PON = -ec * CLOSEDPBE_GAMMAINV;
  double B = CLOSEDPBE_DELTA / (exp(PON) - 1.0);
  double B2 = B * B;
  double T4 = t2 * t2;

  double Q4 = 1.0 + B * t2;
  double Q5 = 1.0 + B * t2 + B2 * T4;

  double H = (CLOSEDPBE_BETA/CLOSEDPBE_DELTA) * log(1.0 + CLOSEDPBE_DELTA * Q4 * t2/Q5);

  // So the correlation energy for pbe is:
  ecpbe = eclda + H;
  //cout << expl(PON) << " " << t2 << endl;
	
  // Now we have to calculate the potential contribution of GGA
  double T6 = T4 * t2;
  double RSTHRD = rs / 3.0;
  double FAC = CLOSEDPBE_DELTA / B + 1.0;
  double BEC = B2 * FAC / CLOSEDPBE_BETA;
	double Q8 = Q5 * Q5 + CLOSEDPBE_DELTA * Q4 * Q5 * t2;
  double Q9 = 1.0 + 2.0 * B * t2;
  double hB = -CLOSEDPBE_BETA * B * T6 * (2.0 + B * t2)/Q8;
  double hRS = -RSTHRD * hB * BEC * ecrs;
  double FACT0 = 2.0 * CLOSEDPBE_DELTA - 6.0 * B;
  double FACT1 = Q5 * Q9 + Q4 * Q9 * Q9;
  double hBT = 2.0 * CLOSEDPBE_BETA * T4 * ((Q4 * Q5 * FACT0 - CLOSEDPBE_DELTA * FACT1)/Q8)/Q8;
  double hRST = RSTHRD * t2 * hBT * BEC * ecrs;
  double hT = 2.0 * CLOSEDPBE_BETA * Q9/Q8;
  double FACT2 = Q4 * Q5 + B * t2 * (Q4 * Q9 + Q5);
  double FACT3 = 2.0 * B * Q5 * Q9 + CLOSEDPBE_DELTA * FACT2;
  double hTT = 4.0 * CLOSEDPBE_BETA * t * (2.0 * B/Q8 -(Q9 * FACT3 / Q8)/Q8);
  double COMM = H + hRS + hRST + t2 * hT/6.0 + 7.0 * t2 * t * hTT/6.0;

  COMM = COMM - UU * hTT - VV * hT;

  // Then, the potential for PBE is:
  vcpbe = vclda + COMM;

  //cout << expbe << " " << Q4 << " " << H << " " << eclda << " " << (double)eclda + H << " " << endl;
	
	//cout << rho << " " << delgrad << " " << rlap << " ret: " << expbe << " " << vxpbe << " " << ecpbe << " " << vcpbe << endl;
}

#define GCORC_A 0.0310907
#define GCORC_A1 0.21370
#define GCORC_B1 7.5957
#define GCORC_B2 3.5876
#define GCORC_B3 1.6382
#define GCORC_B4 0.49294

static void gcorc(double rtrs, double& gg, double& grrs)
{
  double Q0 = -2.0 * GCORC_A * (1.0 + GCORC_A1 * rtrs * rtrs);
  double Q1 = 2.0 * GCORC_A * rtrs * (GCORC_B1 + rtrs * (GCORC_B2 + rtrs * (GCORC_B3 + GCORC_B4 * rtrs)));
  double Q2 = logf(1.0 + 1.0 / Q1);
  gg = Q0 * Q2;
  double Q3 = GCORC_A * (GCORC_B1/rtrs + 2.0 * GCORC_B2 + rtrs * (3.0 * GCORC_B3 + 4.0 * GCORC_B4 * rtrs));
  grrs = -2.0 * GCORC_A * GCORC_A1 * Q2 - Q0 * Q3/(Q1 * (1.0 + Q1));
}

