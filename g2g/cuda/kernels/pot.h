/*
 * function for evaluating exchange correlation density and
 * potential, depends of Iexch, index that says which potential
 * will be used
 * - 1 X alpha
 * - 2 Gunnarson-Lundqvist
 * - 3 Vosko, Wilk and Nusair
 * Exchange part is in all cases the same, for the time being LD
 * Self interaction corrections are used in correlation part
 */

/* pot_kernel constants */
#define POT_ALPHA 		-0.738558766382022447f // -(3/PI)^(1/3)
#define POT_GL 				0.620350490899400087f

#define POT_VOSKO_A1 	0.03109205f
#define POT_VOSKO_B1 	3.72744f
#define POT_VOSKO_C1 	12.9352f
#define POT_VOSKO_X0 	-0.10498f

#define POT_VOSKO_Q 	6.15199066246304849f
#define POT_VOSKO_A16 0.005182008333f
#define POT_VOSKO_A2 	0.015546025f
#define POT_VOSKO_B2 	7.06042f
#define POT_VOSKO_C2 	18.0578f
#define POT_VOSKO_X02 -0.32500f
#define POT_VOSKO_Q2 	4.7309269f
#define POT_VOSKO_A26 0.0025910042f

#define POT_XX0 12.5549141492f // POT_VOSKO_X0 * POT_VOSKO_X0 + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1
#define POT_T6 -0.00836166609762834f // POT_VOSKO_X0 / POT_XX0
#define POT_T4 -0.0311676086789438f // POT_VOSKO_B1 * POT_VOSKO_X0 / POT_XX0
#define POT_VOSKO_2C1 25.8704f // 2 * POT_VOSKO_C1

#define POT_VOSKO_2B1Q 1.21178337371132f // 2 * POT_VOSKO_B1 / POT_VOSKO_Q
#define POT_VOSKO_B2X0Q 1.14352579286644f // 2 * (POT_VOSKO_B1 + 2 * POT_VOSKO_X0) / POT_VOSKO_Q
#define POT_VOSKO_QSQ 37.8469891110325f // POT_VOSKO_Q * POT_VOSKO_Q
#define POT_VOSKO_4B1B1X0 15.4006373696499f // 4.0 * POT_VOSKO_B1 * (1.0f - t6 * (POT_VOSKO_B1 - 2.0f * POT_VOSKO_X0))

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

__device__ void closedpbe(float rho, float agrad, float delgrad, float rlap, float& expbe, float& vxpbe, float& ecpbe, float& vcpbe);
__device__ void gcorc(float rtrs, float& gg, float& grrs);

#include <float.h>

	
template<bool compute_exc, bool compute_y2a, bool lda> __device__ void gpu_pot(float dens, const float4& grad, const float4& hess1, const float4& hess2, float& exc_corr, float& y2a)
{
	// data X alpha
  float ec;

  if (lda) {
    if (dens == 0.0f) {
      if (compute_exc) { exc_corr = 0.0f; }
      if (compute_y2a) y2a = 0.0f;
      return;
    }

    float y = powf(dens, 0.333333333333333333f);  // rho^(1/3)
    float v0 = -0.984745021842697f * y; // -4/3 * (3/PI)^(1/3) * rho^(1/3)

    if (compute_exc) exc_corr = POT_ALPHA * y; // -(3/PI)^(1/3) * rho^(1/3)

    switch(gpu_Iexch) {
      case 1:
      {
        if (compute_y2a) y2a = v0;
      }
      break;
      case 2:
      {
        float rs = POT_GL / y;
        float x1 = rs / 11.4f;
        float vc;

        if (x1 > 1.0f) {
          ec = -0.0333f * (0.5f * x1 - 0.33333333333333f);
          if (compute_y2a) vc = 0.0111f * x1 * 0.5f;
        }
        else {
          float t1 = (1.0f + x1 * x1 * x1);
          float t2 = logf(1.0f + 1.0f / x1);
          float t3 = x1 * x1;
          ec = -0.0333f * (t1 * t2 - t3 + 0.5f * x1 - 0.33333333333333f);
          if (compute_y2a) vc = 0.0111f * x1 * (3.0f * t3 * t2 - t1 / (x1 * (x1 + 1.0f)) - 2.0f * x1 + 0.5f);
        }
        if (compute_y2a) y2a = v0 + ec + vc;
        if (compute_exc) exc_corr += ec;
      }
      break;
      case 3:
      {
        float rs = POT_GL / y;
        float x1 = sqrtf(rs);
        float Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
        float t1 = 2.0f * x1 + POT_VOSKO_B1;
        float t2 = logf(Xx);

        ec = POT_VOSKO_A1 * (2.0f * logf(x1) - t2 + POT_VOSKO_2B1Q * atanf(POT_VOSKO_Q/t1)
          - POT_T4 * (2.0f * logf(x1 - POT_VOSKO_X0) - t2 + POT_VOSKO_B2X0Q * atanf(POT_VOSKO_Q/t1)));

        if (compute_y2a) {
          float vc;
          vc = ec - POT_VOSKO_A16 * x1 * (((POT_VOSKO_B1 * x1 + POT_VOSKO_2C1) / (x1 * Xx)) -
            POT_VOSKO_4B1B1X0 / (t1 * t1 + POT_VOSKO_QSQ) - POT_T4 * (2.0f / (x1 - POT_VOSKO_X0) - t1 / Xx));
          y2a = v0 + vc;
        }
        if (compute_exc) exc_corr += ec;
      }
      break;
    }
  }
  /*** GGA ***/
  else {
    // hess1: xx, yy, zz
    // hess2: xy, xz, yz
    if (dens < 1e-13f) { exc_corr = ec = 0.0f; return; }

    float y = powf(dens, 0.333333333333333333f);  // rho^(1/3)

    float grad2 = grad.x * grad.x + grad.y * grad.y + grad.z * grad.z;
    if (grad2 == 0.0f) grad2 = FLT_MIN;
    float dgrad = sqrtf(grad2);

    float d0 = hess1.x + hess1.y + hess1.z;
    float u0 = ((grad.x * grad.x) * hess1.x + 2.0f * grad.x * grad.y * hess2.x + 2.0f * grad.y * grad.z * hess2.z + 2.0f * grad.x * grad.z * hess2.y +
        (grad.y * grad.y) * hess1.y + (grad.z * grad.z) * hess1.z) / dgrad;

    y2a = 0.0f;

    /** Exchange **/
    if (gpu_Iexch == 4 || gpu_Iexch == 8) {   // Perdew : Phys. Rev B 33 8800 (1986)
      float dens2 = (dens * dens);
      float ckf = 3.0936677f * y;
      float s = dgrad / (2.0f * ckf * dens);

      float fx = 1.0f / 15.0f;
      float s2 = (s * s);
      float s3 = (s * s * s);
      float g0 = 1.0f + 1.296f * s2 + 14.0f * powf(s, 4.0f) + 0.2f * powf(s, 6.0f);
      float F = powf(g0, fx);
      float e = POT_ALPHA * F * y;
      exc_corr = e;

      float t = d0 / (dens * 4.0f * (ckf * ckf));
      float u = u0 / (powf(2.0f * ckf, 3.0f) * dens2);
      //cout << t << " " << u << endl;

      float g2 = 2.592f * s + 56.0f * s3 + 1.2f * powf(s, 5.0f);
      float g3 = 2.592f + 56.0f * s2 + 1.2f * powf(s, 4.0f);
      float g4 = 112.0f * s + 4.8f * s3;
      float dF = fx * F/g0 * g2;
      float dsF = fx * F/g0 * (-14.0f * fx * g3 * g2/g0 + g4);
      y2a = POT_ALPHA * (1.33333333333f * F - t/s * dF - (u-1.3333333333f * s3) * dsF) * y;
    }
    else if (gpu_Iexch >= 5 && gpu_Iexch <= 7) { // Becke  : Phys. Rev A 38 3098 (1988)
      float e0 = POT_ALPHA * y;
      float y2 = dens / 2.0f;
      float r13 = powf(y2, (1.0f / 3.0f));
      float r43 = powf(y2, (4.0f / 3.0f));
      float Xs = dgrad / (2.0f * r43);
      float siper = asinhf(Xs);
      float DN = 1.0f + 6.0f * POT_BETA * Xs * siper;
      float ect = -2.0f * POT_BETA * r43 * Xs * Xs/(DN * dens);
      float e = e0 + ect;
      exc_corr = e;

      // potential
      float v0 = 1.33333333333333f * e0;
      float Fb = 1.0f / DN;
      float XA1 = Xs / sqrtf(1.0f + Xs * Xs);
      float DN1 = 1.0f + Fb * (1.0f - 6.0f * POT_BETA * Xs * XA1);
      float DN2 = 1.0f / (1.0f + Xs * Xs) + 2.0f * Fb * (2.0f - 6.0f * POT_BETA * Xs * XA1);
      float DN3 = siper * (1.0f + 2.0f * Fb) + XA1 * DN2;
      float D02 = d0 / 2.0f;
      float de1 = 1.33333333333333f / (powf(dens, 2.33333333333333f));
      float DGRADx = (grad.x * hess1.x + grad.y * hess2.x + grad.z * hess2.y) / dgrad;
      float GRADXx = powf(2.0f, 0.33333333333333f) * (1.0f / (dens * y) * DGRADx - de1 * grad.x * dgrad);
      float DGRADy = (grad.x * hess2.x + grad.y * hess1.y + grad.z * hess2.z) / dgrad;
      float GRADXy = powf(2.0f, 0.33333333333333f) * (1.0f / (dens * y) * DGRADy - de1 * grad.y * dgrad);
      float DGRADz = (grad.x * hess2.y + grad.y * hess2.z + grad.z * hess1.z) / dgrad;
      float GRADXz = powf(2.0f, 0.33333333333333f) * (1.0f / (dens * y) * DGRADz - de1 * grad.z * dgrad);
      float T1 = grad.x / 2.0f * GRADXx;
      float T2 = grad.y / 2.0f * GRADXy;
      float T3 = grad.z / 2.0f * GRADXz;
      float DN4 = 6.0f * POT_BETA * Fb * (T1 + T2 + T3);
      float DN5 = 1.33333333333333f * r43 * r13 * Xs * Xs;
      float TOT2 = DN5 - D02 * DN1 + DN4 * DN3;
      float vxc = -POT_BETA * Fb/r43 * TOT2;
      //cout << "vars: " << hess1 << endl;
      y2a = v0 + vxc;
    }
    else { // PBE (Iexch 9 complete)
      float dgrad1 = grad.x * grad.x * hess1.x;
      float dgrad2 = grad.y * grad.y * hess1.y;
      float dgrad3 = grad.z * grad.z * hess1.z;
      float dgrad4 = grad.x * grad.y * hess2.x;
      float dgrad5 = grad.x * grad.z * hess2.y;
      float dgrad6 = grad.y * grad.z * hess2.z;
      float delgrad = (dgrad1 + dgrad2 + dgrad3 + 2.0f * (dgrad4 + dgrad5 + dgrad6)) / dgrad;
      float rlap = hess1.x + hess1.y + hess1.z;

      float expbe, vxpbe, ecpbe, vcpbe;
      closedpbe(dens, dgrad, delgrad, rlap, expbe, vxpbe, ecpbe, vcpbe);
      exc_corr = expbe;
      ec = ecpbe;
      y2a = vxpbe + vcpbe;
      return;
    }

    /** Correlation **/
    if (gpu_Iexch >= 4 && gpu_Iexch <= 6) { // Perdew : Phys. Rev B 33 8822 (1986)
      // TODO: hay algun problema con 4 y 5, probablemente este aca
      float dens2 = (dens * dens);
      float rs = POT_GL / y;
      float x1 = sqrtf(rs);
      float Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
      float Xxo = (POT_VOSKO_X0 * POT_VOSKO_X0) + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1;
      float t1 = 2.0f * x1 + POT_VOSKO_B1;
      float t2 = logf(Xx);
      float t3 = atanf(POT_VOSKO_Q/t1);
      float t4 = POT_VOSKO_B1 * POT_VOSKO_X0/Xxo;

      ec = POT_VOSKO_A1 * (2.0f * logf(x1) - t2 + 2.0f * POT_VOSKO_B1/POT_VOSKO_Q * t3 - t4 *(2.0f * logf(x1 - POT_VOSKO_X0) - t2 + 2.0f * (POT_VOSKO_B1 + 2.0f * POT_VOSKO_X0)/POT_VOSKO_Q * t3));
      //cout << "vars: " << 2.0 * log(x1) << " " << t2 << " " << t3 << " " << t4 << endl;

      float t5 = (POT_VOSKO_B1 * x1 + 2.0f * POT_VOSKO_C1)/x1;
      float t6 = POT_VOSKO_X0/Xxo;
      float vc = ec - POT_VOSKO_A16 * x1 * (t5/Xx - 4.0f * POT_VOSKO_B1 / ((t1 * t1)+(POT_VOSKO_Q * POT_VOSKO_Q2)) * (1.0f - t6 * (POT_VOSKO_B1 - 2.0f * POT_VOSKO_X0)) - t4 * (2.0f / (x1 - POT_VOSKO_X0) - t1/Xx));

      if (gpu_Iexch == 6) {
        y2a = y2a + vc;
      }
      else { // ?? citation??
        float rs2 = (rs * rs);
        float Cx1 = 0.002568f + POT_ALF * rs + POT_BET * rs2;
        float Cx2 = 1.0f + POT_GAM * rs + POT_DEL * rs2 + 1.0e4f * POT_BET * (rs * rs * rs);
        float C = 0.001667f + Cx1/Cx2;
        float Cx3 = POT_ALF + 2.0f * POT_BET * rs;
        float Cx4 = POT_GAM + 2.0f * POT_DEL * rs + 3.0e4f * POT_BET * rs2;
        float dC = Cx3/Cx2 - Cx1/(Cx2 * Cx2) * Cx4;
        dC = -0.333333333333333f * dC * POT_GL / (y * dens);
        float phi = 0.0008129082f/C * dgrad/powf(dens, 7.0f / 6.0f);
        float expo = exp(-phi);
        float ex0 = expo * C;
        ec = ec + ex0 * grad2 / (y * dens2);

        float D1 = (2.0f - phi) * d0/dens;
        float phi2 = (phi * phi);
        float D2 = 1.33333333333333333f - 3.666666666666666666f * phi + 1.166666666666666f * phi2;
        D2 = D2 * grad2/dens2;
        float D3 = phi * (phi - 3.0f) * u0/(dens * dgrad);
        float D4 = expo * grad2 / (y * dens) * (phi2 - phi - 1.0f) * dC;
        vc = vc - 1.0f * (ex0 / y * (D1 - D2 + D3) - D4);
        y2a = y2a + vc;
        //cout << ec << " " << y2a << " " << D1 << " " << D2 << " " << D3 << endl;
      }
    }
    else if (gpu_Iexch == 7 || gpu_Iexch == 8) { // Correlation: given by LYP: PRB 37 785 (1988)
      float rom13 = powf(dens, -0.3333333333333f);
      float rom53 = powf(dens, 1.666666666666f);
      float ecro = expf(-POT_CLYP * rom13);
      float f1 = 1.0f / (1.0f + POT_DLYP * rom13);
      float tw = 1.0f / 8.0f * (grad2/dens - d0);
      float term = (tw / 9.0f + d0 / 18.0f) - 2.0f * tw + POT_CF * rom53;
      term = dens + POT_BLYP * (rom13 * rom13) * ecro * term;
      ec = -POT_ALYP * f1 * term/dens;

      // y2a
      float h1 = ecro/rom53;
      float g1 = f1 * h1;
      float tm1 = POT_DLYP3 * (rom13/dens);
      float fp1 = tm1 * (f1 * f1);
      float tm2 = -1.666666666f + POT_CLYP3 * rom13;
      float hp1 = h1 * tm2/dens;
      float gp1 = fp1 * h1 + hp1 * f1;
      float fp2 = tm1 * 2.0f * f1 * (fp1 - 0.6666666666f * f1/dens);
      float tm3 = 1.6666666666f - POT_CLYP3 * 1.3333333333f * rom13;
      float hp2 = hp1 * tm2/dens + h1 * tm3/(dens * dens);
      float gp2 = fp2 * h1 + 2.0f * fp1 * hp1 + hp2 * f1;

      float term3 = -POT_ALYP * (fp1 * dens + f1) - POT_ALYP * POT_BLYP * POT_CF * (gp1 * dens + 8.0f/3.0f * g1) * rom53;
      float term4 = (gp2 * dens * grad2 + gp1 * (3.0f * grad2 + 2.0f * dens * d0) + 4.0f * g1 * d0) * POT_ALYP * POT_BLYP/4.0f;
      float term5 = (3.0f * gp2 * dens * grad2 + gp1 * (5.0f * grad2 + 6.0f * dens * d0) + 4.0f * g1 * d0) * POT_ALYP * POT_BLYP/72.0f;
      //cout <<  term3 << " " << term4 << " " << term5 << endl;

      y2a = y2a + (term3 - term4 - term5);
    }
  }
}

#define CLOSEDPBE_PI32 29.608813203268075856503472999628f
#define CLOSEDPBE_AX -0.738558766382022405884230032680836f
#define CLOSEDPBE_UM 0.2195149727645171f
#define CLOSEDPBE_UK 0.804f
#define CLOSEDPBE_UL 0.273028573090195f // um / uk
#define CLOSEDPBE_GAMMA 0.03109069086965489503494086371273f
#define CLOSEDPBE_GAMMAINV 32.1639684429148f // 1 / gamma
#define CLOSEDPBE_BETA 0.06672455060314922f
#define CLOSEDPBE_DELTA 2.14612633996736f // beta/gamma

__device__ void closedpbe(float rho, float agrad, float delgrad, float rlap, float& expbe, float& vxpbe, float& ecpbe, float& vcpbe)
{
  if (rho < 2e-18f) {
    expbe = vxpbe = ecpbe = vcpbe = 0.0f;
    return;
  }

  float rho2 = rho * rho;
  float rho13 = powf(rho, 1.0f / 3.0f);
  float fk1 = powf(CLOSEDPBE_PI32, 1.0f / 3.0f);
  float fk = fk1 * rho13;

  float twofk = 2.0f * fk;
  float twofk2 = powf(twofk, 2.0f);
  float twofk3 = powf(twofk, 3.0f);

  // S = |grad(rho)|/(2*fk*rho)
  float s = agrad / (twofk * rho);
  float s2 = powf(s, 2.0f);
  float s3 = powf(s, 3.0f);

  // LDA exchange contribution:
  // ex*rho ==> energy, we will calculate ex ==> energy density
  // ex*rho = -(3/4Pi)*(e^2)*(3pi)^2/3*rho^1/3*rho
  // ex*rho = -0.75*(3/Pi)^1/3*rho^4/3
  // ex*rho = ax*rho^4/3
  float exlda = CLOSEDPBE_AX * rho13;

  // In order to calculate the PBE contribution
  // to exchange energy, we have to calculate the
  // enhancement function Fx:
  // Fx = 1+uk -(uk/(1+(um*s^2)/uk)
  // um/uk = ul
  // P0 = 1 + (um*s^2)/uk
  float p0 = 1.0f + CLOSEDPBE_UL * s2;
  float fxpbe = 1.0f + CLOSEDPBE_UK - CLOSEDPBE_UK/p0;

  // exchange pbe energy
  expbe = exlda * fxpbe;

  // Now the potential:
  float v = rlap / (twofk2 * rho);
  float u = (delgrad == 0.0f ? 0.0f : delgrad / (twofk3 * rho2));

  // Calculation of first and second derivatives
  float P2 = p0 * p0;
  float Fs = 2.0f * CLOSEDPBE_UM / P2;

  float F1 = -4.0f * CLOSEDPBE_UL * s * Fs;
  float Fss = F1/p0;

  // Now we calculate the potential Vx
  float vx2 = (4.0f / 3.0f) * fxpbe;
  float vx3 = v * Fs;
  float vx4 = (u - (4.0f / 3.0f) * s3) * Fss;

  vxpbe = exlda * (vx2 - vx4 - vx3);

  // Now we need to calculate the Correlation contribution
  // to the energy
  // ecpbe = eclsd*rho + h*rho
  // first we calculate the lsd contribution to the correlation energy
  // we will use the subroutine GCOR.
  // We need only the  rs (seitz radius) rs = (3/4pi*rho)^1/3
  float pirho = 4.0f * CUDART_PI_F * rho;
  float rs = powf(3.0f / pirho, 1.0f / 3.0f);
  float rtrs = sqrtf(rs);

  float sk = sqrtf(4.0f * fk / CUDART_PI_F);
  float twoks = 2.0f * sk;

  float t = agrad / (twoks * rho);
  float t2 = t * t;

  float twoks2 = twoks * twoks;
  float twoks3 = twoks2 * twoks;

  float UU = (delgrad == 0.0f ? 0.0f : delgrad / (rho2 * twoks3));
  float VV = rlap / (rho * twoks2);

  float ec, eurs;
  gcorc(rtrs, ec, eurs);
	if (ec == 0.0f) ec = FLT_MIN;

  float eclda = ec;
  float ecrs = eurs;
  float vclda = eclda - rs * (1.0f / 3.0f) * ecrs;

  // Now we have to calculate the H function in order to evaluate
  // the GGA contribution to the correlation energy
  float PON = -ec * CLOSEDPBE_GAMMAINV;
  float B = CLOSEDPBE_DELTA / (expf(PON) - 1.0f);
  float B2 = B * B;
  float T4 = t2 * t2;

  float Q4 = 1.0f + B * t2;
  float Q5 = 1.0f + B * t2 + B2 * T4;

  float H = (CLOSEDPBE_BETA/CLOSEDPBE_DELTA) * logf(1.0f + CLOSEDPBE_DELTA * Q4 * t2/Q5);

  // So the correlation energy for pbe is:
  ecpbe = eclda + H;
  //cout << expl(PON) << " " << t2 << endl;

  // Now we have to calculate the potential contribution of GGA
  float T6 = T4 * t2;
  float RSTHRD = rs / 3.0f;
  float FAC = CLOSEDPBE_DELTA / B + 1.0f;
  float BEC = B2 * FAC / CLOSEDPBE_BETA;
	float Q8 = Q5 * Q5 + CLOSEDPBE_DELTA * Q4 * Q5 * t2;
  float Q9 = 1.0f + 2.0f * B * t2;
  float hB = -CLOSEDPBE_BETA * B * T6 * (2.0f + B * t2)/Q8;
  float hRS = -RSTHRD * hB * BEC * ecrs;
  float FACT0 = 2.0f * CLOSEDPBE_DELTA - 6.0f * B;
  float FACT1 = Q5 * Q9 + Q4 * Q9 * Q9;
  float hBT = 2.0f * CLOSEDPBE_BETA * T4 * ((Q4 * Q5 * FACT0 - CLOSEDPBE_DELTA * FACT1)/Q8)/Q8;
  float hRST = RSTHRD * t2 * hBT * BEC * ecrs;
  float hT = 2.0f * CLOSEDPBE_BETA * Q9/Q8;
  float FACT2 = Q4 * Q5 + B * t2 * (Q4 * Q9 + Q5);
  float FACT3 = 2.0f * B * Q5 * Q9 + CLOSEDPBE_DELTA * FACT2;
  float hTT = 4.0f * CLOSEDPBE_BETA * t * (2.0f * B/Q8 -(Q9 * FACT3 / Q8)/Q8);
  float COMM = H + hRS + hRST + t2 * hT/6.0f + 7.0f * t2 * t * hTT/6.0f;

  COMM = COMM - UU * hTT - VV * hT;

  // Then, the potential for PBE is:
  vcpbe = vclda + COMM;

  //cout << expbe << " " << Q4 << " " << H << " " << eclda << " " << (float)eclda + H << " " << endl;

	//cout << rho << " " << delgrad << " " << rlap << " ret: " << expbe << " " << vxpbe << " " << ecpbe << " " << vcpbe << endl;
}

#define GCORC_A 0.0310907f
#define GCORC_A1 0.21370f
#define GCORC_B1 7.5957f
#define GCORC_B2 3.5876f
#define GCORC_B3 1.6382f
#define GCORC_B4 0.49294f

__device__ void gcorc(float rtrs, float& gg, float& grrs)
{
  float Q0 = -2.0f * GCORC_A * (1.0f + GCORC_A1 * rtrs * rtrs);
  float Q1 = 2.0f * GCORC_A * rtrs * (GCORC_B1 + rtrs * (GCORC_B2 + rtrs * (GCORC_B3 + GCORC_B4 * rtrs)));
  float Q2 = logf(1.0f + 1.0f / Q1);
  gg = Q0 * Q2;
  float Q3 = GCORC_A * (GCORC_B1/rtrs + 2.0f * GCORC_B2 + rtrs * (3.0f * GCORC_B3 + 4.0f * GCORC_B4 * rtrs));
  grrs = -2.0f * GCORC_A * GCORC_A1 * Q2 - Q0 * Q3/(Q1 * (1.0f + Q1));
}

