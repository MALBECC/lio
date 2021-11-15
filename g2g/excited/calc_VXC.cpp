#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void calc_VXC(double* dens, double* diff, double* vrho, double* vsigma,
              double* v2rho2, double* v2rhosigma, double* v2sigma2,
              double* dfac, double* pfac, int gamma) {
  // DIFFERENCE RELAXED DENSITY
  double coef[6];

  double dd = diff[0];
  double diff_derX, diff_derY, diff_derZ;
  double dens_derX, dens_derY, dens_derZ;

  diff_derX = diff[1];
  diff_derY = diff[2];
  diff_derZ = diff[3];
  dens_derX = dens[1] * 0.5f;
  dens_derY = dens[2] * 0.5f;
  dens_derZ = dens[3] * 0.5f;

  double grad =
      diff_derX * dens_derX + diff_derY * dens_derY + diff_derZ * dens_derZ;

  // GROUND STATE
  double gdens, gdensG1, gdensG2;
  gdens = vrho[0] + vrho[1];
  gdensG1 = 2.0f * (vsigma[0] * 2.0f + vsigma[1]);
  gdensG2 = vsigma[1] * 2.0f;

  // CONTRACTION
  double gdensX, gdensY, gdensZ;
  gdensX = gdensG1 * dens_derX + gdensG2 * dens_derX;
  gdensY = gdensG1 * dens_derY + gdensG2 * dens_derY;
  gdensZ = gdensG1 * dens_derZ + gdensG2 * dens_derZ;

  // V NON CORE CONTRIBUTION
  double dncV1;
  double grad1, grad2, grad3;
  dncV1 = vrho[0] + vrho[1];
  grad1 = 2.0f * (vsigma[0] * 2.0f + vsigma[1]);
  grad2 = vsigma[1] * 2.0f;
  grad3 = vsigma[1] * 2.0f;

  double Vncdens, VncdensX, VncdensY, VncdensZ;
  double VncdiffX, VncdiffY, VncdiffZ;
  Vncdens = dncV1;
  VncdensX = grad1 * dens_derX + grad2 * dens_derX;
  VncdensY = grad1 * dens_derY + grad2 * dens_derY;
  VncdensZ = grad1 * dens_derZ + grad2 * dens_derZ;

  VncdiffX = grad1 * diff_derX + grad3 * diff_derX;
  VncdiffY = grad1 * diff_derY + grad3 * diff_derY;
  VncdiffZ = grad1 * diff_derZ + grad3 * diff_derZ;

  // V CORE CONTRIBUTION
  coef[0] = dd;
  coef[1] = 2.0f * grad;
  coef[2] = grad;
  coef[3] = dd;
  coef[4] = 2.0f * grad;
  coef[5] = grad;

  double term1, term2, term3;
  term1 = coef[0] * v2rho2[0] * 2.0f;
  term2 = 2.0f * coef[0] * (v2rhosigma[0] * 4.0f + v2rhosigma[1]);
  term3 = coef[0] * v2rhosigma[1] * 2.0f;
  term1 = term1 + coef[1] * (v2rhosigma[0] * 4.0f + v2rhosigma[1]);
  term2 = term2 + 2.0f * coef[1] * (v2sigma2[0] * 8.0f + v2sigma2[1]);
  term3 = term3 + coef[1] * v2sigma2[1] * 2.0f;
  term1 = term1 + coef[2] * v2rhosigma[1] * 2.0f;
  term2 = term2 + 2.0f * coef[2] * v2sigma2[1] * 2.0f;
  term3 = term3 + coef[2] * v2sigma2[1] * 4.0f;
  term1 = term1 + coef[3] * v2rho2[1] * 2.0f;
  term2 = term2 + 2.0f * coef[3] * v2rhosigma[1];
  term3 = term3 + coef[3] * v2rhosigma[1] * 2.0f;
  term1 = term1 + coef[4] * v2rhosigma[1];
  term2 = term2 + 2.0f * coef[4] * v2sigma2[1];
  term3 = term3 + coef[4] * v2sigma2[1] * 2.0f;
  term1 = term1 + coef[5] * v2rhosigma[1] * 2.0f;
  term2 = term2 + 2.0f * coef[5] * v2sigma2[1] * 2.0f;
  term3 = term3 + coef[5] * v2sigma2[1] * 4.0f;

  // CONTRACTION OF VC
  double Vcdiff, VcdiffX, VcdiffY, VcdiffZ;
  Vcdiff = term1;
  VcdiffX = term2 * dens_derX + term3 * dens_derX;
  VcdiffY = term2 * dens_derY + term3 * dens_derY;
  VcdiffZ = term2 * dens_derZ + term3 * dens_derZ;
  // END V CORE

  if (gamma == 1) {
    gdens = 0.0f;
    gdensX = 0.0f;
    gdensY = 0.0f;
    gdensZ = 0.0f;
  }

  // DA
  dfac[0] = gdens + Vcdiff;
  dfac[1] = gdensX + VncdiffX + VcdiffX;
  dfac[2] = gdensY + VncdiffY + VcdiffY;
  dfac[3] = gdensZ + VncdiffZ + VcdiffZ;

  // PA
  pfac[0] = Vncdens;
  pfac[1] = VncdensX;
  pfac[2] = VncdensY;
  pfac[3] = VncdensZ;
}
