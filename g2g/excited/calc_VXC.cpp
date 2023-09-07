#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void calc_VXC(double* dens, double* diff, double vrho,
              double vsigma, double v2rho2, double v2rhosigma,
              double v2sigma2, double* dfac, double* pfac, int gamma)
{
   // GROUND STATE
   double gdens, gdensG;
   gdens  = vrho;
   gdensG = 2.0f*vsigma;

   if (gamma == 1) {
     gdens  = 0.0f;
     gdensG = 0.0f;
   }

   double cruz;
   cruz = dens[1]*diff[1] + dens[2]*diff[2] + dens[3]*diff[3];

   double perm0, derm0;
   perm0 = vrho;
   derm0  = v2rho2*diff[0] + 2.0f*v2rhosigma*cruz;
   derm0 *= 2.0f;

   double temp;
   temp = 2.0f*vsigma;

   double perm1;
   perm1  = v2rhosigma*diff[0] + 2.0f*v2sigma2*cruz;
   perm1 *= 4.0f;

   // Ground + Diff Excited Densities
   dfac[0] = gdens + derm0;
   dfac[1] = gdensG*dens[1] + 2.0f*temp*diff[1] + perm1*dens[1];
   dfac[2] = gdensG*dens[2] + 2.0f*temp*diff[2] + perm1*dens[2];
   dfac[3] = gdensG*dens[3] + 2.0f*temp*diff[3] + perm1*dens[3];

   // Diff Excited Density
   pfac[0] = perm0;
   pfac[1] = temp*dens[1];
   pfac[2] = temp*dens[2];
   pfac[3] = temp*dens[3];
}
