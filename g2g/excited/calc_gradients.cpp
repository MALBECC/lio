#include <iostream>
#include <stdio.h>
#include <string.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"
#include "../libxc/libxcproxy.h"

#include "calc_VXC.h"
#include "calc_FXC.h"

using namespace G2G;
using namespace std;

void calc_gradients(double* dens, double* diff, double* trad,
                    double sigma, double* dfac, double* pfac,
                    double* tfac, int calc_fxc)
{
   memset(dfac,0.0f,4*sizeof(double));
   memset(pfac,0.0f,4*sizeof(double));
   memset(tfac,0.0f,4*sizeof(double));

// LIBXC INITIALIZATION
   const int nspin = XC_UNPOLARIZED;
   const int functionalExchange = fortran_vars.ex_functional_id; //101;
   const int functionalCorrelation = fortran_vars.ec_functional_id; // 130;
   LibxcProxy<double,3> libxcProxy(functionalExchange, functionalCorrelation, nspin, fortran_vars.fexc);

// OUTPUTS FOR LIBXC: 0 = exchange; 1 = correlation
   double* vrho        = (double*)malloc(2*sizeof(double));
   double* vsigma      = (double*)malloc(2*sizeof(double));
   double* v2rho2      = (double*)malloc(2*sizeof(double));
   double* v2rhosigma  = (double*)malloc(2*sizeof(double));
   double* v2sigma2    = (double*)malloc(2*sizeof(double));
   double* v3rho3      = (double*)malloc(2*sizeof(double));
   double* v3rho2sigma = (double*)malloc(2*sizeof(double));
   double* v3rhosigma2 = (double*)malloc(2*sizeof(double));
   double* v3sigma3    = (double*)malloc(2*sizeof(double));

#define libxc_parameter \
   dens, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, \
   v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3
   libxcProxy.terms_derivs(libxc_parameter);
#undef libxc_parameter

#define VXC_parameter \
   dens, diff, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, \
   dfac, pfac, calc_fxc
   calc_VXC(VXC_parameter);
#undef VXC_parameter

   if ( calc_fxc == 0 ) {
   #define FXC_parameter \
      dens, trad, vsigma, v2rho2, v2rhosigma, v2sigma2, \
      v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3, dfac, tfac
      calc_FXC(FXC_parameter);
   #undef FXC_parameter
   }

   free(vrho); vrho = NULL;
   free(vsigma); vsigma = NULL;
   free(v2rho2); v2rho2 = NULL;
   free(v2rhosigma); v2rhosigma = NULL;
   free(v2sigma2); v2sigma2 = NULL;
   free(v3rho3); v3rho3 = NULL;
   free(v3rho2sigma); v3rho2sigma = NULL;
   free(v3rhosigma2); v3rhosigma2 = NULL;
   free(v3sigma3); v3sigma3 = NULL;
}

