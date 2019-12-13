#include <iostream>
#include <stdio.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"
#include "../libxc/libxcproxy.h"

using namespace G2G;
using namespace std; 

void coef_calculator(double pd,double sigma,double pdx,double pdy,double pdz,
                     double td,double tdx,double tdy,double tdz,double* COEF)
{

// LIBXC INITIALIZATION
   const int nspin = XC_POLARIZED;
   const int functionalExchange = fortran_vars.ex_functional_id;//101
   const int functionalCorrelation = fortran_vars.ec_functional_id;//130
   LibxcProxy<double,3> libxcProxy(functionalExchange, functionalCorrelation, nspin, fortran_vars.fexc);

// OUTPUTS FOR LIBXC
   double* vrho        = (double*)malloc(4*sizeof(double));
   double* vsigma      = (double*)malloc(5*sizeof(double));
   double* v2rho2      = (double*)malloc(5*sizeof(double));
   double* v2rhosigma  = (double*)malloc(8*sizeof(double));
   double* v2sigma2    = (double*)malloc(8*sizeof(double));
   double* v3rho3      = (double*)malloc(6*sizeof(double));
   double* v3rho2sigma = (double*)malloc(11*sizeof(double));
   double* v3rhosigma2 = (double*)malloc(14*sizeof(double));
   double* v3sigma3    = (double*)malloc(12*sizeof(double));

// LIBXC CALCULATE DERIVATIVES
   libxcProxy.coefZv(&pd,&sigma,
                     vrho,vsigma,
                     v2rho2,v2rhosigma,v2sigma2,
                     v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3);

   double DUMNV[2],DXV[2],DYV[2],DZV[2],DUMGRV[4],DUMXX[4];
   double C[20];

   DUMNV[0] = DUMNV[1] = td;

   DUMGRV[0]=DUMGRV[1]=tdx*pdx+tdy*pdy+tdz*pdz;
   DUMGRV[2]=DUMGRV[3]=DUMGRV[0];

   DUMXX[0]=DUMXX[1]=tdx*tdx+tdy*tdy+tdz*tdz;
   DUMXX[2]=DUMXX[3]=DUMXX[0];

   C[0]=2.0f*DUMXX[0];
   C[1]=DUMNV[0]*DUMNV[0];
   C[2]=2.0f*DUMNV[0]*DUMGRV[0];
   C[3]=2.0f*DUMGRV[0]*DUMNV[0];
   C[4]=DUMGRV[0]*DUMNV[0];
   C[5]=DUMNV[0]*DUMGRV[0];
   C[6]=4.0f*DUMGRV[0]*DUMGRV[0];
   C[7]=2.0f*DUMGRV[0]*DUMGRV[0];
   C[8]=2.0f*DUMGRV[0]*DUMGRV[0];
   C[9]=DUMGRV[0]*DUMGRV[0];

// -- EXCHANGE
   double XDUMA=0.0f;
   double XDUMAG=0.0f;
   XDUMA=XDUMA+C[0]*v2rhosigma[0];
   XDUMAG=XDUMAG+C[0]*2.0f*v2sigma2[0];
   XDUMA=XDUMA+C[1]*v3rho3[0];
   XDUMAG=XDUMAG+C[1]*2.0f*v3rho2sigma[0];
   XDUMA=XDUMA+C[2]*v3rho2sigma[0];
   XDUMAG=XDUMAG+C[2]*2.0f*v3rhosigma2[0];
   XDUMA=XDUMA+C[3]*v3rho2sigma[0];
   XDUMAG=XDUMAG+C[3]*2.0f*v3rhosigma2[0];
   XDUMA=XDUMA+C[6]*v3rhosigma2[0];
   XDUMAG=XDUMAG+C[6]*2.0f*v3sigma3[0];

// -- CORRELATION
   double CDUMA=0.0f;
   double CDUMAG1=0.0f;
   double CDUMAG2=0.0f;
   CDUMA=CDUMA+C[0]*v2rhosigma[2];
   CDUMAG1=CDUMAG1+C[0]*2.0f*v2sigma2[2];
   CDUMAG2=CDUMAG2+C[0]*v2sigma2[4];
   CDUMA=CDUMA+C[1]*v3rho3[2];
   CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma[2];
   CDUMAG2=CDUMAG2+C[1]*v3rho2sigma[4];
   CDUMA=CDUMA+C[2]*v3rho2sigma[2];
   CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2[2];
   CDUMAG2=CDUMAG2+C[2]*v3rhosigma2[4];
   CDUMA=CDUMA+C[3]*v3rho2sigma[2];
   CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2[2];
   CDUMAG2=CDUMAG2+C[3]*v3rhosigma2[4];
   CDUMA=CDUMA+C[4]*v3rho2sigma[4];
   CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2[4];
   CDUMAG2=CDUMAG2+C[4]*v3rhosigma2[7];
   CDUMA=CDUMA+C[5]*v3rho2sigma[4];
   CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2[4];
   CDUMAG2=CDUMAG2+C[5]*v3rhosigma2[7];
   CDUMA=CDUMA+C[6]*v3rhosigma2[2];
   CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3[2];
   CDUMAG2=CDUMAG2+C[6]*v3sigma3[4];
   CDUMA=CDUMA+C[7]*v3rhosigma2[4];
   CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3[4];
   CDUMAG2=CDUMAG2+C[7]*v3sigma3[7];
   CDUMA=CDUMA+C[8]*v3rhosigma2[4];
   CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3[4];
   CDUMAG2=CDUMAG2+C[8]*v3sigma3[7];
   CDUMA=CDUMA+C[9]*v3rhosigma2[7];
   CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3[7];
   CDUMAG2=CDUMAG2+C[9]*v3sigma3[11];

   C[0]=2.0f*DUMXX[1];
   C[1]=DUMNV[1]*DUMNV[1];
   C[2]=2.0f*DUMNV[1]*DUMGRV[1];
   C[3]=2.0f*DUMGRV[1]*DUMNV[1];
   C[4]=DUMGRV[1]*DUMNV[1];
   C[5]=DUMNV[1]*DUMGRV[1];
   C[6]=4.0f*DUMGRV[1]*DUMGRV[1];
   C[7]=2.0f*DUMGRV[1]*DUMGRV[1];
   C[8]=2.0f*DUMGRV[1]*DUMGRV[1];
   C[9]=DUMGRV[1]*DUMGRV[1];

   CDUMA=CDUMA+C[0]*v2rhosigma[3];
   CDUMAG1=CDUMAG1+C[0]*2.0f*v2sigma2[3];
   CDUMAG2=CDUMAG2+C[0]*v2sigma2[6];
   CDUMA=CDUMA+C[1]*v3rho3[4];
   CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma[8];
   CDUMAG2=CDUMAG2+C[1]*v3rho2sigma[10];
   CDUMA=CDUMA+C[2]*v3rho2sigma[6];
   CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2[9];
   CDUMAG2=CDUMAG2+C[2]*v3rhosigma2[12];
   CDUMA=CDUMA+C[3]*v3rho2sigma[6];
   CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2[9];
   CDUMAG2=CDUMAG2+C[3]*v3rhosigma2[12];
   CDUMA=CDUMA+C[4]*v3rho2sigma[7];
   CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2[10];
   CDUMAG2=CDUMAG2+C[4]*v3rhosigma2[13];
   CDUMA=CDUMA+C[5]*v3rho2sigma[7];
   CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2[10];
   CDUMAG2=CDUMAG2+C[5]*v3rhosigma2[13];
   CDUMA=CDUMA+C[6]*v3rhosigma2[5];
   CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3[5];
   CDUMAG2=CDUMAG2+C[6]*v3sigma3[9];
   CDUMA=CDUMA+C[7]*v3rhosigma2[6];
   CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3[6];
   CDUMAG2=CDUMAG2+C[7]*v3sigma3[10];
   CDUMA=CDUMA+C[8]*v3rhosigma2[6];
   CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3[6];
   CDUMAG2=CDUMAG2+C[8]*v3sigma3[10];
   CDUMA=CDUMA+C[9]*v3rhosigma2[7];
   CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3[7];
   CDUMAG2=CDUMAG2+C[9]*v3sigma3[11];

   C[10]=DUMXX[2];
   C[11]=DUMNV[0]*DUMNV[1];
   C[12]=2.0f*DUMNV[0]*DUMGRV[1];
   C[13]=2.0f*DUMGRV[0]*DUMNV[1];
   C[14]=DUMNV[0]*DUMGRV[3];
   C[15]=DUMGRV[2]*DUMNV[1];
   C[16]=4.0f*DUMGRV[0]*DUMGRV[1];
   C[17]=2.0f*DUMGRV[0]*DUMGRV[3];
   C[18]=2.0f*DUMGRV[2]*DUMGRV[1];
   C[19]=DUMGRV[2]*DUMGRV[3];

   CDUMA=CDUMA+C[10]*v2rhosigma[4];
   CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2[4];
   CDUMAG2=CDUMAG2+C[10]*v2sigma2[7];
   CDUMA=CDUMA+C[11]*v3rho3[3];
   CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma[5];
   CDUMAG2=CDUMAG2+C[11]*v3rho2sigma[7];
   CDUMA=CDUMA+C[12]*v3rho2sigma[3];
   CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2[3];
   CDUMAG2=CDUMAG2+C[12]*v3rhosigma2[6];
   CDUMA=CDUMA+C[13]*v3rho2sigma[5];
   CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2[8];
   CDUMAG2=CDUMAG2+C[13]*v3rhosigma2[10];
   CDUMA=CDUMA+C[14]*v3rho2sigma[4];
   CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2[4];
   CDUMAG2=CDUMAG2+C[14]*v3rhosigma2[7];
   CDUMA=CDUMA+C[15]*v3rho2sigma[7];
   CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2[10];
   CDUMAG2=CDUMAG2+C[15]*v3rhosigma2[13];
   CDUMA=CDUMA+C[16]*v3rhosigma2[3];
   CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3[3];
   CDUMAG2=CDUMAG2+C[16]*v3sigma3[6];
   CDUMA=CDUMA+C[17]*v3rhosigma2[4];
   CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3[4];
   CDUMAG2=CDUMAG2+C[17]*v3sigma3[7];
   CDUMA=CDUMA+C[18]*v3rhosigma2[6];
   CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3[6];
   CDUMAG2=CDUMAG2+C[18]*v3sigma3[10];
   CDUMA=CDUMA+C[19]*v3rhosigma2[7];
   CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3[7];
   CDUMAG2=CDUMAG2+C[19]*v3sigma3[11];

   C[10]=DUMXX[2];
   C[11]=DUMNV[0]*DUMNV[1];
   C[12]=2.0f*DUMNV[0]*DUMGRV[1];
   C[13]=2.0f*DUMGRV[0]*DUMNV[1];
   C[14]=DUMNV[0]*DUMGRV[3];
   C[15]=DUMGRV[2]*DUMNV[1];
   C[16]=4.0f*DUMGRV[0]*DUMGRV[1];
   C[17]=2.0f*DUMGRV[0]*DUMGRV[3];
   C[18]=2.0f*DUMGRV[2]*DUMGRV[1];
   C[19]=DUMGRV[2]*DUMGRV[3];

   CDUMA=CDUMA+C[10]*v2rhosigma[4];
   CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2[4];
   CDUMAG2=CDUMAG2+C[10]*v2sigma2[7];
   CDUMA=CDUMA+C[11]*v3rho3[3];
   CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma[5];
   CDUMAG2=CDUMAG2+C[11]*v3rho2sigma[7];
   CDUMA=CDUMA+C[12]*v3rho2sigma[3];
   CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2[3];
   CDUMAG2=CDUMAG2+C[12]*v3rhosigma2[6];
   CDUMA=CDUMA+C[13]*v3rho2sigma[5];
   CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2[8];
   CDUMAG2=CDUMAG2+C[13]*v3rhosigma2[10];
   CDUMA=CDUMA+C[14]*v3rho2sigma[4];
   CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2[4];
   CDUMAG2=CDUMAG2+C[14]*v3rhosigma2[7];
   CDUMA=CDUMA+C[15]*v3rho2sigma[7];
   CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2[10];
   CDUMAG2=CDUMAG2+C[15]*v3rhosigma2[13];
   CDUMA=CDUMA+C[16]*v3rhosigma2[3];
   CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3[3];
   CDUMAG2=CDUMAG2+C[16]*v3sigma3[6];
   CDUMA=CDUMA+C[17]*v3rhosigma2[4];
   CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3[4];
   CDUMAG2=CDUMAG2+C[17]*v3sigma3[7];
   CDUMA=CDUMA+C[18]*v3rhosigma2[6];
   CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3[6];
   CDUMAG2=CDUMAG2+C[18]*v3sigma3[10];
   CDUMA=CDUMA+C[19]*v3rhosigma2[7];
   CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3[7];
   CDUMAG2=CDUMAG2+C[19]*v3sigma3[11];

   C[0]=2.0f*DUMXX[0];
   C[1]=DUMNV[0]*DUMNV[0];
   C[2]=2.0f*DUMNV[0]*DUMGRV[0];
   C[3]=2.0f*DUMGRV[0]*DUMNV[0];
   C[4]=DUMGRV[0]*DUMNV[0];
   C[5]=DUMNV[0]*DUMGRV[0];
   C[6]=4.0f*DUMGRV[0]*DUMGRV[0];
   C[7]=2.0f*DUMGRV[0]*DUMGRV[0];
   C[8]=2.0f*DUMGRV[0]*DUMGRV[0];
   C[9]=DUMGRV[0]*DUMGRV[0];

// -- CORRELATION
      double CDUMB=0.0f;
      double CDUMBG1=0.0f;
      double CDUMBG2=0.0f;
      CDUMB=CDUMB+C[0]*v2rhosigma[5];
      CDUMBG1=CDUMBG1+C[0]*2.0f*v2sigma2[3];
      CDUMBG2=CDUMBG2+C[0]*v2sigma2[4];
      CDUMB=CDUMB+C[1]*v3rho3[3];
      CDUMBG1=CDUMBG1+C[1]*2.0f*v3rho2sigma[2];
      CDUMBG2=CDUMBG2+C[1]*v3rho2sigma[4];
      CDUMB=CDUMB+C[2]*v3rho2sigma[5];
      CDUMBG1=CDUMBG1+C[2]*2.0f*v3rhosigma2[3];
      CDUMBG2=CDUMBG2+C[2]*v3rhosigma2[4];
      CDUMB=CDUMB+C[3]*v3rho2sigma[5];
      CDUMBG1=CDUMBG1+C[3]*2.0f*v3rhosigma2[3];
      CDUMBG2=CDUMBG2+C[3]*v3rhosigma2[4];
      CDUMB=CDUMB+C[4]*v3rho2sigma[7];
      CDUMBG1=CDUMBG1+C[4]*2.0f*v3rhosigma2[6];
      CDUMBG2=CDUMBG2+C[4]*v3rhosigma2[7];
      CDUMB=CDUMB+C[5]*v3rho2sigma[7];
      CDUMBG1=CDUMBG1+C[5]*2.0f*v3rhosigma2[6];
      CDUMBG2=CDUMBG2+C[5]*v3rhosigma2[7];
      CDUMB=CDUMB+C[6]*v3rhosigma2[8];
      CDUMBG1=CDUMBG1+C[6]*v3sigma3[3];
      CDUMBG2=CDUMBG2+C[6]*v3sigma3[4];
      CDUMB=CDUMB+C[7]*v3rhosigma2[10];
      CDUMBG1=CDUMBG1+C[7]*2.0f*v3sigma3[6];
      CDUMBG2=CDUMBG2+C[7]*v3sigma3[7];
      CDUMB=CDUMB+C[8]*v3rhosigma2[10];
      CDUMBG1=CDUMBG1+C[8]*2.0f*v3sigma3[6];
      CDUMBG2=CDUMBG2+C[8]*2.0f*v3sigma3[7];
      CDUMB=CDUMB+C[9]*v3rhosigma2[13];
      CDUMBG1=CDUMBG1+C[9]*2.0f*v3sigma3[10];
      CDUMBG2=CDUMBG2+C[9]*v3sigma3[11];

      C[0]=2.0f*DUMXX[1];
      C[1]=DUMNV[1]*DUMNV[1];
      C[2]=2.0f*DUMNV[1]*DUMGRV[1];
      C[3]=2.0f*DUMGRV[1]*DUMNV[1];
      C[4]=DUMGRV[1]*DUMNV[1];
      C[5]=DUMNV[1]*DUMGRV[1];
      C[6]=4.0f*DUMGRV[1]*DUMGRV[1];
      C[7]=2.0f*DUMGRV[1]*DUMGRV[1];
      C[8]=2.0f*DUMGRV[1]*DUMGRV[1];
      C[9]=DUMGRV[1]*DUMGRV[1];

// -- EXCHANGE
      double XDUMB=0.0f;
      double XDUMBG=0.0f;
      XDUMB=XDUMB+C[0]*v2rhosigma[1];
      XDUMBG=XDUMBG+C[0]*v2sigma2[1];
      XDUMB=XDUMB+C[1]*v3rho3[1];
      XDUMBG=XDUMBG+C[1]*2.0f*v3rho2sigma[1];
      XDUMB=XDUMB+C[2]*v3rho2sigma[1];
      XDUMBG=XDUMBG+C[2]*2.0f*v3rhosigma2[1];
      XDUMB=XDUMB+C[3]*v3rho2sigma[1];
      XDUMBG=XDUMBG+C[3]*2.0f*v3rhosigma2[1];
      XDUMB=XDUMB+C[6]*v3rhosigma2[1];
      XDUMBG=XDUMBG+C[6]*2.0f*v3sigma3[1];
      CDUMB=CDUMB+C[0]*v2rhosigma[6];
      CDUMBG1=CDUMBG1+C[0]*2.0f*v2sigma2[5];
      CDUMBG2=CDUMBG2+C[0]*v2sigma2[6];
      CDUMB=CDUMB+C[1]*v3rho3[5];
      CDUMBG1=CDUMBG1+C[1]*2.0f*v3rho2sigma[9];
      CDUMBG2=CDUMBG2+C[1]*v3rho2sigma[10];
      CDUMB=CDUMB+C[2]*v3rho2sigma[9];
      CDUMBG1=CDUMBG1+C[2]*2.0f*v3rhosigma2[11];
      CDUMBG2=CDUMBG2+C[2]*v3rhosigma2[12];
      CDUMB=CDUMB+C[3]*v3rho2sigma[9];
      CDUMBG1=CDUMBG1+C[3]*2.0f*v3rhosigma2[11];
      CDUMBG2=CDUMBG2+C[3]*v3rhosigma2[12];
      CDUMB=CDUMB+C[4]*v3rho2sigma[10];
      CDUMBG1=CDUMBG1+C[4]*2.0f*v3rhosigma2[12];
      CDUMBG2=CDUMBG2+C[4]*v3rhosigma2[13];
      CDUMB=CDUMB+C[5]*v3rho2sigma[10];
      CDUMBG1=CDUMBG1+C[5]*2.0f*v3rhosigma2[12];
      CDUMBG2=CDUMBG2+C[5]*v3rhosigma2[13];
      CDUMB=CDUMB+C[6]*v3rhosigma2[11];
      CDUMBG1=CDUMBG1+C[6]*2.0f*v3sigma3[8];
      CDUMBG2=CDUMBG2+C[6]*v3sigma3[9];
      CDUMB=CDUMB+C[7]*v3rhosigma2[12];
      CDUMBG1=CDUMBG1+C[7]*2.0f*v3sigma3[9];
      CDUMBG2=CDUMBG2+C[7]*v3sigma3[10];
      CDUMB=CDUMB+C[8]*v3rhosigma2[12];
      CDUMBG1=CDUMBG1+C[8]*2.0f*v3sigma3[9];
      CDUMBG2=CDUMBG2+C[8]*v3sigma3[10];
      CDUMB=CDUMB+C[9]*v3rhosigma2[13];
      CDUMBG1=CDUMBG1+C[9]*2.0f*v3sigma3[10];
      CDUMBG2=CDUMBG2+C[9]*v3sigma3[11];

      C[10]=DUMXX[2];
      C[11]=DUMNV[0]*DUMNV[1];
      C[12]=2.0f*DUMNV[0]*DUMGRV[1];
      C[13]=2.0f*DUMGRV[0]*DUMNV[1];
      C[14]=DUMNV[0]*DUMGRV[3];
      C[15]=DUMGRV[2]*DUMNV[1];
      C[16]=4.0f*DUMGRV[0]*DUMGRV[1];
      C[17]=2.0f*DUMGRV[0]*DUMGRV[3];
      C[18]=2.0f*DUMGRV[2]*DUMGRV[1];
      C[19]=DUMGRV[2]*DUMGRV[3];

      CDUMB=CDUMB+C[10]*v2rhosigma[7];
      CDUMBG1=CDUMBG1+C[10]*2.0f*v2sigma2[6];
      CDUMBG2=CDUMBG2+C[10]*v2sigma2[7];
      CDUMB=CDUMB+C[11]*v3rho3[4];
      CDUMBG1=CDUMBG1+C[11]*2.0f*v3rho2sigma[6];
      CDUMBG2=CDUMBG2+C[11]*v3rho2sigma[7];
      CDUMB=CDUMB+C[12]*v3rho2sigma[6];
      CDUMBG1=CDUMBG1+C[12]*2.0f*v3rhosigma2[5];
      CDUMBG2=CDUMBG2+C[12]*v3rhosigma2[6];
      CDUMB=CDUMB+C[13]*v3rho2sigma[8];
      CDUMBG1=CDUMBG1+C[13]*2.0f*v3rhosigma2[9];
      CDUMBG2=CDUMBG2+C[13]*v3rhosigma2[10];
      CDUMB=CDUMB+C[14]*v3rho2sigma[7];
      CDUMBG1=CDUMBG1+C[14]*2.0f*v3rhosigma2[6];
      CDUMBG2=CDUMBG2+C[14]*v3rhosigma2[7];
      CDUMB=CDUMB+C[15]*v3rho2sigma[10];
      CDUMBG1=CDUMBG1+C[15]*2.0f*v3rhosigma2[12];
      CDUMBG2=CDUMBG2+C[15]*v3rhosigma2[13];
      CDUMB=CDUMB+C[16]*v3rhosigma2[9];
      CDUMBG1=CDUMBG1+C[16]*2.0f*v3sigma3[5];
      CDUMBG2=CDUMBG2+C[16]*v3sigma3[6];
      CDUMB=CDUMB+C[17]*v3rhosigma2[10];
      CDUMBG1=CDUMBG1+C[17]*2.0f*v3sigma3[6];
      CDUMBG2=CDUMBG2+C[17]*v3sigma3[7];
      CDUMB=CDUMB+C[18]*v3rhosigma2[12];
      CDUMBG1=CDUMBG1+C[18]*2.0f*v3sigma3[9];
      CDUMBG2=CDUMBG2+C[18]*v3sigma3[10];
      CDUMB=CDUMB+C[19]*v3rhosigma2[13];
      CDUMBG1=CDUMBG1+C[19]*2.0f*v3sigma3[10];
      CDUMBG2=CDUMBG2+C[19]*v3sigma3[11];

      CDUMB=CDUMB+C[10]*v2rhosigma[7];
      CDUMBG1=CDUMBG1+C[10]*2.0f*v2sigma2[6];
      CDUMBG2=CDUMBG2+C[10]*v2sigma2[7];
      CDUMB=CDUMB+C[11]*v3rho3[4];
      CDUMBG1=CDUMBG1+C[11]*2.0f*v3rho2sigma[6];
      CDUMBG2=CDUMBG2+C[11]*v3rho2sigma[7];
      CDUMB=CDUMB+C[12]*v3rho2sigma[6];
      CDUMBG1=CDUMBG1+C[12]*2.0f*v3rhosigma2[5];
      CDUMBG2=CDUMBG2+C[12]*v3rhosigma2[6];
      CDUMB=CDUMB+C[13]*v3rho2sigma[8];
      CDUMBG1=CDUMBG1+C[13]*2.0f*v3rhosigma2[9];
      CDUMBG2=CDUMBG2+C[13]*v3rhosigma2[10];
      CDUMB=CDUMB+C[14]*v3rho2sigma[7];
      CDUMBG1=CDUMBG1+C[14]*2.0f*v3rhosigma2[6];
      CDUMBG2=CDUMBG2+C[14]*v3rhosigma2[7];
      CDUMB=CDUMB+C[15]*v3rho2sigma[10];
      CDUMBG1=CDUMBG1+C[15]*2.0f*v3rhosigma2[12];
      CDUMBG2=CDUMBG2+C[15]*v3rhosigma2[13];
      CDUMB=CDUMB+C[16]*v3rhosigma2[9];
      CDUMBG1=CDUMBG1+C[16]*2.0f*v3sigma3[5];
      CDUMBG2=CDUMBG2+C[16]*v3sigma3[6];
      CDUMB=CDUMB+C[17]*v3rhosigma2[10];
      CDUMBG1=CDUMBG1+C[17]*2.0f*v3sigma3[6];
      CDUMBG2=CDUMBG2+C[17]*v3sigma3[7];
      CDUMB=CDUMB+C[18]*v3rhosigma2[12];
      CDUMBG1=CDUMBG1+C[18]*2.0f*v3sigma3[9];
      CDUMBG2=CDUMBG2+C[18]*v3sigma3[10];
      CDUMB=CDUMB+C[19]*v3rhosigma2[13];
      CDUMBG1=CDUMBG1+C[19]*2.0f*v3sigma3[10];
      CDUMBG2=CDUMBG2+C[19]*v3sigma3[11];

// --EXCHANGE
      double XDUMAGEA=0.0f;
      XDUMAGEA=XDUMAGEA+4.0f*DUMNV[0]*v2rhosigma[0];
      XDUMAGEA=XDUMAGEA+8.0f*DUMGRV[1]*v2sigma2[0];
      double XDUMAGEB=0.0f;
      XDUMAGEB=XDUMAGEB+4.0f*DUMNV[1]*v2rhosigma[1];
      XDUMAGEB=XDUMAGEB+8.0f*DUMGRV[1]*v2sigma2[1];

// --CORRELATION
      double CDUMAGEA=0.0f;
      CDUMAGEA=CDUMAGEA+4.0f*DUMNV[0]*v2rhosigma[2];
      CDUMAGEA=CDUMAGEA+8.0f*DUMGRV[0]*v2sigma2[2];
      CDUMAGEA=CDUMAGEA+4.0f*DUMGRV[2]*v2sigma2[4];
      CDUMAGEA=CDUMAGEA+4.0f*DUMNV[1]*v2rhosigma[5];
      CDUMAGEA=CDUMAGEA+8.0f*DUMGRV[1]*v2sigma2[3];
      CDUMAGEA=CDUMAGEA+4.0f*DUMGRV[3]*v2sigma2[4];
      double CDUMAGEB=0.0f;
      CDUMAGEB=CDUMAGEB+4.0f*DUMNV[0]*v2rhosigma[3];
      CDUMAGEB=CDUMAGEB+8.0f*DUMGRV[0]*v2sigma2[3];
      CDUMAGEB=CDUMAGEB+4.0f*DUMGRV[2]*v2sigma2[6];
      CDUMAGEB=CDUMAGEB+4.0f*DUMNV[1]*v2rhosigma[6];
      CDUMAGEB=CDUMAGEB+8.0f*DUMGRV[1]*v2sigma2[5];
      CDUMAGEB=CDUMAGEB+4.0f*DUMGRV[3]*v2sigma2[6];

      double CDUMAGEC=0.0f;
      CDUMAGEC=CDUMAGEC+2.0f*DUMNV[0]*v2rhosigma[4];
      CDUMAGEC=CDUMAGEC+4.0f*DUMGRV[0]*v2sigma2[4];
      CDUMAGEC=CDUMAGEC+2.0f*DUMGRV[2]*v2sigma2[7];
      CDUMAGEC=CDUMAGEC+2.0f*DUMNV[1]*v2rhosigma[7];
      CDUMAGEC=CDUMAGEC+4.0f*DUMGRV[1]*v2sigma2[6];
      CDUMAGEC=CDUMAGEC+2.0f*DUMGRV[3]*v2sigma2[7];

// --CONTRUCTION
      double DUM1A=XDUMA+CDUMA;
      double DUM2A=XDUMAG+CDUMAG1+CDUMAG2;
      double DUM3A=XDUMAGEA+CDUMAGEA+CDUMAGEC;
      COEF[0] = DUM1A;
      COEF[1] = DUM2A;
      COEF[2] = DUM3A;

// FREE MEMORY FROM LIBXC OUTPUT'S
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
