#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void calc_VXC(double* dens, double* diff, double* vrho,
              double* vsigma, double* v2rho2, double* v2rhosigma,
              double* v2sigma2, double* dfac, double* pfac, int gamma)
{
      // DIFFERENCE RELAXED DENSITY
   double C[6];

   double dd = diff[0];
   double ddx, ddy, ddz;
   double pdx, pdy, pdz;

   ddx = diff[1];
   ddy = diff[2];
   ddz = diff[3];
   pdx = dens[1] * 0.5f;
   pdy = dens[2] * 0.5f;
   pdz = dens[3] * 0.5f;

   double DUMGRV[4], DUMXX[4];
   DUMGRV[0]=ddx*pdx+ddy*pdy+ddz*pdz;
   DUMGRV[1]=ddx*pdx+ddy*pdy+ddz*pdz;
   DUMGRV[2]=ddx*pdx+ddy*pdy+ddz*pdz;
   DUMGRV[3]=ddx*pdx+ddy*pdy+ddz*pdz;

   DUMXX[0]=ddx*ddx+ddy*ddy+ddz*ddz;
   DUMXX[1]=ddx*ddx+ddy*ddy+ddz*ddz;
   DUMXX[2]=ddx*ddx+ddy*ddy+ddz*ddz;
   DUMXX[3]=DUMXX[2];

   // GROUND STATE
   double GDUMA,GDUMAG1,GDUMAG2;
   GDUMA=vrho[0]+vrho[1];
   GDUMAG1=2.0f*(vsigma[0]*2.0f+vsigma[1]);
   GDUMAG2=vsigma[1]*2.0f;

   // CONTRACTION
   double GDUMAX,GDUMAY,GDUMAZ;
   GDUMAX=GDUMAG1*pdx+GDUMAG2*pdx;
   GDUMAY=GDUMAG1*pdy+GDUMAG2*pdy;
   GDUMAZ=GDUMAG1*pdz+GDUMAG2*pdz;

   // V NON CORE CONTRIBUTION
   double DUMNV1;
   double DUMGRV1,DUMGRV3,DUMGRV4;
   DUMNV1=vrho[0]+vrho[1];
   DUMGRV1=2.0f*(vsigma[0]*2.0f+vsigma[1]);
   DUMGRV3=vsigma[1]*2.0f;
   DUMGRV4=vsigma[1]*2.0f;

   double VNCDOMA,VNCDOMAX,VNCDOMAY,VNCDOMAZ;
   double VNCDUMAX,VNCDUMAY,VNCDUMAZ;
   VNCDOMA=DUMNV1;
   VNCDOMAX=DUMGRV1*pdx+DUMGRV3*pdx;
   VNCDOMAY=DUMGRV1*pdy+DUMGRV3*pdy;
   VNCDOMAZ=DUMGRV1*pdz+DUMGRV3*pdz;

   VNCDUMAX=DUMGRV1*ddx+DUMGRV4*ddx;
   VNCDUMAY=DUMGRV1*ddy+DUMGRV4*ddy;
   VNCDUMAZ=DUMGRV1*ddz+DUMGRV4*ddz;

   // V CORE CONTRIBUTION
   C[0]=dd;
   C[1]=2.0f*DUMGRV[0];
   C[2]=DUMGRV[2];
   C[3]=dd;
   C[4]=2.0f*DUMGRV[1];
   C[5]=DUMGRV[3];

   double DUMA,DUMAG1,DUMAG2;
   DUMA=C[0]*v2rho2[0]*2.0f;
   DUMAG1=2.0f*C[0]*(v2rhosigma[0]*4.0f+v2rhosigma[1]);
   DUMAG2=C[0]*v2rhosigma[1]*2.0f;
   DUMA=DUMA+C[1]*(v2rhosigma[0]*4.0f+v2rhosigma[1]);
   DUMAG1=DUMAG1+2.0f*C[1]*(v2sigma2[0]*8.0f+v2sigma2[1]);
   DUMAG2=DUMAG2+C[1]*v2sigma2[1]*2.0f;
   DUMA=DUMA+C[2]*v2rhosigma[1]*2.0f;
   DUMAG1=DUMAG1+2.0f*C[2]*v2sigma2[1]*2.0f;
   DUMAG2=DUMAG2+C[2]*v2sigma2[1]*4.0f;
   DUMA=DUMA+C[3]*v2rho2[1]*2.0f;
   DUMAG1=DUMAG1+2.0f*C[3]*v2rhosigma[1];
   DUMAG2=DUMAG2+C[3]*v2rhosigma[1]*2.0f;
   DUMA=DUMA+C[4]*v2rhosigma[1];
   DUMAG1=DUMAG1+2.0f*C[4]*v2sigma2[1];
   DUMAG2=DUMAG2+C[4]*v2sigma2[1]*2.0f;
   DUMA=DUMA+C[5]*v2rhosigma[1]*2.0f;
   DUMAG1=DUMAG1+2.0f*C[5]*v2sigma2[1]*2.0f;
   DUMAG2=DUMAG2+C[5]*v2sigma2[1]*4.0f;

   // CONTRACTION OF VC
   double VCDUMA,VCDUMAX,VCDUMAY,VCDUMAZ;
   VCDUMA=DUMA;
   VCDUMAX=DUMAG1*pdx+DUMAG2*pdx;
   VCDUMAY=DUMAG1*pdy+DUMAG2*pdy;
   VCDUMAZ=DUMAG1*pdz+DUMAG2*pdz;
   // END V CORE

   if (gamma == 1) {
     GDUMA=0.0f;
     GDUMAX=0.0f;
     GDUMAY=0.0f;
     GDUMAZ=0.0f;
   }

   // DA
   dfac[0]=GDUMA+VCDUMA;
   dfac[1]=GDUMAX+VNCDUMAX+VCDUMAX;
   dfac[2]=GDUMAY+VNCDUMAY+VCDUMAY;
   dfac[3]=GDUMAZ+VNCDUMAZ+VCDUMAZ;

   // PA
   pfac[0]=VNCDOMA;
   pfac[1]=VNCDOMAX;
   pfac[2]=VNCDOMAY;
   pfac[3]=VNCDOMAZ;
}

