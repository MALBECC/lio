#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void calc_FXC(double* dens, double* tred, double vsigma, 
              double v2rho2, double v2rhosigma, double v2sigma2,
              double v3rho3, double v3rho2sigma, double v3rhosigma2,
              double v3sigma3, double* dfac, double* tfac)
{
    double cruz, contr;
    cruz  = tred[1]*dens[1] + tred[2]*dens[2] + tred[3]*dens[3];
    contr = tred[1]*tred[1] + tred[2]*tred[2] + tred[3]*tred[3];
    
    double term0;
    term0  = 2.0f * v2rho2 * tred[0] + 4.0f * v2rhosigma * cruz;
    term0 *= 2.0f; // poner tfac[0] = term0

    double derm0, derm1;
    derm0  = tred[0]*tred[0]*v3rho3 + 4.0f*cruz*tred[0]*v3rho2sigma;
    derm0 += 4.0f*cruz*cruz*v3rhosigma2 + 2.0f*contr*v2rhosigma;
    derm0 *= 4.0f; // poner dfac[0] = derm0
    
    derm1  = v3rho2sigma*tred[0]*tred[0] + 4.0f*v3rhosigma2*cruz*tred[0];
    derm1 += 4.0f*v3sigma3*cruz*cruz + 2.0f*v2sigma2*contr;
    derm1 *= 8.0f; // poner dfac[i] = derm1 * dens[i], i=1,2,3

    double temp0, temp1;
    temp0  = 4.0f*v2rhosigma*tred[0] + 8.0f*v2sigma2*cruz;
    temp0 *= 4.0f; // adicionar a dfac[i] = temp0 * tred[i], i=1,2,3

    temp1  = 4.0f*vsigma;
    temp1 *= 4.0f; // poner tfac[i] = 0.5f*(temp0*dens[i]+temp1*tred[i]), i=1,2,3


    // DENSITY FACTOR
    dfac[0] += derm0;
    dfac[1] += derm1*dens[1] + temp0*tred[1];
    dfac[2] += derm1*dens[2] + temp0*tred[2];
    dfac[3] += derm1*dens[3] + temp0*tred[3];

    // TRANSITION DENSITY FACTOR
    tfac[0] = term0;
    tfac[1] = 0.5f*(temp0*dens[1]+temp1*tred[1]);
    tfac[2] = 0.5f*(temp0*dens[2]+temp1*tred[2]);
    tfac[3] = 0.5f*(temp0*dens[3]+temp1*tred[3]);
}

