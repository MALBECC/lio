#include <iostream>

#include <stdio.h>
#include <string.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"

#include "../libint/libintproxy.h"


using namespace G2G;

// CLOSED SHELL
extern "C" void g2g_calculate2e_(double* tao, double* fock)
{
  LIBINTproxy libintproxy;
  double fac = 0.0f;
  if ( fortran_vars.fexc != 1.0f ) {
     fac = 0.25f;
  }
     
  libintproxy.do_CoulombExchange(tao,fock,fac);
}

/*
extern "C" void g2g_exact_exchange_gradient_(double* rho, double* frc)
{
  LIBINTproxy libintproxy;
  libintproxy.do_ExchangeForces(rho,frc);
}
////////////////////////////////////////////////////////////////////////

// OPEN SHELL
extern "C" void g2g_exact_exchange_open_(double* rhoA, double* rhoB,
                                         double* fockA, double* fockB)
{
  LIBINTproxy libintproxy;
  libintproxy.do_exchange(rhoA,rhoB,fockA,fockB);
}
*/
////////////////////////////////////////////////////////////////////////

