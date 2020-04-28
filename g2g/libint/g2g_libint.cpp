#include <iostream>

#include <stdio.h>
#include <string.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"

#include "libintproxy.h"


using namespace G2G;


extern "C" void g2g_libint_init_(double* Cbas, int& recalc, int& idd)
{
   
// INITIALIZATION LIBINT
   LIBINTproxy libintproxy;
   libintproxy.init( // inputs
                    fortran_vars.m, fortran_vars.atoms,
                    &fortran_vars.contractions(0),
                    Cbas, &fortran_vars.a_values(0,0),
                    &fortran_vars.atom_positions_pointer(0,0),
                    &fortran_vars.nucleii(0),
                    fortran_vars.s_funcs, fortran_vars.p_funcs,
                    fortran_vars.d_funcs, recalc, idd);

//  libintproxy.PrintBasis();
}

// CLOSED SHELL
extern "C" void g2g_exact_exchange_(double* rho, double* fock, int* op)
{
  LIBINTproxy libintproxy;
  libintproxy.do_exchange(rho,fock,op);
}

extern "C" void g2g_exact_exchange_gradient_(double* rho, double* frc, int* op)
{
  LIBINTproxy libintproxy;
  libintproxy.do_ExchangeForces(rho,frc,op);
}

// Excited State
extern "C" void g2g_calculate2e_(double* tao, double* fock, int& vecdim)
{
  LIBINTproxy libintproxy;
  libintproxy.do_CoulombExchange(tao,fock,vecdim);
}

// Exact Exchange Gradients in Excited States
extern "C" void g2g_exacgrad_excited_(double* rhoG,double* DiffExc,
                                     double* Xmat,double* fEE)
{
  LIBINTproxy libintproxy;
  libintproxy.do_ExacGradient(rhoG,DiffExc,Xmat,fEE);
}

// Gamma of Coulomb and Exchange
extern "C" void g2g_calcgammcou_(double* rhoG, double* Zmat, double* gamm)
{
  LIBINTproxy libintproxy;
  libintproxy.do_GammaCou(rhoG,Zmat,gamm);
}
////////////////////////////////////////////////////////////////////////

// OPEN SHELL
extern "C" void g2g_exact_exchange_open_(double* rhoA, double* rhoB,
                                         double* fockA, double* fockB)
{
  LIBINTproxy libintproxy;
  libintproxy.do_exchange(rhoA,rhoB,fockA,fockB);
}
////////////////////////////////////////////////////////////////////////

