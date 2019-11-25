#include <iostream>

#include <stdio.h>
#include <string.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"

#include "libintproxy.h"


using namespace G2G;


extern "C" void g2g_libint_init_(double* Cbas)
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
                    fortran_vars.d_funcs);

  libintproxy.PrintBasis();
}

extern "C" void g2g_exact_exchange_(double* rho, double* fock)
{

  std::cout << "BASIS SET LIBINT" << std::endl; // print SET BASIS
  std::copy(begin(fortran_vars.obs), end(fortran_vars.obs),
         std::ostream_iterator<Shell>(std::cout, "\n"));

}
