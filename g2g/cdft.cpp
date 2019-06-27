/* headers */
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <fenv.h>
#include <signal.h>
#include <cassert>
#include "common.h"
#include "init.h"
#include "matrix.h"
using namespace G2G;

/* global variables */
namespace G2G {
  CDFTVars cdft_vars;
}

extern "C" void g2g_cdft_init_(bool& do_c, bool& do_s, uint& regions,
                               uint& max_nat, unsigned int* natoms, 
                               unsigned int* at_list){
  cdft_vars.do_chrg = do_c;
  cdft_vars.do_spin = do_s;
  cdft_vars.regions = regions;

  cdft_vars.natom.resize(cdft_vars.regions);
  for (int i = 0; i < cdft_vars.regions; i++) {
    cdft_vars.natom(i) = natoms[i];
  }

  cdft_vars.atoms.resize(cdft_vars.regions, max_nat);
  // Substracts one from Frotran's indexes.
  for (int i = 0; i < cdft_vars.regions; i++) {
    for (int j = 0; j < cdft_vars.natom(i); j++) {
      cdft_vars.atoms(i,j) = at_list[i +j*max_nat] -1;
    }
  }
  cdft_vars.max_nat = max_nat;

  if (cdft_vars.do_chrg) cdft_vars.Vc.resize(cdft_vars.regions);
  if (cdft_vars.do_spin) cdft_vars.Vs.resize(cdft_vars.regions);
}

extern "C" void g2g_cdft_set_v_(double* Vc, double* Vs) {
  for (int i = 0; i < cdft_vars.regions; i++) {
    if (cdft_vars.do_chrg) cdft_vars.Vc(i) = Vc[i];
    if (cdft_vars.do_spin) cdft_vars.Vs(i) = Vs[i];
  }
}