#include <iostream>

#include <stdlib.h>

using namespace std;


extern "C" void g2g_libint_init_(double* Cbas)
{
  cout << "  Please compile LIO with option: libint=1" << endl;
  exit(-1);
}

extern "C" void g2g_exact_exchange_(double* rho, double* fock)
{
  cout << "  Please compile LIO with option: libint=1" << endl;
  exit(-1);
}

extern "C" void g2g_exact_exchange_gradient_(double* rho, double* frc)
{
  cout << "  Please compile LIO with option: libint=1" << endl;
  exit(-1);
}

extern "C" void g2g_exact_exchange_open_(double* rhoA, double* fockA,
                                         double* rhoB, double* fockB)
{
  cout << "  Please compile LIO with option: libint=1" << endl;
  exit(-1);
}


