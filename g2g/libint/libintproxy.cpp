#include <iostream>
#include "libintproxy.h"

using namespace G2G;
using namespace std;

namespace libint2 {
  int nthreads;
}

int LIBINTproxy::init(int M,uint natoms,uint*ncont,
                         double*cbas,double*a,double*r,uint*nuc,
                         int sfunc,int pfunc,int dfunc)
{
// LIBINT initialization
  cout << " LIBINT initialization" << endl;
  libint2::initialize();
  Shell::do_enforce_unit_normalization(false);

  int err;

// Put coordinates in a object
  // set atoms object
  err = libint_geom(r,natoms);
  if ( err != 0 ) error();

// Basis in LIBINT format
  // this set obs and shell2atom objects
  err = make_basis(atoms,a,cbas,ncont,nuc,sfunc,pfunc,dfunc,M);
  if ( err != 0 ) error();

// First basis in a shell and atom centre in the shell
  // this set shell2bf
  err = map_shell();
  if ( err != 0 ) error();

  return 0;

}

int LIBINTproxy::error()
{
  cout << "Something is wrong in " << "caca" << endl;
  exit(-1);
}

LIBINTproxy::~LIBINTproxy()
{
// LIBINT DEINITIALIZATION
  libint2::finalize();
  vector<Atom>().swap(atoms);
}

int LIBINTproxy::libint_geom(double* r,int natoms)
{
   atoms.resize(natoms);
   int err = 0;

   for(int i=0; i<natoms; i++) {
     atoms[i].x = r[i];
     atoms[i].y = r[i+natoms];
     atoms[i].z = r[i+natoms*2];
   }

   if ( atoms.empty() == true ) 
      err = -1;
   else
      err = 0;

   return err;
}

void LIBINTproxy::PrintBasis()
{
// Check basis normalization with libint
  std::cout << "BASIS SET LIBINT" << std::endl; // print SET BASIS
  std::copy(begin(fortran_vars.obs), end(fortran_vars.obs),
         std::ostream_iterator<Shell>(std::cout, "\n"));
}

int LIBINTproxy::make_basis(
     const vector<Atom>& atoms,double*a,double*c,
     uint*ncont,uint*nuc,int s_func,int p_func,int d_func,int M)
{
   int err = 0;
   int from = 0;
   int to = s_func;

// FOR S FUNCTIONS
   for(int i=from; i<to; i++) {
     int centro = nuc[i]-1;
     int tam = ncont[i];
     vector<double> exp(tam);
     vector<double> coef(tam);
     for(int cont=0; cont<tam; cont++) {
        exp[cont] = a[i+M*cont];
        coef[cont] = c[i+M*cont];
     }
     fortran_vars.obs.push_back(
        {
          exp,
          { 
            {0, false, coef}
          },
          {{atoms[centro].x,atoms[centro].y,atoms[centro].z}}
        }
     );
     fortran_vars.shell2atom.push_back(centro);
     vector<double>().swap(exp);
     vector<double>().swap(coef);
   }

// FOR P FUNCTIONS
   from = s_func;
   to = s_func+p_func*3;
   for(int i=from; i<to; i=i+3) {
      int centro = nuc[i]-1;
      int tam = ncont[i];
      vector<double> exp(tam);
      vector<double> coef(tam);
      for(int cont=0; cont<tam; cont++) {
         exp[cont] = a[i+M*cont];
         coef[cont] = c[i+M*cont];
      }
      fortran_vars.obs.push_back(
         {
           exp,
           {
             {1, false, coef}
           },
           {{atoms[centro].x,atoms[centro].y,atoms[centro].z}}
         }
      );
      fortran_vars.shell2atom.push_back(centro);
      vector<double>().swap(exp);
      vector<double>().swap(coef);
   }

// FOR D FUNCTIONS
   from = s_func+p_func*3;
   to = M;
   for(int i=from; i<to; i=i+6) {
      int centro = nuc[i]-1;
      int tam = ncont[i];
      vector<double> exp(tam);
      vector<double> coef(tam);
      for(int cont=0; cont<tam; cont++) {
         exp[cont] = a[i+M*cont];
         coef[cont] = c[i+M*cont];
      }
      fortran_vars.obs.push_back(
         {
           exp,
           {
             {2, false, coef}
           },
           {{atoms[centro].x,atoms[centro].y,atoms[centro].z}}
         }
      );
      fortran_vars.shell2atom.push_back(centro);
      vector<double>().swap(exp);
      vector<double>().swap(coef);
    }

   err = -1;
   if ( fortran_vars.obs.empty() == false && fortran_vars.shell2atom.empty() == false ) 
      err = 0;

   return err;
}

int LIBINTproxy::map_shell()
{
  int err = 0;
  fortran_vars.shell2bf.reserve(fortran_vars.obs.size());

  int n = 0;
  for (auto shell: fortran_vars.obs) {
    fortran_vars.shell2bf.push_back(n);
    n += shell.size();
  }
  
  err = -1;
  if ( fortran_vars.shell2bf.empty() == false ) 
     err = 0;

  return err;
  
}

