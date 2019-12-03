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
   vector<Shell>().swap(fortran_vars.obs);
   vector<int>().swap(fortran_vars.shell2atom);

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
  vector<int>().swap(fortran_vars.shell2bf);
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

int LIBINTproxy::do_exchange(double* rho, double* fock)
{

   // primero copiar rho en un objeto Matrix con las funciones d
   // en el orden de libint
   // llamar a la rutina que calcula exchange en parallel
   Matrix_E P = order_dfunc_rho(rho,fortran_vars.s_funcs,
                           fortran_vars.p_funcs,fortran_vars.d_funcs,
                           fortran_vars.m);

/*
   for(int ii=0; ii<fortran_vars.m; ii++)
      for(int jj=0;jj<fortran_vars.m; jj++)
         cout << ii << " " << jj << " " << P(ii,jj) << endl;
*/


   Matrix_E F = exchange(fortran_vars.obs, fortran_vars.m, 
                                 fortran_vars.shell2bf, P);

   order_dfunc_fock(fock,F,fortran_vars.s_funcs,
                   fortran_vars.p_funcs,fortran_vars.d_funcs,
                   fortran_vars.m);
}

Matrix_E LIBINTproxy::exchange(vector<Shell>& obs, int M, 
                      vector<int>& shell2bf, Matrix_E& D)
{

   libint2::initialize();
   using libint2::nthreads;

#pragma omp parallel
   nthreads = omp_get_num_threads();

   int nshells = obs.size();
   vector<Matrix_E> G(nthreads,Matrix_E::Zero(M,M));
   double precision = numeric_limits<double>::epsilon();
   
   // SET ENGINE LIBINT
   vector<Engine> engines(nthreads);

   engines[0] = Engine(Operator::coulomb, max_nprim(), max_l(), 0);
   engines[0].set_precision(precision);
   for(int i=1; i<nthreads; i++)
      engines[i] = engines[0];

   auto lambda = [&] (int thread_id) {
      auto& engine = engines[thread_id];
      auto& g = G[thread_id];
      const auto& buf = engine.results();

      for(int s1=0, s1234=0; s1<nshells; ++s1) {
         int bf1_first = shell2bf[s1];
         int n1 = obs[s1].size();

         for(int s2=0; s2<=s1; ++s2) {
            int bf2_first = shell2bf[s2];
            int n2 = obs[s2].size();

            for(int s3=0; s3<=s1; ++s3) {
               int bf3_first = shell2bf[s3];
               int n3 = obs[s3].size();

               int s4_max = (s1 == s3) ? s2 : s3;
               for(int s4=0; s4<=s4_max; ++s4) {
                  if( ( s1234++) % nthreads != thread_id ) continue;
                  int bf4_first = shell2bf[s4];
                  int n4 = obs[s4].size();

                  // compute the permutational degeneracy 
                  // (i.e. # of equivalents) of the given shell set
                  int s12_deg = (s1 == s2) ? 1.0 : 2.0;
                  int s34_deg = (s3 == s4) ? 1.0 : 2.0;
                  int s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
                  int s1234_deg = s12_deg * s34_deg * s12_34_deg;

                  engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                         obs[s1],obs[s2],obs[s3],obs[s4]); // testear esto no se si esta bien

                  const auto* buf_1234 = buf[0];
                  if (buf_1234 == nullptr) continue; // if all integrals screened out

                  for(int f1=0, f1234=0; f1<n1; ++f1) {
                     const int bf1 = f1 + bf1_first;
                     for(int f2=0; f2<n2; ++f2) {
                        const int bf2 = f2 + bf2_first;
                        for(int f3=0; f3<n3; ++f3) {
                           const int bf3 = f3 + bf3_first;
                           for(int f4 = 0; f4<n4; ++f4, ++f1234) {
                              const int bf4 = f4 + bf4_first;

                              const double value = buf_1234[f1234];
                              const double value_scal = value * s1234_deg;
                              g(bf1, bf3) += 0.125 * D(bf2, bf4) * value_scal;
                              g(bf2, bf4) += 0.125 * D(bf1, bf3) * value_scal;
                              g(bf1, bf4) += 0.125 * D(bf2, bf3) * value_scal;
                              g(bf2, bf3) += 0.125 * D(bf1, bf4) * value_scal;
                           }
                        }
                     }
                  } // END f...
               }
            }
         }
      } // END s...

   }; // END lambda function
   
   libint2::parallel_do(lambda);

   // accumulate contributions from all threads
   for (int i=1; i<nthreads; i++) {
       G[0] += G[i];
   }
   Matrix_E GG = 0.5f * ( G[0] + G[0].transpose() );
   return GG;
}

void LIBINTproxy::order_dfunc_fock(double* fock, Matrix_E& F,
                                   int sfunc,int pfunc,int dfunc,
                                   int M)
{
/*
 The order in d functions btween LIO and LIBINT ar differents
 LIO:    XX, XY, YY, ZX, ZY, ZZ
 LIBINT: XX, XY, ZX, YY, ZY, ZZ
*/
   int stot = sfunc;
   int ptot = pfunc*3;
   int dtot = dfunc*6;
   double ele_temp;

   if ( stot+ptot+dtot != M )
      error();

   // Copy block d and change format LIBINT->LIO
   // rows
   for(int ii=stot+ptot; ii<M; ii=ii+6) {
      for(int jj=0; jj<M; jj++) {
         ele_temp   = F(ii+2,jj);
         F(ii+2,jj) = F(ii+3,jj);
         F(ii+3,jj) = ele_temp;
      }
   }

   // cols
   for(int ii=stot+ptot; ii<M; ii=ii+6) {
      for(int jj=0; jj<M; jj++) {
         ele_temp   = F(jj,ii+2);
         F(jj,ii+2) = F(jj,ii+3);
         F(jj,ii+3) = ele_temp;
      }
   }

   for(int ii=0; ii<M; ii++) {
      fock[ii*M+ii] = F(ii,ii);
      for(int jj=0; jj<ii; jj++) {
         fock[ii*M+jj] = F(ii,jj);
         fock[jj*M+ii] = F(jj,ii);
      }
   }

}


Matrix_E LIBINTproxy::order_dfunc_rho(double* dens,int sfunc,
                         int pfunc,int dfunc,int M)
{
/*
 The order in d functions btween LIO and LIBINT ar differents
 LIO:    XX, XY, YY, ZX, ZY, ZZ
 LIBINT: XX, XY, ZX, YY, ZY, ZZ
*/
   int stot = sfunc;
   int ptot = pfunc*3;
   int dtot = dfunc*6;
   double ele_temp;

   Matrix_E DE = Matrix_E::Zero(M,M);

   if ( stot+ptot+dtot != M ) 
      error();

   // Copy All matrix
   for(int ii=0; ii<M; ii++) {
     DE(ii,ii) = dens[ii*M+ii];
     for(int jj=0; jj<ii; jj++) {
       DE(ii,jj) = dens[ii*M+jj];
       DE(jj,ii) = dens[jj*M+ii];
     }
   }
 
   // Copy block d and change format LIO->LIBINT
   // rows
   for(int ii=stot+ptot; ii<M; ii=ii+6) {
      for(int jj=0; jj<M; jj++) {
         ele_temp   = DE(ii+2,jj);
         DE(ii+2,jj)= DE(ii+3,jj);
         DE(ii+3,jj)= ele_temp;
      }
   }
   
   // cols
   for(int ii=stot+ptot; ii<M; ii=ii+6) {
      for(int jj=0; jj<M; jj++) {
         ele_temp    = DE(jj,ii+2);
         DE(jj,ii+2) = DE(jj,ii+3);
         DE(jj,ii+3) = ele_temp;
      }
   }

   return DE;
}

size_t LIBINTproxy::max_nprim() {
  size_t n = 0;
  for (auto shell:fortran_vars.obs)
    n = std::max(shell.nprim(), n);
  return n;
}

int LIBINTproxy::max_l() {
  int l = 0;
  for (auto shell:fortran_vars.obs)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}
