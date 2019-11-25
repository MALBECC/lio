#ifndef LIBINTPROXY_H
#define LIBINTPROXY_H

#include <libint2.hpp>
#include "../init.h"

// namespace LIBINT
using libint2::Atom;
using libint2::BasisSet;
using libint2::Shell;
using libint2::Engine;
using libint2::Operator;
using libint2::BraKet;

// namespace STD
using std::vector;

typedef unsigned int uint;

// NEEDED TO RUN LIBINT IN PARALLEL
namespace libint2 {
extern int nthreads;

/// fires off \c nthreads instances of lambda in parallel
template <typename Lambda>
void parallel_do(Lambda& lambda) {
#ifdef _OPENMP
#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#endif
}
}

class LIBINTproxy
{
private:
       // Variables
       vector<libint2::Shell> obs; // Basis (in libint format)
       vector<int> shell2bf;       // first basis function
       vector<Atom> atoms;         // atoms cordinates
       vector<int> shell2atom;     // atom centre of shell

       // Functions
       int libint_geom(double*,int);

       int make_basis(const
            vector<Atom>&,double*,double*,uint*,
            uint*,int,int,int,int);

       int map_shell();

       int error( );

public:
       LIBINTproxy(int,uint,uint*,double*,double*,
                   double*,uint*,int,int,int); // Constructor

       ~LIBINTproxy(); // Destructor

};

#endif




