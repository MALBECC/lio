#ifndef LIBINTPROXY_H
#define LIBINTPROXY_H

#include <libint2.hpp>
#include "../init.h"
#include <unordered_map>
#include <string>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix_E;

// namespace LIBINT
using libint2::Atom;
using libint2::BasisSet;
using libint2::Shell;
using libint2::Engine;
using libint2::Operator;
using libint2::BraKet;
using libint2::svector;

// namespace STD
using std::vector;
using std::string;

// Precalculated Integrals
using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
using shellpair_data_t = std::vector<std::vector<std::shared_ptr<libint2::ShellPair>>>;

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
       vector<Atom> atoms;         // atoms cordinates

       // Common Functions
       int libint_geom(double*,int);

       int make_basis(const
            vector<Atom>&,double*,double*,uint*,
            uint*,int,int,int,int);

       std::tuple<shellpair_list_t,shellpair_data_t>
       compute_shellpairs(vector<Shell>&,
                          double threshold = 1e-12);

       int map_shell();

       Matrix_E order_dfunc_rho(double*,int,int,int,int);

       void order_dfunc_fock(double*,Matrix_E&,
                             int,int,int,int);

       int error(string);

       size_t max_nprim();

       int max_l();

       // Save Integrals
       int save_ints(vector<Shell>& ,vector<int>&,int);
       long int count_ints(vector<Shell>&);
       template<Operator obtype>
       int save_calculated(vector<Shell>&, vector<int>&, double*, long int&, string);

       // Write Integrals
       int write_ints(vector<Shell>&,vector<int>&,int);
       template<Operator obtype>
       int write_calculated(vector<Shell>& ,vector<int>&, string );

       // Recalculated Integrals
       Matrix_E exchange_method(vector<Shell>&,int,vector<int>&,Matrix_E&,int*);
       template<Operator obtype>
       Matrix_E exchange(vector<Shell>&,int,vector<int>&,Matrix_E&);

       Matrix_E exchange_method_saving(vector<Shell>&,int,vector<int>&,Matrix_E&,int*);
       Matrix_E exchange_saving(vector<Shell>&,int,vector<int>&,double*,Matrix_E&,string);

       Matrix_E exchange_method_reading(vector<Shell>&,int,vector<int>&,Matrix_E&,int*);
       Matrix_E exchange_reading(vector<Shell>&,int,vector<int>&,Matrix_E&,string);

       template<Operator obtype>
       vector<Matrix_E> CoulombExchange(vector<Shell>&,int,vector<int>&,double,int,vector<Matrix_E>&);

       template<Operator obtype>
       vector<Matrix_E> CoulombExchange_saving(vector<Shell>&,int,vector<int>&,double,
                                               int,double*,vector<Matrix_E>&);

       template<Operator obtype>
       vector<Matrix_E> CoulombExchange_reading(vector<Shell>&,int,vector<int>&,double,
                                               int,vector<Matrix_E>&);

       template<Operator obtype>
       vector<Matrix_E> compute_deriv(vector<Shell>&,vector<int>&,vector<int>&,
                              int,int,Matrix_E&);

       template<Operator obtype>
       vector<Matrix_E> compute_gamma(vector<Shell>&,vector<int>&,vector<int>&,
                              int,int,Matrix_E&,double);

       void get_forcesHF(double*, vector<Matrix_E>&,
                  vector<Matrix_E>&, Matrix_E&, Matrix_E&,
                  Matrix_E&, int&, double&);

       // Open shell LR
       template<Operator obtype>
       vector<Matrix_E> CoulombExchange(vector<Shell>&,int,vector<int>&,double,int,vector<Matrix_E>&,vector<Matrix_E>&);

       template<Operator obtype>
       vector<Matrix_E> exchange(vector<Shell>&,int,vector<int>&,Matrix_E&,Matrix_E&);

       template<Operator obtype>
       vector<Matrix_E> compute_deriv(vector<Shell>&,vector<int>&,vector<int>&,
                              int,int,Matrix_E&,Matrix_E&);

public:
       // General routines
       int init(int,uint,uint*,double*,double*,
                   double*,uint*,int,int,int,int,int); // Constructor

       ~LIBINTproxy(); // Destructor

       void PrintBasis(); // Print basis in libint format

       // Closed shell
       int do_exchange(double*, double*, int*); // Energy calc.

       int do_ExchangeForces(double*, double*, int*); // Gradients calc.

       int do_CoulombExchange(double*, double*, int); // Energy calc. Excited

       // Open shell
       int do_exchange(double*, double*, double*, double*, int*); // Energy calc.

       int do_ExchangeForces(double*, double*, double*, int*); // exact exchange GS Gradients calc.

       int do_CoulombExchange(double*, double*, double*, double*, int); // Energy calc. Excited

       // Excited States Gradients with Exact Exchange
       int do_ExacGradient(double*,double*,double*,double*);

       int do_GammaCou(double*, double*, double*);

};

#endif




