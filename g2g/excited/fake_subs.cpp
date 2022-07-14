#include <iostream>

#include <stdlib.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"

using namespace G2G;
using namespace std;

extern "C" void g2g_calculatexc_(double* Tmat,double* Fv)
{
   cout << " " << endl;
   cout << " LIBXC was not installed. Compile LIO with libxc=1 or 2" << endl;
   cout << " " << endl;
   exit(-1);
}

extern "C" void g2g_calculateg_(double* Tmat,double* F,int& DER)
{

}

extern "C" void g2g_calcgradxc_(double* P,double* V, double* F)
{

}

namespace G2G {

void Partition::solve_lr(double* T,double* F) { }

template<class scalar_type> void PointGroupCPU<scalar_type>::
               solve_closed_lr(double* T,HostMatrix<double>& Fock) { }

template <class scalar_type>
void PointGroupCPU<scalar_type>::get_tred_input(
     HostMatrix<scalar_type>& tred_input, HostMatrix<double>& source) const { }

template<class scalar_type> void PointGroupCPU<scalar_type>::
               solve_3rd_der(double* Tmat,HostMatrix<double>& Fock,int DER) { }

template <class scalar_type> void PointGroupCPU<scalar_type>::
               solve_for_exc(double*, double*, G2G::HostMatrix<double>&, int) { }

template <class scalar_type> void PointGroupCPU<scalar_type>::recalc_densGS(
         const scalar_type*, const scalar_type*, const scalar_type*, const scalar_type*,
         HostMatrix<scalar_type>&, HostMatrix<scalar_type>&, int, HostMatrix<scalar_type>&,
         HostMatrix<scalar_type>&) { }

template <class scalar_type> void PointGroupCPU<scalar_type>::recalc_densGS3(
         const scalar_type*, const scalar_type*, const scalar_type*, const scalar_type*,
         HostMatrix<scalar_type>&, HostMatrix<scalar_type>&, HostMatrix<scalar_type>&,
         int, HostMatrix<scalar_type>&, HostMatrix<scalar_type>&, HostMatrix<scalar_type>&) { }

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}

