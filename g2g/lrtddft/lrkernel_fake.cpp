//This is a fake subroutine in case of you run a Linear Response
//calculate without LIBXC and LIBINT or when you compiled with 
// cuda /= 0;
#include <iostream>
#include <stdio.h>
#include <cstdlib>

#include "../init.h"
#include "../partition.h"

using namespace G2G;
extern Partition partition;

//######################################################################
//######################################################################

extern "C" void g2g_linear_response_(double* MatCoef,double* KMat,
                                     double* K_int,double* Cbas,int& dim,
                                     int& NCOlr, int& Nvirt)
{
   printf("Linear Response only works with\n");
   printf("cuda=0 libxc=1 libint=1\n");
   exit(-1);

   fflush(stdout); // NOT BUFFERED
}
namespace G2G {

void Partition::solve_lr(double* K,double* calcK,double* Coef)
{
}
template<class scalar_type> void PointGroupGPU<scalar_type>::
       solve_closed_lr(double* gIntKfxc)
{
}
template<class scalar_type>
void PointGroupCPU<scalar_type>::solve_closed_lr(double* gIntK)
{
}
template<class scalar_type>
void PointGroupCPU<scalar_type>::get_coef_input(HostMatrix<scalar_type>&
               rmm_input,int* vecnum) const
{
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#if GPU_KERNELS
template class PointGroupGPU<double>;
#endif
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#if GPU_KERNELS
template class PointGroupGPU<float>;
#endif
#endif
}
