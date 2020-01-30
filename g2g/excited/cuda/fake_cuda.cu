#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math_constants.h>
#include <string>
#include <vector>
#include <stdio.h>
#include "../../init.h"
#include "../../common.h"
#include "../../partition.h"
#include "../../timer.h"

using namespace std;

namespace G2G {

template<class scalar_type>
void PointGroupGPU<scalar_type>::lr_closed_init() { }

template<class scalar_type>
void PointGroupGPU<scalar_type>::solve_closed_lr(double* T,HostMatrix<double>& Fock) { }

template <class scalar_type>
void PointGroupGPU<scalar_type>::get_tred_input(
     HostMatrix<scalar_type>& tred_input, HostMatrix<double>& source) const { }

template<class scalar_type>
 void PointGroupGPU<scalar_type>::
               solve_3rd_der(double* Tmat,HostMatrix<double>& Fock,int& DER) { }

template<class scalar_type> void PointGroupGPU<scalar_type>::
        solve_for_exc(double*P,double*V,HostMatrix<double>& F) { }


#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupGPU<double>;
#else
template class PointGroup<float>;
template class PointGroupGPU<float>;
#endif
}

