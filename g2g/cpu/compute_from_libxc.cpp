#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../common.h"
#include "../init.h"
#include "../matrix.h"
#include "../partition.h"
using namespace std;

namespace G2G {
template <class scalar_type>
void PointGroupCPU<scalar_type>::compute_rmm_libxc(const uint& group_m, const scalar_type* fv,
               const scalar_type* gxv, const scalar_type* gyv, const scalar_type* gzv,
               const scalar_type& wp, double* coef_a, double* coef_b,
               const G2G::vec_type<scalar_type, 3>& dxyz_a, const G2G::vec_type<scalar_type, 3>& dxyz_b, 
               double* smallFock_a, double* smallFock_b)
{
   double term, term_a, term_b, result;
   for(int i=0; i<group_m; i++) {
     term   = fv[i] * fv[i];
     term_a = 2.0f * fv[i] * ( gxv[i]*dxyz_a.x + gyv[i]*dxyz_a.y + gzv[i]*dxyz_a.z );
     term_b = 2.0f * fv[i] * ( gxv[i]*dxyz_b.x + gyv[i]*dxyz_b.y + gzv[i]*dxyz_b.z );

     // alpha
     result = coef_a[0] * term + 2.0f * coef_a[1] * term_a + coef_a[2] * term_b;
     smallFock_a[i*group_m+i] += result * wp;
     // beta
     result = coef_b[0] * term + 2.0f * coef_b[1] * term_b + coef_b[2] * term_a;
     smallFock_b[i*group_m+i] += result * wp;

     for(int j=0; j<i; j++) {
        term    = fv[i] * fv[j];
        term_a  = fv[j] * ( gxv[i]*dxyz_a.x + gyv[i]*dxyz_a.y + gzv[i]*dxyz_a.z );
        term_a += fv[i] * ( gxv[j]*dxyz_a.x + gyv[j]*dxyz_a.y + gzv[j]*dxyz_a.z );
        term_b  = fv[j] * ( gxv[i]*dxyz_b.x + gyv[i]*dxyz_b.y + gzv[i]*dxyz_b.z );
        term_b += fv[i] * ( gxv[j]*dxyz_b.x + gyv[j]*dxyz_b.y + gzv[j]*dxyz_b.z );

        // alpha
        result = coef_a[0] * term + 2.0f * coef_a[1] * term_a + coef_a[2] * term_b;
        smallFock_a[i*group_m+j] += result * wp;

        // beta
        result = coef_b[0] * term + 2.0f * coef_b[1] * term_b + coef_b[2] * term_a;
        smallFock_b[i*group_m+j] += result * wp;
     } // end j function
   } // end i function
}

template <class scalar_type>
void PointGroupCPU<scalar_type>::compute_forces_libxc(const uint& group_m, const scalar_type& wp, int& local_atoms,
               const scalar_type* fv, const scalar_type* gxv, const scalar_type* gyv, const scalar_type* gzv,
               const scalar_type* hpxv, const scalar_type* hpyv, const scalar_type* hpzv,
               const scalar_type* hixv, const scalar_type* hiyv, const scalar_type* hizv,
               HostMatrix<scalar_type>& rmm_input_a, HostMatrix<scalar_type>& rmm_input_b,
               const G2G::vec_type<scalar_type, 3>& dxyz_a, const G2G::vec_type<scalar_type, 3>& dxyz_b,
               double* coef_a, double* coef_b,
               HostMatrix<scalar_type>& ddx_a, HostMatrix<scalar_type>& ddy_a, HostMatrix<scalar_type>& ddz_a,
               HostMatrix<scalar_type>& ddx_b, HostMatrix<scalar_type>& ddy_b, HostMatrix<scalar_type>& ddz_b,
               double* smallFor_a, double* smallFor_b)
{
   ddx_a.zero(); ddy_a.zero(); ddz_a.zero();
   ddx_b.zero(); ddy_b.zero(); ddz_b.zero();
   double temp[3], Pa, Pb;

   // cantidad de bases, s+p+d (sin tener en cuenta las 3p y 5d
   for (int i = 0, ii = 0; i < this->total_functions_simple(); i++) {
      uint nuc = this->func2local_nuc(ii);// da a que nucleo LOCAL pertenece la base ii
      uint inc_i = this->small_function_type(i);// da cuantas funciones tiene ese tipo: s->1, p->3, d->5
      double grdx_a, grdy_a, grdz_a; grdx_a = grdy_a = grdz_a = 0.0f;
      double grdx_b, grdy_b, grdz_b; grdx_b = grdy_b = grdz_b = 0.0f;

      for (uint k = 0; k < inc_i; k++, ii++) {
         for (uint j = 0; j < group_m; j++) {
            double factor = (ii == j ? 2.0f : 1.0f);
            Pa = rmm_input_a(j,ii) * factor;
            Pb = rmm_input_b(j,ii) * factor;

            // alpha
            temp[0] = 1.0f * coef_a[0] * Pa;
            temp[1] = 2.0f * coef_a[1] * Pa;
            temp[2] = 1.0f * coef_a[2] * Pa;

            grdx_a += temp[0]*gxv[ii]*fv[j];
            grdx_a += temp[1]*( dxyz_a.x*hpxv[ii]*fv[j] + dxyz_a.x*gxv[ii]*gxv[j] );
            grdx_a += temp[1]*( dxyz_a.y*hixv[ii]*fv[j] + dxyz_a.y*gxv[ii]*gyv[j] );
            grdx_a += temp[1]*( dxyz_a.z*hiyv[ii]*fv[j] + dxyz_a.z*gxv[ii]*gzv[j] );
            grdx_a += temp[2]*( dxyz_b.x*hpxv[ii]*fv[j] + dxyz_b.x*gxv[ii]*gxv[j] );
            grdx_a += temp[2]*( dxyz_b.y*hixv[ii]*fv[j] + dxyz_b.y*gxv[ii]*gyv[j] );
            grdx_a += temp[2]*( dxyz_b.z*hiyv[ii]*fv[j] + dxyz_b.z*gxv[ii]*gzv[j] );

            grdy_a += temp[0]*gyv[ii]*fv[j];
            grdy_a += temp[1]*( dxyz_a.x*hixv[ii]*fv[j] + dxyz_a.x*gyv[ii]*gxv[j] );
            grdy_a += temp[1]*( dxyz_a.y*hpyv[ii]*fv[j] + dxyz_a.y*gyv[ii]*gyv[j] );
            grdy_a += temp[1]*( dxyz_a.z*hizv[ii]*fv[j] + dxyz_a.z*gyv[ii]*gzv[j] );
            grdy_a += temp[2]*( dxyz_b.x*hixv[ii]*fv[j] + dxyz_b.x*gyv[ii]*gxv[j] );
            grdy_a += temp[2]*( dxyz_b.y*hpyv[ii]*fv[j] + dxyz_b.y*gyv[ii]*gyv[j] );
            grdy_a += temp[2]*( dxyz_b.z*hizv[ii]*fv[j] + dxyz_b.z*gyv[ii]*gzv[j] );

            grdz_a += temp[0]*gzv[ii]*fv[j];
            grdz_a += temp[1]*( dxyz_a.x*hiyv[ii]*fv[j] + dxyz_a.x*gzv[ii]*gxv[j] );
            grdz_a += temp[1]*( dxyz_a.y*hizv[ii]*fv[j] + dxyz_a.y*gzv[ii]*gyv[j] );
            grdz_a += temp[1]*( dxyz_a.z*hpzv[ii]*fv[j] + dxyz_a.z*gzv[ii]*gzv[j] );
            grdz_a += temp[2]*( dxyz_b.x*hiyv[ii]*fv[j] + dxyz_b.x*gzv[ii]*gxv[j] );
            grdz_a += temp[2]*( dxyz_b.y*hizv[ii]*fv[j] + dxyz_b.y*gzv[ii]*gyv[j] );
            grdz_a += temp[2]*( dxyz_b.z*hpzv[ii]*fv[j] + dxyz_b.z*gzv[ii]*gzv[j] );

            // beta
            temp[0] = 1.0f * coef_b[0] * Pb;
            temp[1] = 2.0f * coef_b[1] * Pb;
            temp[2] = 1.0f * coef_b[2] * Pb;

            grdx_b += temp[0]*gxv[ii]*fv[j];
            grdx_b += temp[1]*( dxyz_b.x*hpxv[ii]*fv[j] + dxyz_b.x*gxv[ii]*gxv[j] );
            grdx_b += temp[1]*( dxyz_b.y*hixv[ii]*fv[j] + dxyz_b.y*gxv[ii]*gyv[j] );
            grdx_b += temp[1]*( dxyz_b.z*hiyv[ii]*fv[j] + dxyz_b.z*gxv[ii]*gzv[j] );
            grdx_b += temp[2]*( dxyz_a.x*hpxv[ii]*fv[j] + dxyz_a.x*gxv[ii]*gxv[j] );
            grdx_b += temp[2]*( dxyz_a.y*hixv[ii]*fv[j] + dxyz_a.y*gxv[ii]*gyv[j] );
            grdx_b += temp[2]*( dxyz_a.z*hiyv[ii]*fv[j] + dxyz_a.z*gxv[ii]*gzv[j] );

            grdy_b += temp[0]*gyv[ii]*fv[j];
            grdy_b += temp[1]*( dxyz_b.x*hixv[ii]*fv[j] + dxyz_b.x*gyv[ii]*gxv[j] );
            grdy_b += temp[1]*( dxyz_b.y*hpyv[ii]*fv[j] + dxyz_b.y*gyv[ii]*gyv[j] );
            grdy_b += temp[1]*( dxyz_b.z*hizv[ii]*fv[j] + dxyz_b.z*gyv[ii]*gzv[j] );
            grdy_b += temp[2]*( dxyz_a.x*hixv[ii]*fv[j] + dxyz_a.x*gyv[ii]*gxv[j] );
            grdy_b += temp[2]*( dxyz_a.y*hpyv[ii]*fv[j] + dxyz_a.y*gyv[ii]*gyv[j] );
            grdy_b += temp[2]*( dxyz_a.z*hizv[ii]*fv[j] + dxyz_a.z*gyv[ii]*gzv[j] );

            grdz_b += temp[0]*gzv[ii]*fv[j];
            grdz_b += temp[1]*( dxyz_b.x*hiyv[ii]*fv[j] + dxyz_b.x*gzv[ii]*gxv[j] );
            grdz_b += temp[1]*( dxyz_b.y*hizv[ii]*fv[j] + dxyz_b.y*gzv[ii]*gyv[j] );
            grdz_b += temp[1]*( dxyz_b.z*hpzv[ii]*fv[j] + dxyz_b.z*gzv[ii]*gzv[j] );
            grdz_b += temp[2]*( dxyz_a.x*hiyv[ii]*fv[j] + dxyz_a.x*gzv[ii]*gxv[j] );
            grdz_b += temp[2]*( dxyz_a.y*hizv[ii]*fv[j] + dxyz_a.y*gzv[ii]*gyv[j] );
            grdz_b += temp[2]*( dxyz_a.z*hpzv[ii]*fv[j] + dxyz_a.z*gzv[ii]*gzv[j] );
         } // end j
      } // end k
      ddx_a(nuc) += grdx_a; ddy_a(nuc) += grdy_a; ddz_a(nuc) += grdz_a;
      ddx_b(nuc) += grdx_b; ddy_b(nuc) += grdy_b; ddz_b(nuc) += grdz_b;
   }
   for (int i = 0; i < local_atoms; i++) {
      // alpha
      smallFor_a[i*3]   -= ddx_a(i) * wp;
      smallFor_a[i*3+1] -= ddy_a(i) * wp;
      smallFor_a[i*3+2] -= ddz_a(i) * wp;

      // beta
      smallFor_b[i*3]   -= ddx_b(i) * wp;
      smallFor_b[i*3+1] -= ddy_b(i) * wp;
      smallFor_b[i*3+2] -= ddz_b(i) * wp;
   }
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}
