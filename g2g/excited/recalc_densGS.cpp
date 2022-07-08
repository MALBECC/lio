#include <iostream>
#include <omp.h>

#include <stdio.h>
#include <string.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"
#include "../libxc/libxcproxy.h"

using namespace G2G;
extern Partition partition;

namespace G2G {

template <class scalar_type>
void PointGroupCPU<scalar_type>::recalc_densGS(const scalar_type* fv, // INPUTS
                                               const scalar_type* gxv, 
                                               const scalar_type* gyv, 
                                               const scalar_type* gzv,
                                               HostMatrix<scalar_type>& gsDens,
                                               HostMatrix<scalar_type>& trDens,
                                               int point,
                                               HostMatrix<scalar_type>& gsOut, // OUTPUTS
                                               HostMatrix<scalar_type>& trOut)
{
   int save = fortran_vars.den_point_save;
   const uint group_m = this->total_functions();
   gsOut(0) = gsOut(1) = gsOut(2) = gsOut(3) = 0.0f;
   trOut(0) = trOut(1) = trOut(2) = trOut(3) = 0.0f;

   // We saved the GS density and derivatives
   if ( save == 0 ) {

      // Transition Density
      for(int i=0;i<group_m;i++) {
         double z3xc, z3yc, z3zc, z; z3xc = z3yc = z3zc = z = 0.0f;
         for(int j=0;j<=i;j++) {
            z += fv[j] * trDens(i,j);
            z3xc += gxv[j] * trDens(i,j);
            z3yc += gyv[j] * trDens(i,j);
            z3zc += gzv[j] * trDens(i,j);
         }
         const double Fi = fv[i];
         const double gx = gxv[i], gy = gyv[i], gz = gzv[i];
         trOut(0) += Fi * z;
         trOut(1) += gx * z + z3xc * Fi;
         trOut(2) += gy * z + z3yc * Fi;
         trOut(3) += gz * z + z3zc * Fi;
      }
  
      // Ground State Density
      gsOut(0) = rho_values(0,point);
      gsOut(1) = rho_values(1,point);
      gsOut(2) = rho_values(2,point);
      gsOut(3) = rho_values(3,point);

   // We have to calculated the GS density and derivatives
   } else {

      // Ground and Transition Densities
      for(int i=0;i<group_m;i++) {
         double z3xc, z3yc, z3zc, z; z3xc = z3yc = z3zc = z = 0.0f;
         double w3xc, w3yc, w3zc, w; w3xc = w3yc = w3zc = w = 0.0f;
         const scalar_type* rm = gsDens.row(i); 
         for(int j=0;j<=i;j++) {
            // Ground state
            const double rmj = rm[j];
            w += fv[j] * rmj;
            w3xc += gxv[j] * rmj;
            w3yc += gyv[j] * rmj;
            w3zc += gzv[j] * rmj;

            // Transition
            z += fv[j] * trDens(i,j);
            z3xc += gxv[j] * trDens(i,j);
            z3yc += gyv[j] * trDens(i,j);
            z3zc += gzv[j] * trDens(i,j);
         }
         const double Fi = fv[i];
         const double gx = gxv[i], gy = gyv[i], gz = gzv[i];
         // Ground state
         gsOut(0) += Fi * w;
         gsOut(1) += gx * w + w3xc * Fi;
         gsOut(2) += gy * w + w3yc * Fi;
         gsOut(3) += gz * w + w3zc * Fi;

         // Transition
         trOut(0) += Fi * z;
         trOut(1) += gx * z + z3xc * Fi;
         trOut(2) += gy * z + z3yc * Fi;
         trOut(3) += gz * z + z3zc * Fi;
      }
   }
}

template <class scalar_type>
void PointGroupCPU<scalar_type>::recalc_densGS3(const scalar_type* fv, // INPUTS
                                               const scalar_type* gxv,
                                               const scalar_type* gyv,
                                               const scalar_type* gzv,
                                               HostMatrix<scalar_type>& gsDens,
                                               HostMatrix<scalar_type>& pdDens,
                                               HostMatrix<scalar_type>& trDens,
                                               int point,
                                               HostMatrix<scalar_type>& gsOut, // OUTPUTS
                                               HostMatrix<scalar_type>& pdOut,
                                               HostMatrix<scalar_type>& trOut)
{
   int save = fortran_vars.den_point_save;
   const uint group_m = this->total_functions();
   gsOut(0) = gsOut(1) = gsOut(2) = gsOut(3) = 0.0f;
   pdOut(0) = pdOut(1) = pdOut(2) = pdOut(3) = 0.0f;
   trOut(0) = trOut(1) = trOut(2) = trOut(3) = 0.0f;

   // We saved the GS density and derivatives
   if ( save == 0 ) {

      // Transition Density
      for(int i=0;i<group_m;i++) {
         double q3xc, q3yc, q3zc, q; q3xc = q3yc = q3zc = q = 0.0f;
         double z3xc, z3yc, z3zc, z; z3xc = z3yc = z3zc = z = 0.0f;
         for(int j=0;j<=i;j++) {
            // Difference Relaxed Excited State Density
            q += fv[j] * pdDens(i,j);
            q3xc += gxv[j] * pdDens(i,j);
            q3yc += gyv[j] * pdDens(i,j);
            q3zc += gzv[j] * pdDens(i,j);

            // Transition Density
            z += fv[j] * trDens(i,j);
            z3xc += gxv[j] * trDens(i,j);
            z3yc += gyv[j] * trDens(i,j);
            z3zc += gzv[j] * trDens(i,j);
         }
         const double Fi = fv[i];
         const double gx = gxv[i], gy = gyv[i], gz = gzv[i];
         // Difference Relaxed Excited State Density
         pdOut(0) += Fi * q;
         pdOut(1) += gx * q + q3xc * Fi;
         pdOut(2) += gy * q + q3yc * Fi;
         pdOut(3) += gz * q + q3zc * Fi;

         // Transition Density
         trOut(0) += Fi * z;
         trOut(1) += gx * z + z3xc * Fi;
         trOut(2) += gy * z + z3yc * Fi;
         trOut(3) += gz * z + z3zc * Fi;
      }

      // Ground State Density
      gsOut(0) = rho_values(0,point);
      gsOut(1) = rho_values(1,point);
      gsOut(2) = rho_values(2,point);
      gsOut(3) = rho_values(3,point);

   // We have to calculated the GS density and derivatives
   } else {

      // Ground and Transition Densities
      for(int i=0;i<group_m;i++) {
         double w3xc, w3yc, w3zc, w; w3xc = w3yc = w3zc = w = 0.0f;
         double q3xc, q3yc, q3zc, q; q3xc = q3yc = q3zc = q = 0.0f;
         double z3xc, z3yc, z3zc, z; z3xc = z3yc = z3zc = z = 0.0f;
         const scalar_type* rm = gsDens.row(i);
         for(int j=0;j<=i;j++) {
            // Ground state Density
            const double rmj = rm[j];
            w += fv[j] * rmj;
            w3xc += gxv[j] * rmj;
            w3yc += gyv[j] * rmj;
            w3zc += gzv[j] * rmj;

            // Difference Relaxed Excited State Density
            q += fv[j] * pdDens(i,j);
            q3xc += gxv[j] * pdDens(i,j);
            q3yc += gyv[j] * pdDens(i,j);
            q3zc += gzv[j] * pdDens(i,j);

            // Transition Density
            z += fv[j] * trDens(i,j);
            z3xc += gxv[j] * trDens(i,j);
            z3yc += gyv[j] * trDens(i,j);
            z3zc += gzv[j] * trDens(i,j);
         }
         const double Fi = fv[i];
         const double gx = gxv[i], gy = gyv[i], gz = gzv[i];
         // Ground state
         gsOut(0) += Fi * w;
         gsOut(1) += gx * w + w3xc * Fi;
         gsOut(2) += gy * w + w3yc * Fi;
         gsOut(3) += gz * w + w3zc * Fi;

         // Difference Relaxed Excited State Density
         pdOut(0) += Fi * q;
         pdOut(1) += gx * q + q3xc * Fi;
         pdOut(2) += gy * q + q3yc * Fi;
         pdOut(3) += gz * q + q3zc * Fi;

         // Transition Density
         trOut(0) += Fi * z;
         trOut(1) += gx * z + z3xc * Fi;
         trOut(2) += gy * z + z3yc * Fi;
         trOut(3) += gz * z + z3zc * Fi;
      }
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
