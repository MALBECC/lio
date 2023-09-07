#include <iostream>
#include <omp.h>

#include <stdio.h>
#include <string.h> 

#include "../common.h"
#include "../init.h"
#include "../partition.h"

using namespace G2G;
extern Partition partition;

//######################################################################
//######################################################################
extern "C" void g2g_saverho_(int& saveRho)
{
   fortran_vars.den_point_save = saveRho;
   if ( saveRho == 0 ) {
      cout << " Saving density and derivatives of Ground State" << endl;
      partition.lr_init();
   } else {
      cout << " Recalculating density and derivatives of Ground State" << endl;
   }


   fflush(stdout); // NOT BUFFERED
}

//######################################################################
//######################################################################

namespace G2G {

void Partition::lr_init()
{

#pragma omp parallel for schedule(static)
    for(uint i=0;i<work.size();i++) {
      for(uint j=0;j<work[i].size();j++) {
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind - cubes.size()]->lr_closed_init();
         } else {
           cubes[ind]->lr_closed_init();
         }
      }
   }
   fflush(stdout);
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################

template<class scalar_type> void PointGroupCPU<scalar_type>::
               lr_closed_init()
{
   const uint group_m = this->total_functions();
   const int npoints = this->points.size();
   bool lda = false;
   bool compute_forces = false;
   compute_functions(compute_forces,!lda);
   HostMatrix<scalar_type> rmm_input(group_m,group_m);
   int* numeros = new int[group_m];
   int M = fortran_vars.m;
   get_rmm_input(rmm_input);
   rho_values.resize(4, npoints);

   for(int point=0; point<npoints; point++) {
      scalar_type pd, tdx, tdy, tdz; pd = tdx = tdy = tdz = 0.0f;
      const scalar_type* fv = function_values.row(point);
      const scalar_type* gxv = gX.row(point);
      const scalar_type* gyv = gY.row(point);
      const scalar_type* gzv = gZ.row(point);
      for(int i=0;i<group_m;i++) {
         double w3xc, w3yc, w3zc, w; w3xc = w3yc = w3zc = w = 0.0f;
         const scalar_type* rm = rmm_input.row(i);
         for(int j=0;j<=i;j++) {
            const double rmj = rm[j];
            w += fv[j] * rmj;
            w3xc += gxv[j] * rmj;
            w3yc += gyv[j] * rmj;
            w3zc += gzv[j] * rmj;
         }
         const double Fi = fv[i];
         const double gx = gxv[i], gy = gyv[i], gz = gzv[i];
         pd += Fi * w;
         tdx += gx * w + w3xc * Fi;
         tdy += gy * w + w3yc * Fi;
         tdz += gz * w + w3zc * Fi;
      }
      // Save Ground State Density and Derivatives
      rho_values(0,point) = pd;
      rho_values(1,point) = tdx;
      rho_values(2,point) = tdy;
      rho_values(3,point) = tdz;
   }  // END points loop
}
//######################################################################
//######################################################################
#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}
