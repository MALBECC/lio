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

extern "C" void g2g_calculateg_(double* Tmat,double* F,int& DER)
{
   partition.solve_Gxc(Tmat,F,DER);
}

namespace G2G {

void Partition::solve_Gxc(double* Tmat,double* F,int DER)
{
   int M = fortran_vars.m;
   std::vector< HostMatrix<double> > Fock_output(G2G::cpu_threads + G2G::gpu_threads);

#pragma omp parallel for num_threads(cpu_threads+gpu_threads) schedule(static)
   for(uint i=0;i<work.size();i++) {
#if GPU_KERNELS
      bool gpu_thread = false;
      if (i >= cpu_threads) {
         gpu_thread = true;
         cudaSetDevice(i - cpu_threads);
      }
#endif
      Fock_output[i].resize(fortran_vars.rmm_output.width,fortran_vars.rmm_output.height);
      Fock_output[i].zero();
      for(uint j=0;j<work[i].size();j++) {
         Timer element;
         element.start_and_sync();
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind-cubes.size()]->solve_3rd_der(Tmat,Fock_output[i],DER);
         } else {
           cubes[ind]->solve_3rd_der(Tmat,Fock_output[i],DER);
         }
#if GPU_KERNELS
        if (gpu_thread) cudaDeviceSynchronize();
#endif
        element.stop_and_sync();
      }
   }
   // Join results of threads
   for(uint i=0; i<work.size(); i++) {
     int index = 0;
     for(uint j=0; j<M; j++) {
       F[j*M+j] += Fock_output[i](index,0);
       index += 1;
       for(uint k=j+1; k<M; k++) {
           F[j*M+k] += Fock_output[i](index,0);
           F[k*M+j] = F[j*M+k];
           index += 1;
       }
     }
   }
   std::vector<HostMatrix<double> >().swap(Fock_output);
} // END solve_Gxc

template<class scalar_type> void PointGroupCPU<scalar_type>::
               solve_3rd_der(double* Tmat,HostMatrix<double>& Fock,int DER)
{
   const uint group_m = this->total_functions();
   const int npoints = this->points.size();
   bool lda = false;
   bool compute_forces = true;
   compute_functions(compute_forces,!lda);
   HostMatrix<scalar_type> rmm_input(group_m,group_m);

   int M = fortran_vars.m;
   get_rmm_input(rmm_input);

   double* zcoef = (double*) malloc(3*sizeof(double));
   double* precond = (double*) malloc(group_m*sizeof(double));
   double* smallFock = (double*) malloc(group_m*group_m*sizeof(double));
   memset(smallFock,0.0f,group_m*group_m*sizeof(double));

// Obtain reduced matrix for this group
   HostMatrix<double> tred(group_m,group_m);
   HostMatrix<double> Tbig(M*(M+1)/2);
   int index = 0;
   for(int row=0; row<M; row++) {
     Tbig(index) = Tmat[row*M+row];
     index += 1;
     for(int col=row+1; col<M; col++) {
       Tbig(index) = Tmat[row*M+col] + Tmat[col*M+row];
       index += 1;
     }
   }
   get_tred_input(tred,Tbig); Tbig.deallocate();

// INITIALIZATION LIBXC
#define libxc_init_param \
  fortran_vars.func_id, fortran_vars.func_coef, fortran_vars.nx_func, \
  fortran_vars.nc_func, fortran_vars.nsr_id, fortran_vars.screen, \
  XC_UNPOLARIZED
  LibxcProxy<scalar_type,3> libxcProxy(libxc_init_param);
#undef libxc_init_param

   HostMatrix<scalar_type> groundD(4);
   HostMatrix<scalar_type> transD(4);

   for(int point=0;point<npoints;point++) {
      scalar_type pd, pdx, pdy, pdz; pd = pdx = pdy = pdz = 0.0f;
      scalar_type red, redx, redy, redz; red = redx = redy = redz = 0.0f;
      const scalar_type* fv = function_values.row(point);
      const scalar_type* gxv = gX.row(point);
      const scalar_type* gyv = gY.row(point);
      const scalar_type* gzv = gZ.row(point);

      // Calculate GS and transition densities and derivatives in the point
      #define recalc_params \
      function_values.row(point), gX.row(point), gY.row(point), gZ.row(point), \
      rmm_input, tred, point, \
      groundD, transD
      recalc_densGS(recalc_params);
      #undef recalc_params

      // Copying outputs
      pd  = groundD(0); pdx = groundD(1); pdy = groundD(2); pdz = groundD(3);
      red = transD(0); redx = transD(1); redy = transD(2); redz = transD(3);

      double sigma = pdx * pdx + pdy * pdy + pdz * pdz;

      // obtain derivatives terms
      if (DER == 2 ) {
          double cruz = (redx * pdx + redy * pdy + redz * pdz) * 0.5f;
          libxcProxy.coefLR(&pd,&sigma,red,cruz,zcoef);
      } else {
          libxcProxy.coefZv(pd,sigma,pdx,pdy,pdz,red,redx,redy,
                            redz, zcoef);
      }
      pdx *= 0.5f; pdy *= 0.5f; pdz *= 0.5f;

      const scalar_type wp = this->points[point].weight;
      double term1, term2, term3, term4, precondii, result;
      for(int i=0; i<group_m; i++) {
        term1 = zcoef[0] * 0.5f * fv[i] + zcoef[1] * pdx * gxv[i];
        term2 = zcoef[1] * pdy * gyv[i] + zcoef[1] * pdz * gzv[i];
        term3 = zcoef[2] * redx * gxv[i] + zcoef[2] * redy * gyv[i];
        term4 = zcoef[2] * redz * gzv[i];
        precond[i] = (term1 + term2 + term3 + term4) * wp;
        precondii = precond[i];
        for(int j=0; j<=i; j++) {
           result = fv[i] * precond[j] + precondii * fv[j];
           smallFock[i*group_m+j] += result;
        }
      }

   }  // END points loop

   const int indexes = this->rmm_bigs.size();
   for (int i = 0; i < indexes; i++) {
      int bi = this->rmm_bigs[i], row = this->rmm_rows[i],
          col = this->rmm_cols[i];
          Fock(bi) += smallFock[col*group_m+row];
   }

   // Free Memory
   free(smallFock); smallFock = NULL;
   free(precond); precond = NULL;
   free(zcoef); zcoef = NULL;
   rmm_input.deallocate(); groundD.deallocate(); 
   tred.deallocate(); transD.deallocate();
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}
