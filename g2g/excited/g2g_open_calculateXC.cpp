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

extern "C" void g2g_open_calculatexc_(double* TmatA, double* TmatB,
                                      double* FvA, double* FvB, int& DER)
{
   partition.solve_lr(TmatA,TmatB,FvA,FvB,DER);
}

namespace G2G {

void Partition::solve_lr(double* Ta,double* Tb,double* Fa,double* Fb,int DER)
{
   int M = fortran_vars.m;
   std::vector< HostMatrix<double> > Fock_output_a(G2G::cpu_threads + G2G::gpu_threads);
   std::vector< HostMatrix<double> > Fock_output_b(G2G::cpu_threads + G2G::gpu_threads);

#pragma omp parallel for num_threads(cpu_threads+gpu_threads) schedule(static)
   for(uint i=0;i<work.size();i++) {
#if GPU_KERNELS
    bool gpu_thread = false;
    if (i >= cpu_threads) {
      gpu_thread = true;
      cudaSetDevice(i - cpu_threads);
    }
#endif
     Fock_output_a[i].resize(fortran_vars.rmm_output_a.width,fortran_vars.rmm_output_a.height);
     Fock_output_a[i].zero();
     Fock_output_b[i].resize(fortran_vars.rmm_output_b.width,fortran_vars.rmm_output_b.height);
     Fock_output_b[i].zero();

     for(uint j=0;j<work[i].size();j++) {
        Timer element;
        element.start_and_sync();
        int ind = work[i][j];
        if(ind >= cubes.size()) {
          spheres[ind - cubes.size()]->solve_open_lr(Ta,Tb,Fock_output_a[i],Fock_output_b[i],DER);
        } else {
          cubes[ind]->solve_open_lr(Ta,Tb,Fock_output_a[i],Fock_output_b[i],DER);
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
       Fa[j*M+j] += Fock_output_a[i](index,0);
       Fb[j*M+j] += Fock_output_b[i](index,0);
       index += 1;
       for(uint k=j+1; k<M; k++) {
           Fa[j*M+k] += Fock_output_a[i](index,0);
           Fa[k*M+j]  = Fa[j*M+k];
           Fb[j*M+k] += Fock_output_b[i](index,0);
           Fb[k*M+j]  = Fb[j*M+k];
           index += 1;
       }
     }
   }
   std::vector<HostMatrix<double>>().swap(Fock_output_a);
   std::vector<HostMatrix<double>>().swap(Fock_output_b);
} // END solve LR

template<class scalar_type> void PointGroupCPU<scalar_type>::solve_open_lr(
                             double* Ta,double* Tb,HostMatrix<double>& FockA,HostMatrix<double>& FockB,int DER)
{
   const uint group_m = this->total_functions();
   const int npoints = this->points.size();
   bool lda = false;
   bool compute_forces = true;
   compute_functions(compute_forces,lda); // get function and gradients
   HostMatrix<scalar_type> rmm_input_a(group_m,group_m), rmm_input_b(group_m, group_m);

   int M = fortran_vars.m;
   get_rmm_input(rmm_input_a,rmm_input_b);

   double* lrCoef_a = (double*) malloc(8*sizeof(double));
   double* lrCoef_b = (double*) malloc(8*sizeof(double));
   double* smallFock_a = (double*)malloc(group_m*group_m*sizeof(double));
   double* smallFock_b = (double*)malloc(group_m*group_m*sizeof(double));
   memset(smallFock_a,0.0f,group_m*group_m*sizeof(double));
   memset(smallFock_b,0.0f,group_m*group_m*sizeof(double));

// Obtain the reduced matrix of excited states for this group
   HostMatrix<double> tred_a(group_m,group_m), tred_b(group_m,group_m);
   HostMatrix<double> Tbig_a(M*(M+1)/2), Tbig_b(M*(M+1)/2);
   int index = 0;
   for(int row=0; row<M; row++) {
     Tbig_a(index) = Ta[row*M+row];
     Tbig_b(index) = Tb[row*M+row];
     index += 1;
     for(int col=row+1; col<M; col++) {
        Tbig_a(index) = Ta[row*M+col] + Ta[col*M+row];
        Tbig_b(index) = Tb[row*M+col] + Tb[col*M+row];
        index += 1;
     }
   }
   get_tred_input(tred_a,Tbig_a); Tbig_a.deallocate();
   get_tred_input(tred_b,Tbig_b); Tbig_b.deallocate();

//LIBXC INITIALIZATION
#define libxc_init_param \
  fortran_vars.func_id, fortran_vars.func_coef, fortran_vars.nx_func, \
  fortran_vars.nc_func, fortran_vars.nsr_id, fortran_vars.screen, \
  XC_POLARIZED
  LibxcProxy<scalar_type,3> libxcProxy(libxc_init_param);
#undef libxc_init_param

   HostMatrix<scalar_type> groundD_a(4), groundD_b(4);
   HostMatrix<scalar_type> transD_a(4), transD_b(4);
   double* pre1a = (double*)malloc(group_m*sizeof(double));
   double* pre2a = (double*)malloc(group_m*sizeof(double));
   double* pre1b = (double*)malloc(group_m*sizeof(double));
   double* pre2b = (double*)malloc(group_m*sizeof(double));

   for(int point=0;point<npoints;point++) {
      const scalar_type* fv = function_values.row(point);
      const scalar_type* gxv = gX.row(point);
      const scalar_type* gyv = gY.row(point);
      const scalar_type* gzv = gZ.row(point);
      groundD_a(0) = groundD_a(1) = groundD_a(2) = groundD_a(3) = 0.0f;
      groundD_b(0) = groundD_b(1) = groundD_b(2) = groundD_b(3) = 0.0f;
      transD_a(0) = transD_a(1) = transD_a(2) = transD_a(3) = 0.0f;
      transD_b(0) = transD_b(1) = transD_b(2) = transD_b(3) = 0.0f;

      // Ground and Transition Densities
      for(int i=0;i<group_m;i++) {
         double z3xc_a, z3yc_a, z3zc_a, z_a; z3xc_a = z3yc_a = z3zc_a = z_a = 0.0f;
         double z3xc_b, z3yc_b, z3zc_b, z_b; z3xc_b = z3yc_b = z3zc_b = z_b = 0.0f;
         double w3xc_a, w3yc_a, w3zc_a, w_a; w3xc_a = w3yc_a = w3zc_a = w_a = 0.0f;
         double w3xc_b, w3yc_b, w3zc_b, w_b; w3xc_b = w3yc_b = w3zc_b = w_b = 0.0f;

         const scalar_type* rma = rmm_input_a.row(i);
         const scalar_type* rmb = rmm_input_b.row(i);
         const scalar_type* tra = tred_a.row(i);
         const scalar_type* trb = tred_b.row(i);

         for(int j=0;j<=i;j++) {
            // Ground state alpha
            const double rmaj = rma[j];
            w_a += fv[j] * rmaj;
            w3xc_a += gxv[j] * rmaj;
            w3yc_a += gyv[j] * rmaj;
            w3zc_a += gzv[j] * rmaj;
            // Ground state beta
            const double rmbj = rmb[j];
            w_b += fv[j] * rmbj;
            w3xc_b += gxv[j] * rmbj;
            w3yc_b += gyv[j] * rmbj;
            w3zc_b += gzv[j] * rmbj;

            // Transition alpha
            const double traj = tred_a(i,j);
            z_a += fv[j] * traj;
            z3xc_a += gxv[j] * traj;
            z3yc_a += gyv[j] * traj;
            z3zc_a += gzv[j] * traj;
            // Transition beta
            const double trbj = tred_b(i,j);
            z_b += fv[j] * trbj;
            z3xc_b += gxv[j] * trbj;
            z3yc_b += gyv[j] * trbj;
            z3zc_b += gzv[j] * trbj;
         }
         const double Fi = fv[i];
         const double gx = gxv[i], gy = gyv[i], gz = gzv[i];
         // Ground state alpha
         groundD_a(0) += Fi * w_a;
         groundD_a(1) += gx * w_a + w3xc_a * Fi;
         groundD_a(2) += gy * w_a + w3yc_a * Fi;
         groundD_a(3) += gz * w_a + w3zc_a * Fi;
         // Ground state beta
         groundD_b(0) += Fi * w_b;
         groundD_b(1) += gx * w_b + w3xc_b * Fi;
         groundD_b(2) += gy * w_b + w3yc_b * Fi;
         groundD_b(3) += gz * w_b + w3zc_b * Fi;

         // Transition alpha
         transD_a(0) += Fi * z_a;
         transD_a(1) += gx * z_a + z3xc_a * Fi;
         transD_a(2) += gy * z_a + z3yc_a * Fi;
         transD_a(3) += gz * z_a + z3zc_a * Fi;
         // Transition beta
         transD_b(0) += Fi * z_b;
         transD_b(1) += gx * z_b + z3xc_b * Fi;
         transD_b(2) += gy * z_b + z3yc_b * Fi;
         transD_b(3) += gz * z_b + z3zc_b * Fi;
      }

      // obtain derivatives terms
      if (DER == 2 ) {
          libxcProxy.coefLR(groundD_a,groundD_b,transD_a,transD_b,lrCoef_a,lrCoef_b);
      } else {
        cout << "Third functional derivative in open shell is not available yet." << endl;
        exit(-1);
      }

      const scalar_type wp = this->points[point].weight;
      double term1, term2, term3, term4, term5, result;
      double fac1, fac2, fac3, fac4;
      double precii_1a, precii_2a, precii_1b, precii_2b;

      for(int ii=0; ii<group_m; ii++) {
         term1 = 0.5f * fv[ii];
         term2 = (gxv[ii]*groundD_a(1) + gyv[ii]*groundD_a(2) + gzv[ii]*groundD_a(3));
         term3 = (gxv[ii]*groundD_b(1) + gyv[ii]*groundD_b(2) + gzv[ii]*groundD_b(3));
         term4 = (gxv[ii]*transD_a(1) + gyv[ii]*transD_a(2) + gzv[ii]*transD_a(3));
         term5 = (gxv[ii]*transD_b(1) + gyv[ii]*transD_b(2) + gzv[ii]*transD_b(3));

         // alpha-alpha
         fac1 = lrCoef_a[0] * term1;
         fac2 = lrCoef_a[1] * term2;
         fac3 = lrCoef_a[2] * term3;
         fac4 = lrCoef_a[3] * term4;
         pre1a[ii] = (fac1+fac2+fac3+fac4) * wp;
         // alpha-beta
         fac1 = lrCoef_b[4] * term1;
         fac2 = lrCoef_b[5] * term2;
         fac3 = lrCoef_b[6] * term3;
         fac4 = lrCoef_b[7] * term5;
         pre2a[ii] = (fac1+fac2+fac3+fac4) * wp;

         precii_1a = pre1a[ii];
         precii_2a = pre2a[ii];

         // beta-beta
         fac1 = lrCoef_b[0] * term1;
         fac2 = lrCoef_b[1] * term3;
         fac3 = lrCoef_b[2] * term2;
         fac4 = lrCoef_b[3] * term5;
         pre1b[ii] = (fac1+fac2+fac3+fac4) * wp;
         // beta-alpha
         fac1 = lrCoef_a[4] * term1;
         fac2 = lrCoef_a[5] * term3;
         fac3 = lrCoef_a[6] * term2;
         fac4 = lrCoef_a[7] * term4;
         pre2b[ii] = (fac1+fac2+fac3+fac4) * wp;

         precii_1b = pre1b[ii];
         precii_2b = pre2b[ii];

         for(int jj=0; jj<=ii; jj++) {
            // alpha
            result = pre1a[jj] * fv[ii] + precii_1a * fv[jj];
            result += pre2a[jj] * fv[ii] + precii_2a * fv[jj];
            smallFock_a[ii*group_m+jj] += result;

            // beta
            result = pre1b[jj] * fv[ii] + precii_1b * fv[jj];
            result += pre2b[jj] * fv[ii] + precii_2b * fv[jj];
            smallFock_b[ii*group_m+jj] += result;
         }
      }
   }  // END points loop

   const int indexes = this->rmm_bigs.size();
   for (int i = 0; i < indexes; i++) {
      int bi = this->rmm_bigs[i], row = this->rmm_rows[i],
          col = this->rmm_cols[i];
          FockA(bi) += smallFock_a[col*group_m+row];
          FockB(bi) += smallFock_b[col*group_m+row];
   }

   // Free Memory
   free(smallFock_a); smallFock_a  = NULL;
   free(smallFock_b); smallFock_b  = NULL;
   free(lrCoef_a); lrCoef_a = NULL;
   free(lrCoef_b); lrCoef_b = NULL;
   free(pre1a); pre1a = NULL; free(pre2a); pre2a = NULL;
   free(pre1b); pre1b = NULL; free(pre2b); pre2b = NULL;
   tred_a.deallocate(); groundD_a.deallocate(); 
   tred_b.deallocate(); groundD_b.deallocate(); 
   transD_a.deallocate(); rmm_input_a.deallocate();
   transD_b.deallocate(); rmm_input_b.deallocate();
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}
