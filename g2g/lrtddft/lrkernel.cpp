#include <iostream>
#include <xc.h>
#include <omp.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"
#include "../matrix.h"

#include "lrkernel.h"
#include "eri.h"
#include "couplingform.h"

using namespace G2G;
extern Partition partition;

//######################################################################
//######################################################################
#if FULL_DOUBLE
extern "C" g2g_linear_response_(double* MatCoef,double* Kfxc,double* KcMat,
                                double* Kxc_int,double* Kc_int,double* Cbas,
                                int& dim)
#else
extern "C" g2g_linear_response_(float* MatCoef,float* Kfxc,float* KcMat,
                                float* Kxc_int,float* Kc_int,float* Cbas,
                                int& dim)
#endif
{
   int M = fortran_vars.m;
   int nco = fortran_vars.nco;
   fortran_vars.dim = dim;
   int s_func = fortran_vars.s_funcs;
   int p_func = fortran_vars.p_funcs;
   int d_func = fortran_vars.d_funcs;
#if FULL_DOUBLE
   fortran_vars.MatCoef = FortranMatrix<double>(MatCoef,M,M,M);
   FortranMatrix<double> Kxc=FortranMatrix<double>(Kfxc,dim,dim,dim,dim);
   FortranMatrix<double> Kc=FortranMatrix<double>(KcMat,dim,dim,dim,dim);
   double* aContr = &fortran_vars.a_values(0,0);
   double* pos = &fortran_vars.atom_positions_pointer(0,0);
#else
   fortran_vars.MatCoef = FortranMatrix<float>(MatCoef,M,M,M);
   FortranMatrix<float> Kxc=FortranMatrix<float>(Kfxc,dim,dim,dim,dim);
   FortranMatrix<float> Kc=FortranMatrix<float>(KcMat,dim,dim,dim,dim);
   float* aContr = &fortran_vars.a_values(0,0);
   float* pos = &fortran_vars.atom_positions_pointer(0,0);
#endif
   double timeI, timeF;

   timeI = omp_get_wtime();
   eri(Kc_int,M,nco,dim,fortran_vars.atoms,ncont,Cbas,aContr,pos,nuc,
       s_func,p_func,d_func);
   timeF = omp_get_wtime();
   printf("ERI SUBROUTINE %f\n",timeF-timeI);
   timeI = timeF = 0.0;

   timeI = omp_get_wtime();
   partition.solve_lr(Kxc,Kc,Kxc_int,Kc_int);
   timeF = omp_get_wtime();
   printf("SOLVE_LR SUBROUTINE %f\n",timeF-timeI);

   fflush(stdout); // NOT BUFFERED
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################
namespace G2G {
#if FULL_DOUBLE
void Partition::solve_lr(FortranMatrix<double>& Kxc,FortranMatrix<double>& Kc,
                         double* calcKfxc,double* calcKc)
#else
void Partition::solve_lr(FortranMatrix<float>& Kxc,FortranMatrix<float>& Kc,
                         float* calcKfxc,float* calcKc)
#endif
{
   xc_func_type funcx, funcc;
   switch (fortran_vars.iexch) {
      case 9: {
             xc_func_init(&funcx,XC_GGA_X_PBE,XC_UNPOLARIZED);
             xc_func_init(&funcc,XC_GGA_X_PBE,XC_UNPOLARIZED);
      } break;
      default: {
             std::cout << "LINEAR RESPONSE ONLY WORKS WITH PBE" << std::endl;
      } break;
   }
   double timeI, timeF;
   timeI = omp_get_wtime();
#pragma omp parallel for schedule(static)
   for(uint i=0;i<work.size();i++) {
      for(uint j=0;j<work[i].size();j++) {
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind-cubes.size()]->solve_closed_lr(funcx,funcc,calcKfxc);
         } else {
           cubes[ind]->solve_closed_lr(funcx,funcc,calcKfxc);
         }
      }
   }
   timeF = omp_get_wtime();
   printf("SOLVE_CLOSED_LR SUBROUTINE %f\n",timeF-timeI);
   xc_func_end(&funcx);
   xc_func_end(&funcc);

   CouplingForm(Kxc,Kc,calcKfxc,calcKc);
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################
#if FULL_DOUBLE
template<class scalar_type> void PointGroupCPU<scalar_type>::
       solve_closed_lr(xc_func_type funcx,xc_func_type funcc, double* gIntKfxc)
#else
template<class scalar_type> void PointGroupCPU<scalar_type>::
       solve_closed_lr(xc_func_type funcx,xc_func_type funcc, float* gIntKfxc)
#endif
{
   const uint group_m = this->total_functions();
   const int npoints = this->points.size();
   const int iexch = fortran_vars.iexch;
   bool lda = false;
   bool compute_forces = true;
   compute_functions(compute_forces,!lda);
   HostMatrix<scalar_type> rmm_input(group_m,group_m);
   int* numeros = new int[group_m];
   int NCO = fortran_vars.nco;
   int M = fortran_vars.m;
   int Nvirt = M - NCO;
   int M3, M2;
   M2 = M*M;
   M3 = M2*M;
   get_rmm_input(rmm_input);
   get_coef_input(rmm_input,numeros);
#if FULL_DOUBLE
   double fxc, fx, fc, dens, trans1, trans2, scratch;
#else
   float fxc, fx, fc, dens, trans1, trans2, scratch;
#endif
   fxc = fx = fc = dens = trans1 = trans2 = scratch = 0.0;
   
   for(int point=0;point<npoints;point++) {
      scalar_type pd, tdx, tdy, tdz; pd = tdx = tdy = tdz = 0.0;
      const scalar_type* fv = function_values.row(point);
      const scalar_type* gxv = gX.row(point);
      const scalar_type* gyv = gY.row(point);
      const scalar_type* gzv = gZ.row(point);
      for(int i=0;i<group_m;i++) {
         scalar_type w3xc, w3yc, w3zc, w; w3xc = w3yc = w3zc = w = 0;
         const scalar_type* rm = rmm_input.row(i);
         for(int j=0;j<=i;j++) {
            const scalar_type rmj = rm[j];
            w += fv[j] * rmj;
            w3xc += gxv[j] * rmj;
            w3yc += gyv[j] * rmj;
            w3zc += gzv[j] * rmj;
         }
         const scalar_type Fi = fv[i];
         const scalar_type gx = gxv[i], gy = gyv[i], gz = gzv[i];
         pd += Fi * w;
         tdx += gx * w + w3xc * Fi;
         tdy += gy * w + w3yc * Fi;
         tdz += gz * w + w3zc * Fi;
      }
      if(fabs(pd) < 1e-5) {
        const vec_type3 dxyz(tdx, tdy, tdz);
#if FULL_DOUBLE
        double grad = dxyz.x * dxyz.x + dxyz.y * dxyz.y + dxyz.z * dxyz.z;
        double sigma = sqrt(grad);
        double v2rhosigma,v2sigma2;
#else
        float grad = dxyz.x * dxyz.x + dxyz.y * dxyz.y + dxyz.z * dxyz.z;
        float sigma = sqrt(grad);
        float v2rhosigma,v2sigma2;
#endif
        if(grad < MIN_PRECISION) grad = (scalar_type)MIN_PRECISION;
        dens = pd;
        xc_gga_fxc(&funcx,1,&dens,&sigma,&fx,&v2rhosigma,&v2sigma2);
        xc_gga_fxc(&funcc,1,&dens,&sigma,&fc,&v2rhosigma,&v2sigma2);
        fxc = fx + fc;
        const scalar_type wp = this->points[point].weight;
        int row1, col1, row2, col2;
        for(int ic=0;ic<group_m;ic++) {
           row1 = numeros[ic];
           for(int jc=0;jc<=ic;jc++) {
              col1 = numeros[jc];
              trans1 = fv[ic]*fv[jc]*fxc*wp;
              if(fabs(trans1) > 1e-10) {
                 for(int kc=0;kc<group_m;kc++) {
                    row2 = numeros[kc];
                    for(int lc=0;lc<=kc;lc++) {
                       col2 = numeros[lc];
                       trans2 = trans1*fv[kc]*fv[lc];
                       gIntKfxc[row1*M3+col1*M2+row2*M+col2] += trans2;
                       scratch=gIntKfxc[row1*M3+col1*M2+row2*M+col2];
                       gIntKfxc[row1*M3+col1*M2+col2*M+row2] = scratch;
                       gIntKfxc[col1*M3+row1*M2+row2*M+col2] = scratch;
                       gIntKfxc[col1*M3+row1*M2+col2*M+row2] = scratch;
                    }
                 }
              } // END IF TRANS1
           }
        }
      } // END IF DENSITY
   }  // END points loop

   delete numeros;
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################
template<class scalar_type>
#if FULL_DOUBLE
void PointGroupCPU<scalar_type>::get_coef_input(HostMatrix<scalar_type>&
               rmm_input,int* vecnum, FortranMatrix<double>& Coef) const
#else
void PointGroupCPU<scalar_type>::get_coef_input(HostMatrix<scalar_type>&
               rmm_input,int* vecnum, FortranMatrix<float>& Coef) const
#endif
{
   const int indexes = this->rmm_bigs.size();
   std::vector<int> row;
   for(int i = 0; i < indexes; i++) {
      int bi = this->rmm_bigs[i];
      for(int l=1; l<=fortran_vars.m; l++) {
         for(int k=1; k<=(l-1); k++) {
            int idx = l + (2*fortran_vars.m-k)*(k-1)/2;
            if(bi+1 == idx) {
              row.push_back(l-1);
              row.push_back(k-1);
            }
         }
         int idx = l+(2*fortran_vars.m-l)*(l-1)/2;
         if(bi+1 == idx) {
           row.push_back(l-1);
         }
      }
   }
// delete of repeated indexes
   int numeros[row.size()];
   int r = 0, i, k;
   for(i=0,k=0;i<row.size(),k<row.size();i++,k++) {
      numeros[i]=row[r];r++;
      for(int j =i-1; j>=0;j--) {
        if(numeros[i] == numeros[j])
        {
           i--;break;
        }
     }
  }
  for(int wx=0; wx<i; wx++)
        vecnum[wx] = numeros[wx];
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################
template<class scalar_type>
void PointGroupCPU<scalar_type>::get_coef_input
          (HostMatrix<scalar_type>& rmm_input,int* nume) const
{
  get_coef_input(rmm_input, nume, fortran_vars.MatCoef);
}
//######################################################################
//######################################################################

template class PointGroup<double>;
template class PointGroup<float>;
template class PointGroupCPU<double>;
template class PointGroupCPU<float>;
}
