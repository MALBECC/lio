#include <iostream>
#include <xc.h>
#include <omp.h>

#include <stdio.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"

#include "eri.h"
#include "couplingform.h"
#include "../libxc/libxcproxy.h"

#define DENSMIN 1e-5

using namespace G2G;
extern Partition partition;

//######################################################################
//######################################################################
extern "C" void g2g_linear_response_(double* MatCoef,double* KMat,
                                     double* K_int,double* Cbas,int& dim,
                                     int& NCOlr, int& Nvirt)
{
   int M = fortran_vars.m;
   int nco = fortran_vars.nco;
   int s_func = fortran_vars.s_funcs;
   int p_func = fortran_vars.p_funcs;
   int d_func = fortran_vars.d_funcs;
   double* aContr = &fortran_vars.a_values(0,0);
   double* pos = &fortran_vars.atom_positions_pointer(0,0);
   double timeI, timeF;

   fortran_vars.dim = dim;
   fortran_vars.nvirt = Nvirt;
   fortran_vars.ncolr = NCOlr;

   uint* ncont = &fortran_vars.contractions(0);
   uint* nuc = &fortran_vars.nucleii(0);

   timeI = omp_get_wtime();
   eri(K_int,M,fortran_vars.atoms,ncont,Cbas,aContr,pos,nuc,
       s_func,p_func,d_func);
   timeF = omp_get_wtime();

   printf("ERI SUBROUTINE %f\n",timeF-timeI);
   timeI = timeF = 0.0;

   timeI = omp_get_wtime();
   partition.solve_lr(KMat,K_int,MatCoef);
   timeF = omp_get_wtime();
   printf("SOLVE_LR SUBROUTINE %f\n",timeF-timeI);

   fflush(stdout); // NOT BUFFERED

}
//######################################################################
//######################################################################

//######################################################################
//######################################################################

namespace G2G {

void Partition::solve_lr(double* K,double* calcK,double* Coef)
{

   int M = fortran_vars.m;
   int NCO = fortran_vars.ncolr;
   int Nvirt = fortran_vars.nvirt;
   double timeI, timeF;
   timeI = omp_get_wtime();
#pragma omp parallel for schedule(static)
   for(uint i=0;i<work.size();i++) {
      for(uint j=0;j<work[i].size();j++) {
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind-cubes.size()]->solve_closed_lr(calcK);
         } else {
           cubes[ind]->solve_closed_lr(calcK);
         }
      }
   }
   timeF = omp_get_wtime();
   printf("SOLVE_CLOSED_LR SUBROUTINE %f\n",timeF-timeI);

   timeI = omp_get_wtime();
   CouplingForm(K,calcK,M,NCO,Nvirt,Coef);
   timeF = omp_get_wtime();
   printf("COUPLINGFORM SUBROUTINE %f\n",timeF-timeI);

}
//######################################################################
//######################################################################

//######################################################################
//######################################################################

template<class scalar_type> void PointGroupCPU<scalar_type>::
               solve_closed_lr(double* gIntK)
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
   double trans1, trans2, scratch;
   trans1 = trans2 = scratch = 0.0;

//LIBXC INITIALIZATION
   const int nspin = XC_UNPOLARIZED;
   const int functionalExchange = fortran_vars.ex_functional_id; //101;
   const int functionalCorrelation = fortran_vars.ec_functional_id; // 130;
   LibxcProxy<scalar_type,3> libxcProxy(functionalExchange, functionalCorrelation, nspin);
   
   for(int point=0;point<npoints;point++) {
      scalar_type pd, fxc, tdx, tdy, tdz; pd = fxc = tdx = tdy = tdz = 0.0;
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
      if(fabs(pd) < DENSMIN) {
        const vec_type3 dxyz(tdx, tdy, tdz);
        double grad = dxyz.x * dxyz.x + dxyz.y * dxyz.y + dxyz.z * dxyz.z;
        scalar_type sigma = sqrt(grad);
        scalar_type v2rhosigma,v2sigma2; v2rhosigma = v2sigma2 = 0.0;

        if(grad < MIN_PRECISION) grad = (scalar_type)MIN_PRECISION;
        libxcProxy.doGGA(pd,sigma,&fxc,v2rhosigma,v2sigma2);
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
                       gIntK[row1*M3+col1*M2+row2*M+col2] += trans2;
                       scratch=gIntK[row1*M3+col1*M2+row2*M+col2];
                       gIntK[row1*M3+col1*M2+col2*M+row2] = scratch;
                       gIntK[col1*M3+row1*M2+row2*M+col2] = scratch;
                       gIntK[col1*M3+row1*M2+col2*M+row2] = scratch;
                    }
                 }
              } // END IF TRANS1
           }
        }
      } // END IF DENSITY
   }  // END points loop

   delete[] numeros;
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################
template<class scalar_type>
void PointGroupCPU<scalar_type>::get_coef_input(HostMatrix<scalar_type>&
               rmm_input,int* vecnum) const
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

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}
