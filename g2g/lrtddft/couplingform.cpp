#include <iostream>
#include <omp.h>
#include "../init.h"

#if FULL_DOUBLE
void CouplingForm(FortranMatrix<double>& Kxc2,FortranMatrix<double>& Kc2,
                  double* Kxc4,double* Kc4)
#else
void CouplingForm(FortranMatrix<float>& Kxc2,FortranMatrix<float>& Kc2,
                  float* Kxc4,float* Kc4)
#endif
{
   int row, col, M3, M2, M, NCO, Nvirt;
#if FULL_DOUBLE
   double final_xc, final_c, Coef_MO1, CoefMO2;
#else
   float final_xc, final_c, Coef_MO1, CoefMO2;
#endif

   final_xc = final_c = Coef_MO1 = CoefMO2 = 0.0;
   M = fortran_vars.m;
   NCO = fortran_vars.nco;
   Nvirt = M - NCO;
   M2 = M*M;
   M3 = M2*M;
   row = col = 0;

   double timeI, timeF;
   timeI = omp_get_wtime();
   for(int i_MO = NCO-1; i_MO >= 0; i_MO--) {
     for(int a_MO = 0; a_MO < Nvirt; a_MO++) {
       int fuk = a_MO;
       for(int j_MO = i_MO; j_MO >= 0; j_MO--) {
         for(int b_MO = fuk; b_MO < Nvirt; b_MO++) {
#pragma omp parallel for private(Coef_MO1,Coef_MO2), reduction(+:final_c,final_xc)
           for(int i = 0; i < M; i++) {
             for(int j = 0; j < M; j++) {
               Coef_MO1 = fortran_vars.MatCoef(i,a_MO+NCO) * fortran_vars.MatCoef(j,i_MO);
               if(fabs(Coef_MO1) > 1e-10) {
                 for(int k = 0; k < M; k++) {
                   for(int l = 0; l < M; l++) {
                     Coef_MO2 = Coef_MO1 * fortran_vars.MatCoef(k,b_MO+NCO)*fortran_vars.MatCoef(l,j_MO);
                     final_xc += Coef_MO2 * Kxc4[i*M3+j*M2+k*M+l];
                     final_c += Coef_MO2 * Kc4[i*M3+j*M2+k*M+l];
                   }
                 }
               } // if COEF_MO1
             }
           } // Fin basis
           Kxc2(row,col) += final_xc;
           Kc2(row,col) += final_c;
           final_xc = 0.0;
           final_c = 0.0;
           col += 1;
         }
         fuk = 0;
       }
       row += 1;
       col = row;
     }
   }
   timeF = omp_get_wtime();
   printf("COUPLINGFORM SUBROUTINE %f\n",timeF-timeI);
}
