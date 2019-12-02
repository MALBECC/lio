#include <iostream>
#include <cstdlib>
#include <math.h>

#define DENSMIN 1e-10

using namespace std;

void CouplingForm(double* K2,double* K4,int& M, int& NCO,
                  int& Nvirt,double* MatCoef)
{
//MatCoef is used as transpose matrix
   int a_MO, b_MO, j_MO, i_MO;
   int row, col, M3, M2, NCOC, dim;
   double finalk, Coef_MO1, Coef_MO2;

   finalk = Coef_MO1 = Coef_MO2 = 0.0;
   NCOC = NCO - 1;
   dim = NCO * Nvirt;
   M2 = M*M;
   M3 = M2*M;
   row = col = 0;
     
   for(row=0;row<dim;row++) {
      a_MO = row % Nvirt + NCO;
      i_MO = NCOC - row / Nvirt;
      for(col=row;col<dim;col++) {
         b_MO = col % Nvirt + NCO;
         j_MO = NCOC - col / Nvirt;
#pragma omp parallel for private(Coef_MO1,Coef_MO2), reduction(+:finalk)
         for(int i = 0; i < M; i++) {
            for(int j = 0; j < M; j++) {
               Coef_MO1 = MatCoef[M*a_MO+i] * MatCoef[M*i_MO+j];
               if(fabs(Coef_MO1) > DENSMIN) {
                 for(int k = 0; k < M; k++) {
                    for(int l = 0; l < M; l++) {
                       Coef_MO2 = Coef_MO1 * MatCoef[M*b_MO+k] * MatCoef[M*j_MO+l];
                       finalk += Coef_MO2 * K4[i*M3+j*M2+k*M+l];
                    }
                 }
               } // if COEF_MO1
            }
         } // Fin basis
         K2[row*dim+col] += finalk;
         finalk = 0.0;
      }
   }
}
