#ifndef _COUPLINGFORM_H
#define _COUPLINGFORM_H
 #if FULL_DOUBLE
 void CouplingForm(FortranMatrix<double>&,FortranMatrix<double>&,
                   double*,double*)
 #else
 void CouplingForm(FortranMatrix<float>&,FortranMatrix<float>&,
                   float*,float*)
 #endif
#endif
 
