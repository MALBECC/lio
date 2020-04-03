#include <iostream>


void cioverlap(double*,double*,double*,double*,
               double*,double*,int*,double*,double*,int&,int&,int&,int&,
               int& ndets);

using namespace std;
extern "C" void g2g_cioverlap_(double* wfunc, double* wfunc_old,
        double* coef, double* coef_old, double* Sbig, double* sigma,
        int* kind_coupling, double* phases, double* phases_old,  
        int& M, int& NCO, int& Nvirt, int& nstates, int& ndets)
{
   cioverlap(wfunc,wfunc_old,coef,coef_old,Sbig,sigma,
             kind_coupling,phases,phases_old,
             M,NCO,Nvirt,nstates,ndets);
}

