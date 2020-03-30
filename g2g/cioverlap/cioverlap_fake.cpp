#include <iostream>


using namespace std;
extern "C" void g2g_cioverlap_(double* wfunc, double* wfunc_old,
        double* coef, double* coef_old, double* Sbig, double* sigma,
        int* kind_coupling, double* phases, double* phases_old,  
        int& M, int& NCO, int& Nvirt, int& nstates, int& ndets)
{
   cout << " In order to run TSH dynamics, you have to compile LIO " ;
   cout << "with tsh=1 option" << endl;
   exit(-1);
}

