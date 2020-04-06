#include <iostream>

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "init.h"
#include "partition.h"

using namespace G2G;


#if USE_LIBXC
#include "functional.h"

extern "C" void g2g_extern_functional_(int& main_id, bool* externFunc,
                              int* HF, double* HF_fac, double* screen)
{
   cout << " " << endl;
   cout << " Extern Functional Module " << endl;

   // Free Memory
   if ( fortran_vars.func_id != NULL ) {
      free(fortran_vars.func_id); fortran_vars.func_id = NULL;
   }
   if ( fortran_vars.func_coef != NULL ) {
      free(fortran_vars.func_coef); fortran_vars.func_coef = NULL;
   }
   if ( fortran_vars.HF != NULL ) {
      free(fortran_vars.HF); fortran_vars.HF = NULL;
   }
   if ( fortran_vars.HF_fac != NULL ) {
      free(fortran_vars.HF_fac); fortran_vars.HF_fac = NULL;
   }
   
   switch (main_id) {
      // PBE
      case 101:
           set_pbe(HF,HF_fac,screen); break;

      case 130:
           set_pbe(HF,HF_fac,screen); break;

      // PBE0
      case 406:
           set_pbe0(HF,HF_fac,screen); break;

      // B3LYP
      case 402:
           set_b3lyp(HF,HF_fac,screen); break;

      // CAM-B3LYP
      case 433:
           set_cam_b3lyp(HF,HF_fac,screen); break;

      // LC-WPBE
      case 478:
           set_lc_wpbe(HF,HF_fac,screen); break;

      // LC-BLYP
      case 400:
           set_lc_blyp(HF,HF_fac,screen); break;
   
      default:
           cout << "The Functional id " << main_id << " doesn't implemented yet" << endl;
           exit(-1); break;
   }
   cout << " " << endl;
}
#else
extern "C" void g2g_extern_functional_(int& main_id, bool* externFunc,
                              int* HF, double* HF_fac, double* screen)
{
   cout << " " << endl;
   cout << " Extern Functional Module " << endl;

   if ( main_id == 406 ) {
      fortran_vars.fexc = 0.75f;
      *HF = 1; *HF_fac = 0.25f; *screen = 0.0f;
   } else {
      cout << "In order to use external Functional you need to recompile ";
      cout << "LIO with libxc=1 or 2" << endl;
      exit(-1);
   }
}
#endif
