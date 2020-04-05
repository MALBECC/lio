#include <xc.h>

int print_info(int func_id) 
{
   xc_func_type *func=(xc_func_type*)malloc(sizeof(xc_func_type));
   if(xc_func_init(func, func_id, XC_UNPOLARIZED) != 0) {
     fprintf(stderr, "Functional %d not found\n", func_id);
     return -1;
   }

   int value = 0;
   printf("  The functional '%s' is ", func->info->name);
   switch (func->info->kind) {
       case (XC_EXCHANGE):
           printf("an exchange functional");
       break;
       case (XC_CORRELATION):
           printf("a correlation functional");
       break;
       case (XC_EXCHANGE_CORRELATION):
           printf("an exchange-correlation functional");
       break;
       case (XC_KINETIC):
           printf("a kinetic energy functional");
       break;
       default:
           value = -1;
           printf("of unknown kind");
       break;
   }

   printf(",\n  it belongs to the ");
   switch (func->info->family) {
      case (XC_FAMILY_LDA):
          printf("LDA"); break;
      case (XC_FAMILY_GGA):
          printf("GGA"); break;
      case (XC_FAMILY_HYB_GGA):
          printf("Hybrid GGA"); break;
      case (XC_FAMILY_MGGA):
          printf("MGGA"); break;
      case (XC_FAMILY_HYB_MGGA):
          printf("Hybrid MGGA");; break;
      default:
          value = -1;
          printf("Family Unknown"); break;
   }

   printf("' family and is defined in the reference(s):\n");
   for (int ii = 0; func->info->refs[ii] != NULL; ii++) {
       printf ("  [%d] %s\n", ii+1, func->info->refs[ii]->ref);
   }
   xc_func_end(func);

   return value;
}


void set_pbe(int* HF, double* HF_fac, double* screen)
{
   int err = -1;
   fortran_vars.nx_func=1;
   fortran_vars.nc_func=1;
   fortran_vars.func_id  = (int*) malloc(sizeof(int)*2);
   fortran_vars.func_coef= (double*) malloc(sizeof(double)*2);

   // FUNCT
   fortran_vars.func_id[0] = XC_GGA_X_PBE;
   fortran_vars.func_id[1] = XC_GGA_C_PBE;
   fortran_vars.func_coef[0]=1.0f; // X
   fortran_vars.func_coef[1]=1.0f; // C
   fortran_vars.nsr_id = -1;
   fortran_vars.screen = -1.0f;

   *HF = 0; *HF_fac = 0.0f; *screen = 0.0f;
   err = print_info(XC_GGA_X_PBE);
   if ( err != 0 ) exit(-1);
   err = print_info(XC_GGA_C_PBE);
   if ( err != 0 ) exit(-1);
}

void set_pbe0(int* HF, double* HF_fac, double* screen)
{
   int err = -1;
   fortran_vars.nx_func=1;
   fortran_vars.nc_func=1;
   fortran_vars.func_id = (int*) malloc(sizeof(int)*2);
   fortran_vars.func_coef= (double*) malloc(sizeof(double)*2);

   // FUNCT
   fortran_vars.func_id[0] = XC_GGA_X_PBE;
   fortran_vars.func_id[1] = XC_GGA_C_PBE;
   fortran_vars.func_coef[0]=0.75f; // X
   fortran_vars.func_coef[1]=1.0f; // C
   fortran_vars.nsr_id = -1;
   fortran_vars.screen = -1.0f;

   *HF = 1; *HF_fac = 0.25f; *screen = 0.0f;
   err = print_info(XC_HYB_GGA_XC_PBEH);
   if ( err != 0 ) exit(-1);
}

void set_b3lyp(int* HF, double* HF_fac, double* screen)
{
   int err = -1;
   fortran_vars.nx_func=2;
   fortran_vars.nc_func=2;
   fortran_vars.func_id = (int*) malloc(sizeof(int)*4);
   fortran_vars.func_coef= (double*) malloc(sizeof(double)*4);
   
   // FUNCT
   fortran_vars.func_id[0] = XC_LDA_X;
   fortran_vars.func_id[1] = XC_GGA_X_B88;
   fortran_vars.func_id[2] = XC_LDA_C_VWN_RPA;
   fortran_vars.func_id[3] = XC_GGA_C_LYP;

   fortran_vars.func_coef[0]=0.08f; // X
   fortran_vars.func_coef[1]=0.72f; // X
   fortran_vars.func_coef[2]=0.19f; // C
   fortran_vars.func_coef[3]=0.81f; // C

   fortran_vars.nsr_id = -1;
   fortran_vars.screen = -1.0f;

   *HF = 1; *HF_fac = 0.20f; *screen = 0.0f;
   err = print_info(XC_HYB_GGA_XC_B3LYP);
   if ( err != 0 ) exit(-1);
}  

void set_cam_b3lyp(int* HF, double* HF_fac, double* screen)
{
   int err = -1;
   fortran_vars.nx_func=2;
   fortran_vars.nc_func=2;
   fortran_vars.func_id = (int*) malloc(sizeof(int)*4);
   fortran_vars.func_coef= (double*) malloc(sizeof(double)*4);

   // FUNCT
   fortran_vars.func_id[0] = XC_GGA_X_B88;
   fortran_vars.func_id[1] = XC_GGA_X_ITYH;
   fortran_vars.func_id[2] = XC_LDA_C_VWN;
   fortran_vars.func_id[3] = XC_GGA_C_LYP;

   fortran_vars.func_coef[0]=0.35f; // X
   fortran_vars.func_coef[1]=0.46f; // X
   fortran_vars.func_coef[2]=0.19f; // C
   fortran_vars.func_coef[3]=0.81f; // C

   fortran_vars.nsr_id = 1;
   fortran_vars.screen = 0.33f;

   *HF = 2; *HF_fac = 0.65f; *screen = 0.33f;
   err = print_info(XC_HYB_GGA_XC_CAM_B3LYP);
   if ( err != 0 ) exit(-1);
}
