using namespace G2G;
#include "gpu_vxc.h"
#include "gpu_fxc.h"

template<class T, bool compute_energy, bool compute_factor, bool lda>
void gpu_calc_gradients(uint npoints,T* dens, T* tred, T* diff,
                        G2G::vec_type<T,WIDTH>* dxyz,
                        G2G::vec_type<T,WIDTH>* tredxyz,
                        G2G::vec_type<T,WIDTH>* diffxyz, 
                        T* ddum, T* vdum, T* pdum,
                        G2G::vec_type<T,WIDTH>* ddum_xyz,
                        G2G::vec_type<T,WIDTH>* vdum_xyz,
                        G2G::vec_type<T,WIDTH>* pdum_xyz, 
                        int calc_fxc)
{
// LIBXC INITIALIZATION
  fortran_vars.fexc = fortran_vars.func_coef[0];
#define libxc_init_param \
  fortran_vars.func_id, fortran_vars.func_coef, fortran_vars.nx_func, \
  fortran_vars.nc_func, fortran_vars.nsr_id, fortran_vars.screen, \
  XC_UNPOLARIZED
  LibxcProxy_cuda<T,4> libxcProxy_cuda(libxc_init_param);
#undef libxc_init_param
 
// LIBXC VARIABLES
   CudaMatrix<T> vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3;
   vrho.resize(COALESCED_DIMENSION(npoints));
   vsigma.resize(COALESCED_DIMENSION(npoints));
   v2rho2.resize(COALESCED_DIMENSION(npoints));
   v2rhosigma.resize(COALESCED_DIMENSION(npoints));
   v2sigma2.resize(COALESCED_DIMENSION(npoints));
   v3rho3.resize(COALESCED_DIMENSION(npoints));
   v3rho2sigma.resize(COALESCED_DIMENSION(npoints));
   v3rhosigma2.resize(COALESCED_DIMENSION(npoints));
   v3sigma3.resize(COALESCED_DIMENSION(npoints));

// CALL LIBXC
#define libxc_parameter \
   &libxcProxy_cuda, npoints, dens, dxyz, vrho.data, vsigma.data, v2rho2.data, v2rhosigma.data, \
   v2sigma2.data, v3rho3.data, v3rho2sigma.data, v3rhosigma2.data, v3sigma3.data
   libxc_gpu_derivs<T, true, true, false>(libxc_parameter);
#undef libxc_parameter

   int threadsPerBlock = 256;
   int blocksPerGrid = (npoints + threadsPerBlock - 1) / threadsPerBlock;
#define VXC_parameter \
   npoints, diff, dxyz, diffxyz, vrho.data, vsigma.data, v2rho2.data, \
   v2rhosigma.data, v2sigma2.data, ddum, pdum, ddum_xyz, pdum_xyz, calc_fxc
   gpu_calc_VXC<T,true,true,false><<<blocksPerGrid,threadsPerBlock>>>(VXC_parameter);
#undef VXC_parameter

   if ( calc_fxc == 0 ) {
      #define FXC_parameter \
         npoints, dens, tred, dxyz, tredxyz, vsigma.data, v2rho2.data, v2rhosigma.data, \
         v2sigma2.data, v3rho3.data, v3rho2sigma.data, v3rhosigma2.data, v3sigma3.data, \
         ddum, vdum, ddum_xyz, vdum_xyz
         gpu_calc_FXC<T,true,true,false><<<blocksPerGrid,threadsPerBlock>>>(FXC_parameter);
      #undef FXC_parameter
   }
   vrho.deallocate(); vsigma.deallocate();
   v2rho2.deallocate(), v2rhosigma.deallocate(); v2sigma2.deallocate();
   v3rho3.deallocate(), v3rho2sigma.deallocate(), v3rhosigma2.deallocate(), v3sigma3.deallocate();

/*
   // DEBUG //
   T* rho_cpu = (T*) malloc(sizeof(T) * npoints);
   cudaMemcpy(rho_cpu,diff,npoints*(sizeof(T)),cudaMemcpyDeviceToHost);

   HostMatrix< vec_type<T,4> > var_cpu;
   var_cpu.resize(COALESCED_DIMENSION(npoints)); var_cpu.zero();
   cudaMemcpy(var_cpu.data,diffxyz,npoints*(sizeof(vec_type<T,4>)),cudaMemcpyDeviceToHost);

   for(int i=0; i<npoints; i++) {
      printf("%f %f %f %f\n",rho_cpu[i], var_cpu(i).x, var_cpu(i).y, var_cpu(i).z);

   }
   exit(-1);
*/
}

