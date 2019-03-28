#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math_constants.h>
#include <float.h>
#include <string>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>

#include "../../../g2g/common.h"
#include "../../../g2g/init.h"
#include "../../../g2g/cuda/cuda_extra.h"
#include "../../../g2g/matrix.h"
#include "../../../g2g/timer.h"
#include "../../../g2g/partition.h"
#include "../../../g2g/scalar_vector_types.h"
#include "../../../g2g/global_memory_pool.h"
#include "../../../g2g/libxc/libxcproxy.h"
#include "../../../g2g/libxc/libxc_accumulate_point.h"

#include "../../unit_test/commons/test_input.h"

//////////////////////////////////////
//// CALC_GGA
#define POT_ALPHA     ((double)-0.738558766382022447)
#define POT_GL        ((double)0.620350490899400087)

#define POT_VOSKO_A1  ((double)0.03109205)
#define POT_VOSKO_B1  ((double)3.72744)
#define POT_VOSKO_C1  ((double)12.9352)
#define POT_VOSKO_X0  ((double)-0.10498)
#define POT_VOSKO_Q   ((double)6.15199066246304849)
#define POT_VOSKO_A16 ((double)0.005182008333)
#define POT_VOSKO_Q2  ((double)4.7309269)

#define POT_ALYP  ((double)0.04918)
#define POT_BLYP  ((double)0.132)
#define POT_CLYP  ((double)0.2533)
#define POT_CLYP3 ((double)0.0844333333)
#define POT_DLYP  ((double)0.349)
#define POT_DLYP3 ((double)0.116333333)
#define POT_CF    ((double)2.87123400018819)
#define POT_BETA  ((double)0.0042)

#define POT_ALF ((double)0.023266)
#define POT_BET ((double)7.389)
#define POT_GAM ((double)8.723)
#define POT_DEL ((double)0.472)

//////////////////////////////////////
//// CALC_GGACS
//
// Esta funciona la utilizamos para simular llamadas
// a Lio y comparar con los resultados que obtenemos
// desde Libxc.
template<class T, int iexch, unsigned int width>  __device__
void calc_ggaCS( T dens, 
                 const G2G::vec_type<T,width>& grad,
                 const G2G::vec_type<T,width>& hess1,
                 const G2G::vec_type<T,width>& hess2,
                 T& ex, 
                 T& ec, 
                 T& y2a)
{
   // hess1: xx, yy, zz  || hess2: xy, xz, yz
   const T MINIMUM_DENSITY_VALUE = 1e-13f;
   if (dens < MINIMUM_DENSITY_VALUE) { ex = ec = 0; return; }

   T y     = cbrt( (T)dens );  // rho^(1/3)
   T grad2 = grad.x * grad.x + grad.y * grad.y + grad.z * grad.z;
   if (grad2 == 0) grad2 = (T)FLT_MIN;
   T dgrad = sqrt(grad2);

   T d0 = hess1.x + hess1.y + hess1.z;
   T u0 = ((grad.x * grad.x) * hess1.x 
                  + 2.0 * grad.x * grad.y * hess2.x 
                  + 2.0 * grad.y * grad.z * hess2.z 
                  + 2.0 * grad.x * grad.z * hess2.y 
                  + (grad.y * grad.y) * hess1.y 
                  + (grad.z * grad.z) * hess1.z) / dgrad;
   y2a = 0;

   // Exchange - Perdew : Phys. Rev B 33 8800 (1986)
   if (iexch == 4 || iexch == 8) {
      T dens2 = (dens * dens);
      T ckf   = (T)3.0936677 * y;
      T s     = dgrad / ((T)2.0 * ckf * dens);

      T fx = (1.0 / 15.0);
      T s2 = (s * s);
      T s3 = (s * s * s);
      T g0 = 1.0 + 1.296 * s2 + 14.0 * pow(s, 4) + 0.2 * pow(s, 6);
      T F  = pow(g0, fx);
      T e  = POT_ALPHA * F * y;
      ex = e;

      T t = d0 / (dens * 4.0 * (ckf * ckf));
      T u = u0 / (pow( (T)2.0 * ckf, 3) * dens2);

      T g2  = 2.592 * s + 56.0 * s3 + 1.2 * pow(s, 5);
      T g3  = 2.592     + 56.0 * s2 + 1.2 * pow(s, 4);
      T g4  = 112.0 * s + 4.8  * s3;
      T dF  = fx * F/g0 * g2;
      T dsF = fx * F/g0 * (-14.0 * fx * g3 * g2/g0 + g4);

      y2a = POT_ALPHA * y * (1.33333333333 * F - t/s * dF - (u-1.3333333333 * s3) * dsF);
   } else if (iexch >= 5 && iexch <= 7) { // Becke  : Phys. Rev A 38 3098 (1988)
      T e0 = POT_ALPHA * y;
      T y2 = dens / 2.0;
      T r13 = cbrt( y2 );
      T r43 = cbrt( pow(y2, 4) );
      T Xs = dgrad / (2.0 * r43);
      T siper = asinh(Xs);
      T DN = 1.0 + 6.0 * POT_BETA * Xs * siper;
      T ect = -2.0 * POT_BETA * r43 * Xs * Xs/(DN * dens);
      T e = e0 + ect;
      ex = e;

      // Potential
      T v0 = 1.33333333333333 * e0;
      T Fb = 1.0 / DN;
      T XA1 = Xs / sqrt(1.0 + Xs * Xs);
      T DN1 = 1.0 + Fb * (1.0 - 6.0 * POT_BETA * Xs * XA1);
      T DN2 = 1.0 / (1.0 + Xs * Xs) + 2.0 * Fb * (2.0 - 6.0 * POT_BETA * Xs * XA1);
      T DN3 = siper * (1.0 + 2.0 * Fb) + XA1 * DN2;
      T D02 = d0 / 2.0;
      T de1 = 1.33333333333333 / (cbrt(pow((T)dens,7) ) );

      T DGRADx = (grad.x * hess1.x + grad.y * hess2.x + grad.z * hess2.y) / dgrad;
      T GRADXx = cbrt( (T) 2.0 ) * (1.0 / (dens * y) * DGRADx - de1 * grad.x * dgrad);
      T DGRADy = (grad.x * hess2.x + grad.y * hess1.y + grad.z * hess2.z) / dgrad;
      T GRADXy = cbrt( (T) 2.0 ) * (1.0 / (dens * y) * DGRADy - de1 * grad.y * dgrad);
      T DGRADz = (grad.x * hess2.y + grad.y * hess2.z + grad.z * hess1.z) / dgrad;
      T GRADXz = cbrt( (T) 2.0 ) * (1.0 / (dens * y) * DGRADz - de1 * grad.z * dgrad);

      T T1   = grad.x / 2.0 * GRADXx;
      T T2   = grad.y / 2.0 * GRADXy;
      T T3   = grad.z / 2.0 * GRADXz;
      T DN4  = 6.0 * POT_BETA * Fb * (T1 + T2 + T3);
      T DN5  = 1.33333333333333 * r43 * r13 * Xs * Xs;
      T TOT2 = DN5 - D02 * DN1 + DN4 * DN3;

      T vxc = -POT_BETA * Fb/r43 * TOT2;
      y2a = v0 + vxc;
   } else { // PBE
      //T dgrad2 = grad.y * grad.y * hess1.y;
      //T dgrad3 = grad.z * grad.z * hess1.z;
      //T dgrad4 = grad.x * grad.y * hess2.x;
      //T dgrad5 = grad.x * grad.z * hess2.y;
      //T dgrad6 = grad.y * grad.z * hess2.z;
      //T delgrad = (dgrad1 + dgrad2 + dgrad3 + 2 * (dgrad4 + dgrad5 + dgrad6)) / dgrad;
      //T rlap = hess1.x + hess1.y + hess1.z;

      //T expbe, vxpbe, ecpbe, vcpbe;
      //pbeCS(dens, dgrad, delgrad, rlap, expbe, vxpbe, ecpbe, vcpbe);

      //ex  = expbe;
      //ec  = ecpbe;
      //y2a = vxpbe + vcpbe;
      return;
   }

   // Correlation - Perdew : Phys. Rev B 33 8822 (1986)
   if (iexch >= 4 && iexch <= 6) {
      // TO-DO: hay algun problema con 4 y 5, probablemente este aca
      T dens2 = (dens * dens);
      T rs  = POT_GL / y;
      T x1  = sqrt(rs);
      T Xx  = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
      T Xxo = (POT_VOSKO_X0 * POT_VOSKO_X0) 
                      + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1;
  
      T t1 = 2.0 * x1 + POT_VOSKO_B1;
      T t2 = log(Xx);
      T t3 = atan(POT_VOSKO_Q/t1);
      T t4 = POT_VOSKO_B1 * POT_VOSKO_X0/Xxo;

      ec = POT_VOSKO_A1 * ( 2.0 * log(x1) - t2 
           + 2.0 * POT_VOSKO_B1/POT_VOSKO_Q * t3
           - t4 *(2.0 * log(x1 - POT_VOSKO_X0) - t2 
           + 2.0 * (POT_VOSKO_B1 + 2.0 * POT_VOSKO_X0) / POT_VOSKO_Q * t3));

      T t5 = (POT_VOSKO_B1 * x1 + 2.0 * POT_VOSKO_C1) / x1;
      T t6 = POT_VOSKO_X0 / Xxo;
      T vc = ec - POT_VOSKO_A16 * x1 * 
                   ( t5/Xx - 4.0 * POT_VOSKO_B1 / ((t1 * t1)+(POT_VOSKO_Q * POT_VOSKO_Q2)) 
                   * (1.0 - t6 * (POT_VOSKO_B1 - 2.0 * POT_VOSKO_X0)) 
                   - t4 * (2.0 / (x1 - POT_VOSKO_X0) - t1/Xx));

      if (iexch == 6) {
         y2a = y2a + vc;
      } else {
         T rs2 = (rs * rs);
         T Cx1 = 0.002568f + POT_ALF * rs + POT_BET * rs2;
         T Cx2 = 1.0f + POT_GAM * rs + POT_DEL * rs2 + 1.0e4 * POT_BET * (rs * rs * rs);
         T C   = 0.001667 + Cx1/Cx2;
         T Cx3 = POT_ALF + 2.0f * POT_BET * rs;
         T Cx4 = POT_GAM + 2.0f * POT_DEL * rs + 3.0e4 * POT_BET * rs2;
         T dC  = Cx3/Cx2 - Cx1/(Cx2 * Cx2) * Cx4;
         dC = -0.333333333333333f * dC * POT_GL / (y * dens);

         T phi  = 0.0008129082f/C * dgrad/pow((T)dens, (T)(7.0f/6.0f));
         T expo = exp(-phi);
         T ex0  = expo * C;

         ec = ec + ex0 * grad2 / (y * dens2);

         T D1   = (2.0f - phi) * d0/dens;
         T phi2 = (phi * phi);
         T D2   = 1.33333333333333333f - 3.666666666666666666f * phi + 1.166666666666666f * phi2;
         D2 = D2 * grad2/dens2;
         T D3 = phi * (phi - 3.0f) * u0/(dens * dgrad);
         T D4 = expo * grad2 / (y * dens) * (phi2 - phi - 1.0f) * dC;

         vc = vc - 1.0 * (ex0 / y * (D1 - D2 + D3) - D4);
         y2a = y2a + vc;
      }
   } else if (iexch == 7 || iexch == 8) { // Correlation - LYP: PRB 37 785 (1988)
      T rom13 = 1 / cbrt( dens );
      T rom53 = cbrt( pow(dens, 5) );
      T ecro  = expf(-POT_CLYP * rom13);
      T f1    = 1.0f / (1.0f + POT_DLYP * rom13);
      T tw    = 1.0f / 8.0f * (grad2/dens - d0);
      T term  = (tw / 9.0f + d0 / 18.0f) - 2.0f * tw + POT_CF * rom53;
      term = dens + POT_BLYP * (rom13 * rom13) * ecro * term;

      ec = -POT_ALYP * f1 * term/dens;

      T h1 = ecro/rom53;
      T g1 = f1 * h1;
      T tm1 = POT_DLYP3 * (rom13/dens);
      T fp1 = tm1 * (f1 * f1);
      T tm2 = -1.666666666f + POT_CLYP3 * rom13;
      T hp1 = h1 * tm2/dens;
      T gp1 = fp1 * h1 + hp1 * f1;
      T fp2 = tm1 * 2.0f * f1 * (fp1 - 0.6666666666f * f1/dens);
      T tm3 = 1.6666666666f - POT_CLYP3 * 1.3333333333f * rom13;
      T hp2 = hp1 * tm2/dens + h1 * tm3/(dens * dens);
      T gp2 = fp2 * h1 + 2.0f * fp1 * hp1 + hp2 * f1;

      T term3 = -POT_ALYP * (fp1 * dens + f1) 
                          -POT_ALYP * POT_BLYP * POT_CF * (gp1 * dens + 8.0f/3.0f * g1) * rom53;
      T term4 = (gp2 * dens * grad2 + gp1 * (3.0f * grad2 + 2.0f * dens * d0)
                          + 4.0f * g1 * d0) * POT_ALYP * POT_BLYP/4.0f;
      T term5 = (3.0f * gp2 * dens * grad2 + gp1 * (5.0f * grad2 + 6.0f * dens * d0)
                          + 4.0f * g1 * d0) * POT_ALYP * POT_BLYP/72.0f;

      y2a = y2a + (term3 - term4 - term5);
   }
}

template<class T, unsigned int width> __device__
void calc_ggaCS_in( T dens, 
                    const G2G::vec_type<T,width>& grad,
                    const G2G::vec_type<T,width>& hess1,
                    const G2G::vec_type<T,width>& hess2,
                    T& ex, 
                    T& ec, 
                    T& y2a,
                    const int iexch)
{
   switch(iexch) {
      case 0: return calc_ggaCS<T, 0, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      case 1: return calc_ggaCS<T, 1, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      case 2: return calc_ggaCS<T, 2, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      case 3: return calc_ggaCS<T, 3, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      case 4: return calc_ggaCS<T, 4, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      case 5: return calc_ggaCS<T, 5, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      case 6: return calc_ggaCS<T, 6, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      case 7: return calc_ggaCS<T, 7, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      case 8: return calc_ggaCS<T, 8, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      case 9: return calc_ggaCS<T, 9, width>(dens, grad, hess1, hess2, ex, ec, y2a);
      default: assert(false);
   }
}


//////////////////////////////////////
// KERNEL FOR ACCUMULATE_POINT
//
template<class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_accumulate_point (T* const energy, T* const factor, 
		    const T* const point_weights,
            	    uint points, int block_height, T* partial_density, 
		    G2G::vec_type<T,4>* dxyz,
                    G2G::vec_type<T,4>* dd1, 
		    G2G::vec_type<T,4>* dd2) {

  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;
  //printf("point:%i\n", point);

  T point_weight = 0.0f;
  T y2a, exc_corr, exc_c, exc_x;

  T _partial_density(0.0f);
  G2G::vec_type<T,4> _dxyz, _dd1, _dd2;

  _dxyz = _dd1 = _dd2 = G2G::vec_type<T,4>(0.0f,0.0f,0.0f,0.0f);

  bool valid_thread = (point < points);
  if (valid_thread)
    point_weight = point_weights[point];

  if (valid_thread) {
    for(int j =0 ; j<block_height; j++) {
      const int this_row = j*points+point;

      _partial_density += partial_density[this_row];
      _dxyz += dxyz[this_row];
      _dd1 += dd1[this_row];
      _dd2 += dd2[this_row];
     }
  }

  calc_ggaCS_in<T, 4> (_partial_density, _dxyz, _dd1, _dd2, exc_x, exc_c, y2a, 9);
  exc_corr = exc_x + exc_c;

  if (compute_energy && valid_thread){
    energy[point] = (_partial_density * point_weight) * exc_corr;
    //printf("exc_corr:%e\n", exc_corr);
    //printf("_partial_density:%e\n", _partial_density);
    //printf("point_weight:%e\n", point_weight);
    //printf("energy[%i]:%e\n", point, energy[point]);
  }

  if (compute_factor && valid_thread){
    factor[point] = point_weight * y2a;
    //printf("y2a:%e\n", y2a);
    //printf("factor[%i]:%e\n", point, point_weight);
  }

}


//////////////////////////////////////
//// TESTS

/////////////////////////////////////////////////////
// Test: accumulate_data_for_libxc_test0001
//
void gpu_accumulate_point_test0001()
{
    printf("** gpu_accumulate_point_test0001 **\n");

    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu;

    dxyz_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu.resize(COALESCED_DIMENSION(number_of_points),m);

    // Now the arrays for energy, factors, point_weight and partial_density
    double *energy_gpu = NULL;
    double *factor_gpu = NULL;
    double *point_weights_gpu = NULL;
    double *partial_density_gpu = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);
    err = cudaMalloc((void**)&energy_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu, 0, size);
    cudaMemset(factor_gpu, 0, size);
    cudaMemset(point_weights_gpu, 1, size);
    cudaMemset(partial_density_gpu, 1, size);

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    // Call the CUDA KERNEL
    gpu_accumulate_point<double,true, true, false><<<blocksPerGrid, threadsPerBlock>>> 
		    (energy_gpu, factor_gpu, 
		    point_weights_gpu,
            	    numElements, 1, partial_density_gpu, 
		    dxyz_gpu.data,
                    dd1_gpu.data, 
		    dd2_gpu.data);


    // Allocate the host input vectors
    double *energy_cpu = (double *)malloc(size);
    double *factor_cpu = (double *)malloc(size);
    double *point_weights_cpu = (double *)malloc(size);
    double *partial_density_cpu = (double *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu, energy_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu, factor_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(point_weights_cpu, point_weights_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_cpu, partial_density_gpu, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Before print the results, get the data from cuda
    G2G::vec_type<double,4>* dxyz_cpu;
    G2G::vec_type<double,4>* dd1_cpu;
    G2G::vec_type<double,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint cpu_vecs_size = number_of_points * sizeof(G2G::vec_type<double,4>);
    dxyz_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);
    dd1_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);
    dd2_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    print_accumulate_point_data (dxyz_cpu, dd1_cpu, dd2_cpu, energy_cpu, 
	factor_cpu, point_weights_cpu, 
	partial_density_cpu, number_of_points);

    free(energy_cpu);
    free(factor_cpu);
    free(point_weights_cpu);
    free(partial_density_cpu);

    free(dxyz_cpu);
    free(dd1_cpu);
    free(dd2_cpu);

    cudaFree(energy_gpu);
    cudaFree(factor_gpu);
    cudaFree(point_weights_gpu);
    cudaFree(partial_density_gpu);

}

/////////////////////////////////////////////////////
// Test: accumulate_data_for_libxc_test0001
//
void accumulate_data_for_libxc_test0001()
{
    printf("** accumulate_data_for_libxc_test0001 **\n");

    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in;

    // Accumulate
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_accum;

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);

    dxyz_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);

    // Now the arrays for energy, factors, point_weight and partial_density
    double *point_weights_gpu_in = NULL;
    double *partial_density_gpu_in = NULL;
    // Accum
    double *partial_density_gpu_accum = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);
    err = cudaMalloc((void**)&point_weights_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_accum, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_accum!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(point_weights_gpu_in, 1, size);
    cudaMemset(partial_density_gpu_in, 1, size);
    cudaMemset(partial_density_gpu_accum, 1, size);

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    uint block_height = 1;

    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

    // Call the CUDA KERNEL
    gpu_accumulate_point_for_libxc<double,true,true,false><<<blocksPerGrid,threadsPerBlock>>> (point_weights_gpu_in,
	number_of_points, block_height,
	partial_density_gpu_in, dxyz_gpu_in.data, dd1_gpu_in.data, dd2_gpu_in.data,
	partial_density_gpu_accum, dxyz_gpu_accum.data, dd1_gpu_accum.data, dd2_gpu_accum.data);

    // Free memory
    cudaFree (point_weights_gpu_in);
    cudaFree (partial_density_gpu_in);
    cudaFree (partial_density_gpu_accum);

}

/////////////////////////////////////////////////////
// Test: accumulate_data_for_libxc_test0002
//
void accumulate_data_for_libxc_test0002()
{
    printf("** accumulate_data_for_libxc_test0002 **\n");

    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in;

    // Accumulate
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_accum;

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);

    dxyz_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);

    // Now the arrays for energy, factors, point_weight and partial_density
    double *point_weights_gpu_in = NULL;
    double *partial_density_gpu_in = NULL;
    // Accum
    double *partial_density_gpu_accum = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);
    err = cudaMalloc((void**)&point_weights_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_accum, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_accum!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(point_weights_gpu_in, 1, size);
    cudaMemset(partial_density_gpu_in, 1, size);
    cudaMemset(partial_density_gpu_accum, 1, size);

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    uint block_height = 1;

    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

    // Call the CUDA KERNEL
    gpu_accumulate_point_for_libxc<double,true,true,false><<<blocksPerGrid,threadsPerBlock>>> (point_weights_gpu_in,
	number_of_points, block_height,
	partial_density_gpu_in, dxyz_gpu_in.data, dd1_gpu_in.data, dd2_gpu_in.data,
	partial_density_gpu_accum, dxyz_gpu_accum.data, dd1_gpu_accum.data, dd2_gpu_accum.data);

    // With the accumulate parameters we can use
    // the data for libxc.
    G2G::vec_type<double,4>* dxyz_cpu;
    G2G::vec_type<double,4>* dd1_cpu;
    G2G::vec_type<double,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint cpu_vecs_size = number_of_points * sizeof(G2G::vec_type<double,4>);
    dxyz_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);
    dd1_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);
    dd2_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu_accum.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu_accum from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu_accum.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu_accum from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu_accum.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu_accum from device to host!\n");
        exit(EXIT_FAILURE);
    }

    print_accumulate_point_data (dxyz_cpu, dd1_cpu, dd2_cpu, NULL, 
	NULL, NULL, NULL, number_of_points);

    // Free memory
    free(dxyz_cpu);
    free(dd1_cpu);
    free(dd2_cpu);

    cudaFree (point_weights_gpu_in);
    cudaFree (partial_density_gpu_in);
    cudaFree (partial_density_gpu_accum);

}

/////////////////////////////////////////////////////
// Test: accumulate_data_for_libxc_test0003
//
void accumulate_data_for_libxc_test0003()
{
    printf("** accumulate_data_for_libxc_test0003 **\n");
#if FULL_DOUBLE
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in;

    // Accumulate
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_accum;

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);

    dxyz_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);

    // Now the arrays for energy, factors, point_weight and partial_density
    double *point_weights_gpu_in = NULL;
    double *partial_density_gpu_in = NULL;
    // Accum
    double *accumulated_density_gpu = NULL;

    // Now the arrays for energy, factors, point_weight and partial_density
    double *energy_gpu_in = NULL;
    double *factor_gpu_in = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);

    err = cudaMalloc((void**)&energy_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_in!\n");
    }

    err = cudaMalloc((void**)&accumulated_density_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_accum!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu_in, 0, size);
    cudaMemset(factor_gpu_in, 0, size);
    cudaMemset(point_weights_gpu_in, 1, size);
    cudaMemset(partial_density_gpu_in, 1, size);
    cudaMemset(accumulated_density_gpu, 1, size);

    // Launch the CUDA Kernel
    int numElements = number_of_points;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    uint block_height = 1;

    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

    // Call the CUDA KERNEL
    gpu_accumulate_point_for_libxc<double,true,true,false><<<blocksPerGrid,threadsPerBlock>>> (point_weights_gpu_in,
	number_of_points, block_height,
	partial_density_gpu_in, dxyz_gpu_in.data, dd1_gpu_in.data, dd2_gpu_in.data,
	accumulated_density_gpu, dxyz_gpu_accum.data, dd1_gpu_accum.data, dd2_gpu_accum.data);

    // Create the libxcproxy
    const int nspin = 1;
    const int functionalExchange = 101;
    const int functionalCorrelation = 130;
    LibxcProxy<double,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    // Calculate exc_corr and y2a
    libxc_exchange_correlation_cpu<double, true, true, false> (&libxcProxy,
	energy_gpu_in,
	factor_gpu_in,
	number_of_points,
	accumulated_density_gpu,
	dxyz_gpu_accum.data,
        dd1_gpu_accum.data,
	dd2_gpu_accum.data);

    // TODO: queda unir los resultados.
    // con la funcion gpu_accumulate_energy_and_factor_from_libxc(...);

    // Free memory
    cudaFree(energy_gpu_in);
    cudaFree(factor_gpu_in);
    cudaFree(point_weights_gpu_in);
    cudaFree(partial_density_gpu_in);
    cudaFree(accumulated_density_gpu);
#endif
}


/////////////////////////////////////////////////////
// Test: accumulate_data_for_libxc_test0004
//
void accumulate_data_for_libxc_test0004()
{
    printf("** accumulate_data_for_libxc_test0004 **\n");
#if FULL_DOUBLE
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in;

    // Accumulate
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_accum;

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);

    dxyz_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);

    // Now the arrays for energy, factors, point_weight and partial_density
    double *point_weights_gpu_in = NULL;
    double *partial_density_gpu_in = NULL;
    // Accum
    double *accumulated_density_gpu = NULL;

    // Now the arrays for energy, factors, point_weight and partial_density
    double *energy_gpu_in = NULL;
    double *factor_gpu_in = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);

    err = cudaMalloc((void**)&energy_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_in!\n");
    }

    err = cudaMalloc((void**)&accumulated_density_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_accum!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu_in, 0, size);
    cudaMemset(factor_gpu_in, 0, size);
    cudaMemset(point_weights_gpu_in, 1, size);
    cudaMemset(partial_density_gpu_in, 1, size);
    cudaMemset(accumulated_density_gpu, 1, size);

    // Launch the CUDA Kernel
    int numElements = number_of_points;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    uint block_height = 1;

    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

    // Call the CUDA KERNEL
    gpu_accumulate_point_for_libxc<double,true,true,false><<<blocksPerGrid, threadsPerBlock>>> (point_weights_gpu_in,
	number_of_points, block_height,
	partial_density_gpu_in, dxyz_gpu_in.data, dd1_gpu_in.data, dd2_gpu_in.data,
	accumulated_density_gpu, dxyz_gpu_accum.data, dd1_gpu_accum.data, dd2_gpu_accum.data);

    // Create the libxcproxy
    const int nspin = 1;
    const int functionalExchange = 101;
    const int functionalCorrelation = 130;
    LibxcProxy<double,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    // Calculate exc_corr and y2a
    libxc_exchange_correlation_cpu<double, true, true, false> (&libxcProxy,
	energy_gpu_in,
	factor_gpu_in,
	number_of_points,
	accumulated_density_gpu,
	dxyz_gpu_accum.data,
        dd1_gpu_accum.data,
	dd2_gpu_accum.data);

    // Join the results.
    gpu_accumulate_energy_and_forces_from_libxc<double, true, true, false><<<blocksPerGrid, threadsPerBlock>>> (
	energy_gpu_in,
	factor_gpu_in,
	point_weights_gpu_in,
	number_of_points,
	accumulated_density_gpu);

    // Free memory
    cudaFree(energy_gpu_in);
    cudaFree(factor_gpu_in);
    cudaFree(point_weights_gpu_in);
    cudaFree(partial_density_gpu_in);
    cudaFree(accumulated_density_gpu);
#endif
}

/////////////////////////////////////////////////////
// Test: accumulate_data_for_libxc_test0005
//
void accumulate_data_for_libxc_test0005()
{
    printf("** accumulate_data_for_libxc_test0005 **\n");
#if FULL_DOUBLE
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in;

    // Accumulate
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_accum;

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);

    dxyz_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);

    // Now the arrays for energy, factors, point_weight and partial_density
    double *point_weights_gpu_in = NULL;
    double *partial_density_gpu_in = NULL;
    // Accum
    double *accumulated_density_gpu = NULL;

    // Now the arrays for energy, factors, point_weight and partial_density
    double *energy_gpu_in = NULL;
    double *factor_gpu_in = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);

    err = cudaMalloc((void**)&energy_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_in!\n");
    }

    err = cudaMalloc((void**)&accumulated_density_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_accum!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu_in, 0, size);
    cudaMemset(factor_gpu_in, 0, size);
    cudaMemset(point_weights_gpu_in, 1, size);
    cudaMemset(partial_density_gpu_in, 1, size);
    cudaMemset(accumulated_density_gpu, 1, size);

    /////////////////////////
    // Create the libxcproxy
    const int nspin = 1;
    const int functionalExchange = 101;
    const int functionalCorrelation = 130;
    LibxcProxy<double,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    /////////////////////////
    // Calculate exc_corr and y2a
    libxc_exchange_correlation_cpu<double, true, true, false> (&libxcProxy,
	energy_gpu_in,
	factor_gpu_in,
	number_of_points,
	accumulated_density_gpu,
	dxyz_gpu_accum.data,
        dd1_gpu_accum.data,
	dd2_gpu_accum.data);

    /////////////////////////
    // Check and print the results.
    // Copy back the results before print.
    // Allocate the host input vectors
    double *energy_cpu = (double *)malloc(size);
    double *factor_cpu = (double *)malloc(size);
    double *point_weights_cpu = (double *)malloc(size);
    double *partial_density_cpu = (double *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu, energy_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu, factor_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(point_weights_cpu, point_weights_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_cpu, partial_density_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Before print the results, get the data from cuda
    G2G::vec_type<double,4>* dxyz_cpu;
    G2G::vec_type<double,4>* dd1_cpu;
    G2G::vec_type<double,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint cpu_vecs_size = number_of_points * sizeof(G2G::vec_type<double,4>);
    dxyz_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);
    dd1_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);
    dd2_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu_in.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu_in.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu_in.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    print_accumulate_point_data (dxyz_cpu, dd1_cpu, dd2_cpu, energy_cpu, 
	factor_cpu, NULL, NULL, number_of_points);


    ////////////////////////////
    // Free Memory
    free(energy_cpu);
    free(factor_cpu);
    free(point_weights_cpu);
    free(partial_density_cpu);

    free(dxyz_cpu);
    free(dd1_cpu);
    free(dd2_cpu);

    cudaFree(energy_gpu_in);
    cudaFree(factor_gpu_in);
    cudaFree(point_weights_gpu_in);
    cudaFree(partial_density_gpu_in);
    cudaFree(accumulated_density_gpu);
#endif
}

/////////////////////////////////////////////////////
// Test: accumulate_data_for_libxc_test0006
//
void accumulate_data_for_libxc_test0006()
{
    printf("** accumulate_data_for_libxc_test0006() **\n");
#if FULL_DOUBLE
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in;

    // Accumulate
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_accum;
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_accum;

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);

    dxyz_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);

    // Now the arrays for energy, factors, point_weight and partial_density
    double *point_weights_gpu_in = NULL;
    double *partial_density_gpu_in = NULL;
    // Accum
    double *accumulated_density_gpu = NULL;

    // Now the arrays for energy, factors, point_weight and partial_density
    double *energy_gpu_in = NULL;
    double *factor_gpu_in = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);

    err = cudaMalloc((void**)&energy_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    err = cudaMalloc((void**)&point_weights_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_in!\n");
    }

    err = cudaMalloc((void**)&accumulated_density_gpu, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_accum!\n");
    }

    // Set the cuda array values to a default value.
    cudaMemset(energy_gpu_in, 0, size);
    cudaMemset(factor_gpu_in, 0, size);
    cudaMemset(point_weights_gpu_in, 1, size);
    cudaMemset(partial_density_gpu_in, 1, size);
    cudaMemset(accumulated_density_gpu, 1, size);

    // Launch the CUDA Kernel
    int numElements = number_of_points;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    uint block_height = 1;

    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

    //////////////////////////////////////////////
    // Call the CUDA KERNEL for accumulate point
    gpu_accumulate_point_for_libxc<double,true,true,false><<<blocksPerGrid, threadsPerBlock>>> (point_weights_gpu_in,
	number_of_points, block_height,
	partial_density_gpu_in, dxyz_gpu_in.data, dd1_gpu_in.data, dd2_gpu_in.data,
	accumulated_density_gpu, dxyz_gpu_accum.data, dd1_gpu_accum.data, dd2_gpu_accum.data);

    // Create the libxcproxy
    const int nspin = 1;
    const int functionalExchange = 101;
    const int functionalCorrelation = 130;
    LibxcProxy<double,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    //////////////////////////////////////////////////
    // Calculate exc_corr and y2a in CPU with LIBXC
    libxc_exchange_correlation_cpu<double, true, true, false> (&libxcProxy,
	energy_gpu_in,
	factor_gpu_in,
	number_of_points,
	accumulated_density_gpu,
	dxyz_gpu_accum.data,
        dd1_gpu_accum.data,
	dd2_gpu_accum.data);

    ////////////////////////////////
    // Join the results for LIO
    gpu_accumulate_energy_and_forces_from_libxc<double, true, true, false><<<blocksPerGrid, threadsPerBlock>>> (
	energy_gpu_in,
	factor_gpu_in,
	point_weights_gpu_in,
	number_of_points,
	accumulated_density_gpu);

    ///////////////////////////////////////
    // Check and print the results.
    // Copy back the results before print.
    // Allocate the host input vectors
    double *energy_cpu = (double *)malloc(size);
    double *factor_cpu = (double *)malloc(size);
    double *point_weights_cpu = (double *)malloc(size);
    double *partial_density_cpu = (double *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu, energy_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu, factor_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(point_weights_cpu, point_weights_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector point_weights_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(partial_density_cpu, partial_density_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector partial_density_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    // Before print the results, get the data from cuda
    G2G::vec_type<double,4>* dxyz_cpu;
    G2G::vec_type<double,4>* dd1_cpu;
    G2G::vec_type<double,4>* dd2_cpu;

    // Alloc memory in the host for the gpu data
    uint cpu_vecs_size = number_of_points * sizeof(G2G::vec_type<double,4>);
    dxyz_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);
    dd1_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);
    dd2_cpu = (G2G::vec_type<double,4> *)malloc(cpu_vecs_size);

    // Copy data from device to host.
    err = cudaMemcpy(dxyz_cpu, dxyz_gpu_in.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dxyz_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd1_cpu, dd1_gpu_in.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd1_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(dd2_cpu, dd2_gpu_in.data, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector dd2_gpu from device to host!\n");
        exit(EXIT_FAILURE);
    }

    print_accumulate_point_data (dxyz_cpu, dd1_cpu, dd2_cpu, energy_cpu, 
	factor_cpu, point_weights_cpu, 
	partial_density_cpu, number_of_points);


    ////////////////////////////
    // Free Memory
    free(energy_cpu);
    free(factor_cpu);
    free(point_weights_cpu);
    free(partial_density_cpu);

    free(dxyz_cpu);
    free(dd1_cpu);
    free(dd2_cpu);

    cudaFree(energy_gpu_in);
    cudaFree(factor_gpu_in);
    cudaFree(point_weights_gpu_in);
    cudaFree(partial_density_gpu_in);
    cudaFree(accumulated_density_gpu);
#endif
}

///////////////////////////////////////////////////
// Cuda Matrix Helpers
//
G2G::HostMatrix< G2G::vec_type<double,4> > createMatrix(int size)
{
    printf("matrix_test0009()\n");
    //typedef G2G::vec_type<float,4> vec_type4;
    //G2G::CudaMatrix<T> function_values;
    //G2G::CudaMatrix<vec_type4> gradient_values;
    //G2G::CudaMatrix<vec_type4> hessian_values_transposed;

    //G2G::CudaMatrix< G2G::vec_type<double,4> > aCudaMatrix;
    G2G::HostMatrix< G2G::vec_type<double,4> > aHostMatrix(size, 1);

    //G2G::vec_type<float,4> one(1,1,1,1);
    G2G::vec_type<float,4> zero(0,0,0,0);

    for (int i=0; i<5; i++){
	aHostMatrix(i).x = zero.x;
	aHostMatrix(i).y = zero.y;
	aHostMatrix(i).z = zero.z;
    }

    //aCudaMatrix = aHostMatrix;
    return aHostMatrix;
}


/////////////////////////////////////////////////////
// Test: accumulate_data_for_libxc_test0007
//
void accumulate_data_for_libxc_test0007() {
    printf("accumulate_data_for_libxc_test0007()\n");
#if FULL_DOUBLE
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in = createMatrix(number_of_points);
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in = createMatrix(number_of_points);
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in = createMatrix(number_of_points);

    // Accumulate
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_accum = createMatrix(number_of_points);
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_accum = createMatrix(number_of_points);
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_accum = createMatrix(number_of_points);

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);

    dxyz_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_accum.resize(COALESCED_DIMENSION(number_of_points),m);

    // Now the arrays for energy, factors, point_weight and partial_density
    double *point_weights_gpu_in = NULL;
    double *partial_density_gpu_in = NULL;
    // Accum
    double *partial_density_gpu_accum = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);
    err = cudaMalloc((void**)&point_weights_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_accum, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_accum!\n");
    }

    // Now the arrays for energy, factors
    double *energy_gpu_in = NULL;
    double *factor_gpu_in = NULL;
    err = cudaMalloc((void**)&energy_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    ///////////////////////
    // Set data
    double partial_densities_cpu[10] = {1.016692e-33,2.333626e-34,8.367814e-34,6.744978e-35,4.493371e-36,4.396106e-37,1.908333e-34,4.848228e-35,7.228556e-34,1.717567e-38};
    double point_weights_cpu[10] = {0.000000e+00,0.000000e+00,6.356219e-06,3.324887e-04,3.143648e-02,3.212402e-01,1.299464e-05,7.277725e-04,0.000000e+00,2.066700e+00};
    double input[number_of_points];
    for (int i=0; i<number_of_points; i++) {
	input[i]=0.001*i;
    }

    cudaMemcpy(point_weights_gpu_in, point_weights_cpu, size, cudaMemcpyHostToDevice);
    cudaMemcpy(partial_density_gpu_in, partial_densities_cpu, size, cudaMemcpyHostToDevice);

    // Create the libxcproxy
    const int nspin = 1;
    const int functionalExchange = 101;
    const int functionalCorrelation = 130;
    LibxcProxy<double,4> libxcProxy(functionalExchange, functionalCorrelation, nspin);

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    uint block_height = 1;

    /////////////////////////////////
    // LIBXC VERSION
    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    gpu_accumulate_point_for_libxc<double,true,true,false><<<blocksPerGrid,threadsPerBlock>>> (point_weights_gpu_in,
	number_of_points, block_height,
	partial_density_gpu_in, dxyz_gpu_in.data, dd1_gpu_in.data, dd2_gpu_in.data,
	partial_density_gpu_accum, dxyz_gpu_accum.data, dd1_gpu_accum.data, dd2_gpu_accum.data);

    // Calculate exc_corr and y2a
    libxc_exchange_correlation_cpu<double, true, true, false> (&libxcProxy,
	energy_gpu_in,
	factor_gpu_in,
	number_of_points,
	partial_density_gpu_accum,
	dxyz_gpu_accum.data,
        dd1_gpu_accum.data,
	dd2_gpu_accum.data);

    // Join the results.
    gpu_accumulate_energy_and_forces_from_libxc<double, true, true, false><<<blocksPerGrid, threadsPerBlock>>> (
	energy_gpu_in,
	factor_gpu_in,
	point_weights_gpu_in,
	number_of_points,
	partial_density_gpu_accum);

    ///////////////////////////
    // Print libxc results
    // Allocate the host input vectors
    double *energy_cpu = (double *)malloc(size);
    double *factor_cpu = (double *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu, energy_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu_in from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu, factor_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu_in from device to host!\n");
        exit(EXIT_FAILURE);
    }

    print_accumulate_point_data (NULL, NULL, NULL, energy_cpu,
	factor_cpu, NULL, NULL, number_of_points);

    ////////////////////////////////////////
    // LIO VERSION
    // Now the arrays for energy, factors
    double *energy_gpu_in2 = NULL;
    double *factor_gpu_in2 = NULL;
    err = cudaMalloc((void**)&energy_gpu_in2, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu_in2!\n");
    }

    err = cudaMalloc((void**)&factor_gpu_in2, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu_in2!\n");
    }

    gpu_accumulate_point<double, true, true, false><<<blocksPerGrid,threadsPerBlock>>> (
	energy_gpu_in2, 
	factor_gpu_in2,
	point_weights_gpu_in,
	number_of_points,
	block_height,
	partial_density_gpu_in,
	dxyz_gpu_in.data,
	dd1_gpu_in.data,
	dd2_gpu_in.data);

    ///////////////////////////////////////
    // Print LIO results.
    // Copy back the results before print.
    // Allocate the host input vectors
    double *energy_cpu2 = (double *)malloc(size);
    double *factor_cpu2 = (double *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu2, energy_gpu_in2, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu_in2 from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu2, factor_gpu_in2, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu_in2 from device to host!\n");
        exit(EXIT_FAILURE);
    }

    print_accumulate_point_data (NULL, NULL, NULL, energy_cpu2,
	factor_cpu2, NULL, NULL, number_of_points);

    ////////////////////////////
    // Free Memory CPU
    free(energy_cpu);
    free(factor_cpu);
    free(energy_cpu2);
    free(factor_cpu2);

    ///////////////////////////
    // Free memory GPU
    cudaFree (point_weights_gpu_in);
    cudaFree (partial_density_gpu_in);
    cudaFree (partial_density_gpu_accum);
    cudaFree (energy_gpu_in);
    cudaFree (factor_gpu_in);
    cudaFree (energy_gpu_in2);
    cudaFree (factor_gpu_in2);

#endif
}



//////////////////////////////////////////////////
// Test 08
// accumulate_data_for_libxc_test0008
//
void accumulate_data_for_libxc_test0008() {
    printf("accumulate_data_for_libxc_test0008()\n");
    cudaError_t err = cudaSuccess;
    uint n = 5;
    uint m = 5;
    uint number_of_points = n+m;

    // Input
    G2G::CudaMatrix< G2G::vec_type<double,4> > dxyz_gpu_in = createMatrix(number_of_points);
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd1_gpu_in = createMatrix(number_of_points);
    G2G::CudaMatrix< G2G::vec_type<double,4> > dd2_gpu_in = createMatrix(number_of_points);

    dxyz_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd1_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);
    dd2_gpu_in.resize(COALESCED_DIMENSION(number_of_points),m);

    // Now the arrays for energy, factors, point_weight and partial_density
    double *point_weights_gpu_in = NULL;
    double *partial_density_gpu_in = NULL;

    // Create the arrays in CUDA memory.
    uint size = number_of_points * sizeof(double);
    err = cudaMalloc((void**)&point_weights_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector point_weights_gpu_in!\n");
    }

    err = cudaMalloc((void**)&partial_density_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector partial_density_gpu_in!\n");
    }

    // Now the arrays for energy, factors
    double *energy_gpu_in = NULL;
    double *factor_gpu_in = NULL;
    err = cudaMalloc((void**)&energy_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector energy_gpu!\n");
    }

    err = cudaMalloc((void**)&factor_gpu_in, size);
    if (err != cudaSuccess)
    {
	printf("Failed to allocate vector factor_gpu!\n");
    }

    ///////////////////////
    // Set random data
    double input[number_of_points];
    for (int i=0; i<number_of_points; i++) {
	input[i]=0.001*i;
    }
    cudaMemcpy(energy_gpu_in, input, size, cudaMemcpyHostToDevice);
    cudaMemcpy(factor_gpu_in, input, size, cudaMemcpyHostToDevice);
    cudaMemcpy(point_weights_gpu_in, input, size, cudaMemcpyHostToDevice);
    cudaMemcpy(partial_density_gpu_in, input, size, cudaMemcpyHostToDevice);

    ////////////////////////////////////////
    // LIO VERSION

    // Launch the CUDA Kernel
    int numElements = n+m;
    int threadsPerBlock = 32;
    int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
    uint block_height = 1;

    gpu_accumulate_point<double, true, true, false><<<blocksPerGrid,threadsPerBlock>>> (
	energy_gpu_in, 
	factor_gpu_in,
	point_weights_gpu_in,
	number_of_points,
	block_height,
	partial_density_gpu_in,
	dxyz_gpu_in.data,
	dd1_gpu_in.data,
	dd2_gpu_in.data);

    ///////////////////////////////////////
    // Print LIO results.
    // Copy back the results before print.
    // Allocate the host input vectors
    double *energy_cpu = (double *)malloc(size);
    double *factor_cpu = (double *)malloc(size);

    // Copy the vectors from gpu to cpu
    // Be aware that energy_gpu can be NULL.
    err = cudaMemcpy(energy_cpu, energy_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector energy_gpu_in from device to host!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(factor_cpu, factor_gpu_in, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        printf("Failed to copy vector factor_gpu_in from device to host!\n");
        exit(EXIT_FAILURE);
    }

    print_accumulate_point_data (NULL, NULL, NULL, energy_cpu,
	factor_cpu, NULL, NULL, number_of_points);

    ////////////////////////////
    // Free Memory CPU
    free(energy_cpu);
    free(factor_cpu);

    ///////////////////////////
    // Free memory GPU
    cudaFree (point_weights_gpu_in);
    cudaFree (partial_density_gpu_in);
    cudaFree (energy_gpu_in);
    cudaFree (factor_gpu_in);

}


/////////////////////////////////////
//// MAIN

int main(int argc, char **argv)
{
    printf("****************************\n");
    printf("** Accumulate Point CPU test  **\n");
    printf("****************************\n");

    try {
	gpu_accumulate_point_test0001();
	accumulate_data_for_libxc_test0001();
	accumulate_data_for_libxc_test0002();
	accumulate_data_for_libxc_test0003();
	accumulate_data_for_libxc_test0004();
	accumulate_data_for_libxc_test0005();
	accumulate_data_for_libxc_test0006();
	accumulate_data_for_libxc_test0007();
	accumulate_data_for_libxc_test0008();
    } catch (int e) {
	printf("An exception occurred. Exception Nr. %u \n");
	exit (EXIT_FAILURE);
    }

    printf("*************************\n");
    printf("**      Test End       **\n");
    printf("*************************\n");

    return 0;
}