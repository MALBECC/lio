/* -*- mode: c -*- */
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math_constants.h>
#include <string>
#include <vector>
#include <stdio.h>
#include "../../init.h"
#include "../../common.h"
#include "../../partition.h"
#include "../../timer.h"

#if USE_LIBXC
 #include "../../libxc/libxc_accumulate_point.h"
#endif

#define POINTS_BLOCK_SIZE 256
#include "accumulate_values.h"

#include "gpu_calc_gradients.h"
#include "gpu_partial_forces.h"

using namespace std;

namespace G2G {
#if FULL_DOUBLE
texture<int2, 2, cudaReadModeElementType> rmm_gpu_for_tex;
texture<int2, 2, cudaReadModeElementType> tred_gpu_for_tex;
texture<int2, 2, cudaReadModeElementType> diff_gpu_for_tex;
#else
texture<float, 2, cudaReadModeElementType> rmm_gpu_for_tex;
texture<float, 2, cudaReadModeElementType> tred_gpu_for_tex;
texture<float, 2, cudaReadModeElementType> diff_gpu_for_tex;
#endif

#include "../../cuda/kernels/transpose.h"
#include "ES_compute_for_partial.h"

template<class scalar_type> void PointGroupGPU<scalar_type>::
        solve_for_exc(double*P,double*V,HostMatrix<double>& F,int MET)
{
   uint group_m = this->total_functions();
   bool lda = false;
   bool compute_forces = false;
   compute_functions(compute_forces, !lda);

// Point weights on CPU and GPU
   CudaMatrix<scalar_type> point_weights_gpu;
   HostMatrix<scalar_type> point_weights_cpu(this->number_of_points, 1);
   uint i = 0;
   for (vector<Point>::const_iterator p = this->points.begin(); p != this->points.end(); ++p, ++i) {
     point_weights_cpu(i) = p->weight;
   }
   point_weights_gpu = point_weights_cpu;
   point_weights_cpu.deallocate();

// Variables to kernels: gpu_compute
   dim3 threadBlock, threadGrid;
   const int block_height= divUp(group_m, 2*DENSITY_BLOCK_SIZE);
   threadBlock = dim3(DENSITY_BLOCK_SIZE,1,1); // Hay que asegurarse que la cantidad de funciones este en rango
   threadGrid = dim3(this->number_of_points,block_height,1);

// Partial Transition Density GPU
   CudaMatrix<scalar_type> partial_tred_gpu;
   CudaMatrix< vec_type<scalar_type,4> > tredxyz_gpu;
   partial_tred_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
   tredxyz_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height);

// Partial Difference Density GPU
   CudaMatrix<scalar_type> partial_diff_gpu;
   CudaMatrix< vec_type<scalar_type,4> > diffxyz_gpu;
   partial_diff_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
   diffxyz_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height);

// Accumulate Values
   CudaMatrix<scalar_type> tred_accum_gpu;
   CudaMatrix< vec_type<scalar_type,4> > tredxyz_accum_gpu;
   tred_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
   tredxyz_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));

   CudaMatrix<scalar_type> diff_accum_gpu;
   CudaMatrix< vec_type<scalar_type,4> > diffxyz_accum_gpu;
   diff_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
   diffxyz_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));

// Variables to kernels: accumulate density
   const dim3 threadGrid_accumulate(divUp(this->number_of_points,DENSITY_ACCUM_BLOCK_SIZE),1,1);
   const dim3 threadBlock_accumulate(DENSITY_ACCUM_BLOCK_SIZE,1,1);

// Transpose functions and gradients
   #define BLOCK_DIM 16
   int transposed_width = COALESCED_DIMENSION(this->number_of_points);
   dim3 transpose_grid(transposed_width / BLOCK_DIM, divUp((group_m),BLOCK_DIM), 1);
   dim3 transpose_threads(BLOCK_DIM, BLOCK_DIM, 1);

   CudaMatrix<scalar_type> function_values_transposed;
   CudaMatrix<vec_type<scalar_type,4> > gradient_values_transposed;
   function_values_transposed.resize(group_m, COALESCED_DIMENSION(this->number_of_points));
   gradient_values_transposed.resize( group_m,COALESCED_DIMENSION(this->number_of_points));

   transpose<<<transpose_grid, transpose_threads>>> (function_values_transposed.data,
       function_values.data, COALESCED_DIMENSION(this->number_of_points), group_m);

   transpose<<<transpose_grid, transpose_threads>>> (gradient_values_transposed.data,
       gradient_values.data, COALESCED_DIMENSION(this->number_of_points), group_m );

   int M = fortran_vars.m;
// FORM reduce Transition density
   HostMatrix<scalar_type> tred_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);

// FORM reduce Difference density
   HostMatrix<scalar_type> diff_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);

   HostMatrix<double> Pbig(M*(M+1)/2);
   HostMatrix<double> Vbig(M*(M+1)/2);
   int index = 0;
   int row, col;
   for(row=0;row<M;row++) {
     Pbig(index) = P[row*M+row];
     Vbig(index) = V[row*M+row];
     index += 1;
     for(col=row+1;col<M;col++) {
        Pbig(index) = P[row*M+col] + P[col*M+row];
        Vbig(index) = V[row*M+col] + V[col*M+row];
        index += 1;
     }
   }
   get_tred_input(tred_cpu,Vbig);
   get_tred_input(diff_cpu,Pbig);

// We put zero
   for (uint i=0; i<(group_m+DENSITY_BLOCK_SIZE); i++)
   {
     for(uint j=0; j<COALESCED_DIMENSION(group_m); j++)
     {
       if((i>=group_m) || (j>=group_m) || (j > i))
       {
         tred_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
         diff_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
       }
     }
   }

// Form Bind Texture
   cudaArray* cuArraytred;
   cudaMallocArray(&cuArraytred, &tred_gpu_for_tex.channelDesc, tred_cpu.width, tred_cpu.height);
   cudaMemcpyToArray(cuArraytred,0,0,tred_cpu.data,sizeof(scalar_type)*tred_cpu.width*tred_cpu.height,cudaMemcpyHostToDevice);
   cudaBindTextureToArray(tred_gpu_for_tex, cuArraytred);
   cudaArray* cuArraydiff;
   cudaMallocArray(&cuArraydiff, &diff_gpu_for_tex.channelDesc, diff_cpu.width, diff_cpu.height);
   cudaMemcpyToArray(cuArraydiff,0,0,diff_cpu.data,sizeof(scalar_type)*diff_cpu.width*diff_cpu.height,cudaMemcpyHostToDevice);
   cudaBindTextureToArray(diff_gpu_for_tex, cuArraydiff);

   tred_cpu.deallocate(); diff_cpu.deallocate();

   // If we saved the GS density and derivatives
   if (fortran_vars.den_point_save == 0) {
      // CALCULATE PARTIAL DENSITIES
      #define compden_parameter \
         point_weights_gpu.data,this->number_of_points,function_values_transposed.data,\
         group_m,gradient_values_transposed.data, partial_tred_gpu.data,tredxyz_gpu.data, \
         partial_diff_gpu.data, diffxyz_gpu.data
      ES_compute_for_partial<scalar_type,true,true,false><<<threadGrid, threadBlock>>>(compden_parameter);

      // ACCUMULATE DENSITIES
      #define accumden_parameter \
         this->number_of_points, block_height, partial_tred_gpu.data, \
         partial_diff_gpu.data, tredxyz_gpu.data, diffxyz_gpu.data, \
         tred_accum_gpu.data, diff_accum_gpu.data, tredxyz_accum_gpu.data, diffxyz_accum_gpu.data
      accumulate_values<scalar_type,true,true,false><<<threadGrid_accumulate,threadBlock_accumulate>>>(accumden_parameter);
      #undef compute_parameters
      #undef accumulate_parameters

   // When we did not save the GS density and derivatives
   }  else {
      // Ground State Density GPU
      CudaMatrix<scalar_type> partial_densities_gpu;
      CudaMatrix< vec_type<scalar_type,4> > dxyz_gpu;
      partial_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
      dxyz_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height);

      // Accumulate Values
      rmm_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
      dxyz_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));

      // FORM reduce GS density
      HostMatrix<scalar_type> rmm_cpuX(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
      get_rmm_input(rmm_cpuX);

      // ponemos ceros fuera de group_m
      for (uint i=0; i<(group_m+DENSITY_BLOCK_SIZE); i++)
      {
        for(uint j=0; j<COALESCED_DIMENSION(group_m); j++)
        {
          if((i>=group_m) || (j>=group_m) || (j > i))
          {
            rmm_cpuX.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
          }
        }
      }
      // Form Bind Textures
      cudaArray* cuArrayrmm;
      cudaMallocArray(&cuArrayrmm, &rmm_gpu_for_tex.channelDesc, rmm_cpuX.width, rmm_cpuX.height);
      cudaMemcpyToArray(cuArrayrmm,0,0,rmm_cpuX.data,sizeof(scalar_type)*rmm_cpuX.width*rmm_cpuX.height,cudaMemcpyHostToDevice);
      cudaBindTextureToArray(rmm_gpu_for_tex, cuArrayrmm);
      rmm_cpuX.deallocate();
      // CALCULATE PARTIAL DENSITIES
      #define compden_parameter \
         point_weights_gpu.data,this->number_of_points,function_values_transposed.data,\
         group_m,gradient_values_transposed.data, partial_tred_gpu.data,tredxyz_gpu.data, \
         partial_diff_gpu.data, diffxyz_gpu.data, partial_densities_gpu.data, dxyz_gpu.data
      ES_compute_for_partial<scalar_type,true,true,false><<<threadGrid, threadBlock>>>(compden_parameter);

      // ACCUMULATE DENSITIES
      #define accumden_parameter \
         this->number_of_points, block_height, \
         partial_densities_gpu.data, partial_tred_gpu.data, partial_diff_gpu.data, \
         dxyz_gpu.data, tredxyz_gpu.data, diffxyz_gpu.data, \
         rmm_accum_gpu.data,tred_accum_gpu.data, diff_accum_gpu.data, \
         dxyz_accum_gpu.data,tredxyz_accum_gpu.data, diffxyz_accum_gpu.data
      accumulate_values<scalar_type,true,true,false><<<threadGrid_accumulate,threadBlock_accumulate>>>(accumden_parameter);
      #undef compute_parameters
      #undef accumulate_parameters

      cudaUnbindTexture(rmm_gpu_for_tex);
      cudaFreeArray(cuArrayrmm);
      partial_densities_gpu.deallocate();
      dxyz_gpu.deallocate();
   } // end density save

   partial_tred_gpu.deallocate(); partial_diff_gpu.deallocate();
   tredxyz_gpu.deallocate(); diffxyz_gpu.deallocate();

// OBTAIN GRADIENTS TERMS
   CudaMatrix<scalar_type> gdens, ddens, tdens;
   CudaMatrix< vec_type<scalar_type,4> > gdens_xyz, ddens_xyz, tdens_xyz;
   gdens.resize(COALESCED_DIMENSION(this->number_of_points));
   ddens.resize(COALESCED_DIMENSION(this->number_of_points));
   tdens.resize(COALESCED_DIMENSION(this->number_of_points)); tdens.zero();
   gdens_xyz.resize(COALESCED_DIMENSION(this->number_of_points));
   ddens_xyz.resize(COALESCED_DIMENSION(this->number_of_points));
   tdens_xyz.resize(COALESCED_DIMENSION(this->number_of_points)); tdens_xyz.zero();

   gpu_calc_gradients<scalar_type,true,true,false>(this->number_of_points,rmm_accum_gpu.data,tred_accum_gpu.data,diff_accum_gpu.data,
                      dxyz_accum_gpu.data, tredxyz_accum_gpu.data, diffxyz_accum_gpu.data,
                      gdens.data, tdens.data, ddens.data, gdens_xyz.data, tdens_xyz.data, ddens_xyz.data, MET);


// FORCES CALCULATION
   CudaMatrix<scalar_type> mat_dens_gpu(group_m,group_m);
   CudaMatrix<scalar_type> mat_diff_gpu(group_m,group_m);
   CudaMatrix<scalar_type> mat_tred_gpu(group_m,group_m);
   HostMatrix<scalar_type> rmm_cpu;
   rmm_cpu.resize(group_m,group_m);  get_rmm_input(rmm_cpu);
   diff_cpu.resize(group_m,group_m); get_tred_input(diff_cpu,Pbig);
   tred_cpu.resize(group_m,group_m); get_tred_input(tred_cpu,Vbig);
   Pbig.deallocate(); Vbig.deallocate();

   // Copy Matrices
   int dim_mat = group_m*group_m;
   cudaMemcpy(mat_dens_gpu.data,rmm_cpu.data,sizeof(scalar_type)*dim_mat,cudaMemcpyHostToDevice);
   cudaMemcpy(mat_tred_gpu.data,tred_cpu.data,sizeof(scalar_type)*dim_mat,cudaMemcpyHostToDevice);
   cudaMemcpy(mat_diff_gpu.data,diff_cpu.data,sizeof(scalar_type)*dim_mat,cudaMemcpyHostToDevice);
   rmm_cpu.deallocate(); diff_cpu.deallocate(); tred_cpu.deallocate();

   int block_for = divUp(this->number_of_points, POINTS_BLOCK_SIZE);
   CudaMatrix< vec_type<scalar_type,4> > forces_basis;
   forces_basis.resize(group_m,block_for); forces_basis.zero();
   CudaMatrix< vec_type<scalar_type,4> > forces_basis_accum;
   forces_basis_accum.resize(group_m); forces_basis_accum.zero();

   threadGrid  = dim3(group_m,1,1);
   dim3 threadBlock_partial = dim3(block_for,1,1);
   dim3 threadBlock_accum = dim3(block_for,1,1);

#define partial_forces \
   function_values_transposed.data, gradient_values_transposed.data, hessian_values_transposed.data, \
   mat_dens_gpu.data, mat_tred_gpu.data, mat_diff_gpu.data, gdens.data, tdens.data, ddens.data, \
   gdens_xyz.data, tdens_xyz.data, ddens_xyz.data, forces_basis.data, group_m, this->number_of_points, \
   point_weights_gpu.data, block_for
   
   gpu_partial_forces<scalar_type,true,true,false><<<threadGrid,threadBlock_partial>>>(partial_forces);
   gpu_accum_forces<scalar_type,true,true,false><<<threadGrid,threadBlock_accum>>>(forces_basis.data,
                                                   forces_basis_accum.data, block_for, group_m);
#undef partial_forces

   //cudaThreadSynchronize();
   HostMatrix< vec_type<scalar_type,4> > forces_basis_cpu;
   forces_basis_cpu.resize(group_m); forces_basis_cpu.zero();
   cudaMemcpy(forces_basis_cpu.data,forces_basis_accum.data,group_m*sizeof(vec_type<scalar_type,4>), cudaMemcpyDeviceToHost);
   forces_basis.deallocate(); forces_basis_accum.deallocate();

   int local_atoms = this->total_nucleii(); // da cuantos atomos locales hay en el grupo
   HostMatrix<scalar_type> ddx, ddy, ddz;
   ddx.resize(local_atoms, 1); ddx.zero();
   ddy.resize(local_atoms, 1); ddy.zero();
   ddz.resize(local_atoms, 1); ddz.zero();

   for(int i=0; i<group_m; i++) {
      uint nuc = this->func2local_nuc(i);
      ddx(nuc) -= forces_basis_cpu(i).x;
      ddy(nuc) -= forces_basis_cpu(i).y;
      ddz(nuc) -= forces_basis_cpu(i).z;
   }
   forces_basis_cpu.deallocate();

   for (int i = 0; i < local_atoms; i++) {
      uint global_atom = this->local2global_nuc[i]; // da el indice del atomo LOCAL al GLOBAL
      F(global_atom,0) += ddx(i);
      F(global_atom,1) += ddy(i);
      F(global_atom,2) += ddz(i);
   }
   ddx.deallocate(); ddy.deallocate(), ddz.deallocate();
   gdens.deallocate(); tdens.deallocate(); ddens.deallocate();
   gdens_xyz.deallocate(); tdens_xyz.deallocate(); ddens_xyz.deallocate();

// Free Texture and Memory
   cudaUnbindTexture(tred_gpu_for_tex);
   cudaUnbindTexture(diff_gpu_for_tex);
   cudaFreeArray(cuArraytred);
   cudaFreeArray(cuArraydiff);
   mat_dens_gpu.deallocate(); mat_diff_gpu.deallocate(); 
   mat_tred_gpu.deallocate(); diff_accum_gpu.deallocate(); 
   tred_accum_gpu.deallocate(); diffxyz_accum_gpu.deallocate(); 
   tredxyz_accum_gpu.deallocate(); point_weights_gpu.deallocate(); 
   function_values_transposed.deallocate(); gradient_values_transposed.deallocate();
   if (fortran_vars.den_point_save != 0) {
      rmm_accum_gpu.deallocate();
      dxyz_accum_gpu.deallocate();
   }
}
#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupGPU<double>;
#else
template class PointGroup<float>;
template class PointGroupGPU<float>;
#endif
}
