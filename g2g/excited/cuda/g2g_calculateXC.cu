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

#include "accumulate_values.h"

using namespace std;

namespace G2G {
#if FULL_DOUBLE
texture<int2, 2, cudaReadModeElementType> rmm_gpu_tex;
texture<int2, 2, cudaReadModeElementType> tred_gpu_tex;
#else
texture<float, 2, cudaReadModeElementType> rmm_gpu_tex;
texture<float, 2, cudaReadModeElementType> tred_gpu_tex;
#endif

#include "../../cuda/kernels/transpose.h"
#include "obtain_fock_cuda.h"
#include "obtain_terms.h"

// comienzan las nuevas
#include "GS_compute_partial.h"
#include "ES_compute_partial.h"

template<class scalar_type>
void PointGroupGPU<scalar_type>::solve_closed_lr(double* T, HostMatrix<double>& Fock)
{
   cudaError_t err = cudaSuccess;
   Timers timers;

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

// Accumulate Values
   CudaMatrix<scalar_type> tred_accum_gpu;
   CudaMatrix< vec_type<scalar_type,4> > tredxyz_accum_gpu;
   tred_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
   tredxyz_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));

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

// FORM reduce Transition density
   HostMatrix<scalar_type> tred_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
   int M = fortran_vars.m;

   HostMatrix<double> Tbig(M*(M+1)/2);
   int index = 0;
   int row, col;
   for(row=0;row<M;row++) {
     Tbig(index) = T[row*M+row];
     index += 1;
     for(col=row+1;col<M;col++) {
        Tbig(index) = T[row*M+col] + T[col*M+row];
        index += 1;
     }
   }
   get_tred_input(tred_cpu,Tbig);
   Tbig.deallocate();

// ponemos ceros fuera de group_m
   for (uint i=0; i<(group_m+DENSITY_BLOCK_SIZE); i++)
   {
     for(uint j=0; j<COALESCED_DIMENSION(group_m); j++)
     {
       if((i>=group_m) || (j>=group_m) || (j > i))
       {
         tred_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
       }
     }
   }

// Form Bind Textures
   cudaArray* cuArraytred;
   cudaMallocArray(&cuArraytred, &tred_gpu_tex.channelDesc, tred_cpu.width, tred_cpu.height);
   cudaMemcpyToArray(cuArraytred,0,0,tred_cpu.data,sizeof(scalar_type)*tred_cpu.width*tred_cpu.height,cudaMemcpyHostToDevice);
   cudaBindTextureToArray(tred_gpu_tex, cuArraytred);
   tred_cpu.deallocate();

// CALCULATE PARTIAL DENSITIES
#define compden_parameter \
   this->number_of_points,function_values_transposed.data,group_m,gradient_values_transposed.data,\
   partial_tred_gpu.data,tredxyz_gpu.data
   ES_compute_partial<scalar_type,true,true,false><<<threadGrid, threadBlock>>>(compden_parameter);

// ACCUMULATE DENSITIES
#define accumden_parameter \
   this->number_of_points, block_height, partial_tred_gpu.data, tredxyz_gpu.data, \
   tred_accum_gpu.data, tredxyz_accum_gpu.data
   accumulate_values<scalar_type,true,true,false><<<threadGrid_accumulate,threadBlock_accumulate>>>(accumden_parameter);

#undef compute_parameters
#undef accumulate_parameters

// LIBXC INITIALIZATION
   const int nspin = XC_UNPOLARIZED;
   const int functionalExchange = fortran_vars.ex_functional_id + 1000; // 1101;
   const int functionalCorrelation = fortran_vars.ec_functional_id + 1000; // 1130;
   LibxcProxy_cuda<scalar_type,4> libxcProxy_cuda(functionalExchange, functionalCorrelation, nspin, fortran_vars.fexc);

   CudaMatrix<scalar_type> lrCoef_gpu;
   lrCoef_gpu.resize(COALESCED_DIMENSION(this->number_of_points));

// DEFINE OUTPUTS
   CudaMatrix< vec_type<scalar_type,4> > Txyz;
   CudaMatrix< vec_type<scalar_type,4> > Dxyz;
   Txyz.resize(COALESCED_DIMENSION(this->number_of_points));
   Dxyz.resize(COALESCED_DIMENSION(this->number_of_points));

   libxc_gpu_coefLR<scalar_type, true, true, false>(&libxcProxy_cuda,this->number_of_points,
               rmm_accum_gpu.data,tred_accum_gpu.data,dxyz_accum_gpu.data,
               tredxyz_accum_gpu.data,
               // Outputs
               Dxyz.data, Txyz.data, lrCoef_gpu.data);

// CALCULATE TERMS
   CudaMatrix<scalar_type> terms_lr;
   terms_lr.resize(group_m,COALESCED_DIMENSION(this->number_of_points));
   threadGrid  = dim3(this->number_of_points);
   threadBlock = dim3(group_m);

#define compute_parameters \
   this->number_of_points, group_m, point_weights_gpu.data, function_values_transposed.data, \
   gradient_values_transposed.data, lrCoef_gpu.data, Dxyz.data, \
   Txyz.data, terms_lr.data

   gpu_obtain_term<scalar_type,false,true,false><<<threadGrid,threadBlock>>>(compute_parameters);
#undef compute_parameters

// CALCULATE FOCK
   const int block_size = divUp(group_m,RMM_BLOCK_SIZE_XY);
   const int MM = block_size * (block_size+1) / 2;
   threadGrid = dim3(MM,1,1);
   threadBlock = dim3(RMM_BLOCK_SIZE_XY,RMM_BLOCK_SIZE_XY);

   CudaMatrix<scalar_type> smallFock_gpu(group_m,group_m);
   smallFock_gpu.zero();
#define compute_parameters \
   this->number_of_points, group_m, point_weights_gpu.data, function_values_transposed.data, terms_lr.data, \
   smallFock_gpu.data

   gpu_obtain_fock<scalar_type,false,true,false><<<threadGrid,threadBlock>>>(compute_parameters);
#undef compute_parameters

// Obtain the global fock
   HostMatrix<scalar_type> smallFock(smallFock_gpu);
   smallFock_gpu.deallocate();
   this->add_rmm_output(smallFock,Fock);

// Free Memory
   smallFock.deallocate();
   cudaUnbindTexture(tred_gpu_tex);
   cudaFreeArray(cuArraytred);
   Txyz.deallocate();
   Dxyz.deallocate();
   partial_tred_gpu.deallocate();
   tredxyz_gpu.deallocate();
   tred_accum_gpu.deallocate();
   tredxyz_accum_gpu.deallocate();
   lrCoef_gpu.deallocate();
}

template <class scalar_type>
void PointGroupGPU<scalar_type>::get_tred_input(
    HostMatrix<scalar_type>& tred_input, HostMatrix<double>& source) const
{
  tred_input.zero();
  const int indexes = this->rmm_bigs.size();
  for (int i = 0; i < indexes; i++) {
    int ii = this->rmm_rows[i], jj = this->rmm_cols[i], bi = this->rmm_bigs[i];
    tred_input(ii, jj) = tred_input(jj, ii) = (scalar_type)source(bi);
  }

}

template<class scalar_type> void PointGroupGPU<scalar_type>::
               lr_closed_init()
{
   cudaError_t err = cudaSuccess;
   uint group_m = this->total_functions();
   bool lda = false;
   bool compute_forces = false;

   compute_functions(compute_forces, !lda);

// Variables to kernels: gpu_compute
   dim3 threadBlock, threadGrid;
   const int block_height= divUp(group_m, 2*DENSITY_BLOCK_SIZE);
   threadBlock = dim3(DENSITY_BLOCK_SIZE,1,1); // Hay que asegurarse que la cantidad de funciones este en rango
   threadGrid = dim3(this->number_of_points,block_height,1);

// Ground State Density GPU
   CudaMatrix<scalar_type> partial_densities_gpu;
   CudaMatrix< vec_type<scalar_type,4> > dxyz_gpu;
   partial_densities_gpu.resize(COALESCED_DIMENSION(this->number_of_points), block_height);
   dxyz_gpu.resize(COALESCED_DIMENSION(this->number_of_points),block_height);

// Accumulate Values
   rmm_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));
   dxyz_accum_gpu.resize(COALESCED_DIMENSION(this->number_of_points));

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

// FORM reduce GS density
   HostMatrix<scalar_type> rmm_cpu(COALESCED_DIMENSION(group_m), group_m+DENSITY_BLOCK_SIZE);
   get_rmm_input(rmm_cpu);

// ponemos ceros fuera de group_m
   for (uint i=0; i<(group_m+DENSITY_BLOCK_SIZE); i++)
   {
     for(uint j=0; j<COALESCED_DIMENSION(group_m); j++)
     {
       if((i>=group_m) || (j>=group_m) || (j > i))
       {
         rmm_cpu.data[COALESCED_DIMENSION(group_m)*i+j]=0.0f;
       }
     }
   }

// Form Bind Textures
   cudaArray* cuArrayrmm;
   cudaMallocArray(&cuArrayrmm, &rmm_gpu_tex.channelDesc, rmm_cpu.width, rmm_cpu.height);
   cudaMemcpyToArray(cuArrayrmm,0,0,rmm_cpu.data,sizeof(scalar_type)*rmm_cpu.width*rmm_cpu.height,cudaMemcpyHostToDevice);
   cudaBindTextureToArray(rmm_gpu_tex, cuArrayrmm);
   rmm_cpu.deallocate();

// CALCULATE PARTIAL DENSITIES
#define compden_parameter \
   this->number_of_points,function_values_transposed.data,group_m,gradient_values_transposed.data, \
   partial_densities_gpu.data,dxyz_gpu.data
   GS_compute_partial<scalar_type,true,true,false><<<threadGrid, threadBlock>>>(compden_parameter);

// ACCUMULATE DENSITIES
#define accumden_parameter \
   this->number_of_points, block_height, partial_densities_gpu.data, dxyz_gpu.data, \
   rmm_accum_gpu.data, dxyz_accum_gpu.data
   accumulate_values<scalar_type,true,true,false><<<threadGrid_accumulate,threadBlock_accumulate>>>(accumden_parameter);

#undef compute_parameters
#undef accumulate_parameters

// FREE MEMORY
   cudaUnbindTexture(rmm_gpu_tex);
   cudaFreeArray(cuArrayrmm);
   partial_densities_gpu.deallocate();
   dxyz_gpu.deallocate();
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupGPU<double>;
#else
template class PointGroup<float>;
template class PointGroupGPU<float>;
#endif
}

