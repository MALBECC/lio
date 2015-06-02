#include <fstream>
#include <iostream>
#include <vector>

#include "../../common.h"
#include "../../cuda/cuda_extra.h"
#include "../../init.h"
#include "../../matrix.h"
#include "../../timer.h"
#include "../../scalar_vector_types.h"

#include "../qmmm_integral.h"
#include "../aint_init.h"
#include "../aint_common.h"

namespace G2G
{
#include "gpu_vars/g2g_gpu_variables.h"
}

using std::cout;
using std::vector;
using std::endl;

namespace AINT
{

#include "gpu_vars/os_gpu_variables.h"
#include "gpu_vars/qmmm_gpu_variables.h"

#include "kernels/os_util.h"
#include "kernels/qmmm_forces.h"
#include "kernels/qmmm_energy.h"

//
// Send the MM atom positions and charges to the device
//
template <class scalar_type>
bool QMMMIntegral<scalar_type>::load_clatoms( void )
{

    cudaMemcpyToSymbol(gpu_clatoms, &integral_vars.clatoms, sizeof(integral_vars.clatoms), 0, cudaMemcpyHostToDevice);

    G2G::HostMatrix<G2G::vec_type<scalar_type,3> > clatom_pos_cpu(integral_vars.clatoms,1);
    G2G::HostMatrix<scalar_type> clatom_chg_cpu(integral_vars.clatoms,1);
    for (uint i = 0; i < integral_vars.clatoms; i++) {
      clatom_pos_cpu(i) = G2G::vec_type<scalar_type,3>(integral_vars.clatom_positions(i).x,integral_vars.clatom_positions(i).y,integral_vars.clatom_positions(i).z);
      clatom_chg_cpu(i) = integral_vars.clatom_charges(i);
    }
    size_t gpu_size = clatom_pos_cpu.bytes() + clatom_chg_cpu.bytes();
    float mb_size = (float)gpu_size / 1048576.0f;
    cout << "QM/MM input size: " << mb_size << " MB" << endl;

    clatom_pos_dev = clatom_pos_cpu;
    clatom_chg_dev = clatom_chg_cpu;

    cudaAssertNoError("QMMMIntegral::load_clatoms");

    return true;
}

//
// Allocate output arrays on device
//
template <class scalar_type>
bool QMMMIntegral<scalar_type>::alloc_output( void )
{
    // Currently, each block in the kernel reduces its output, so the output arrays have length (# block)
    uint partial_out_size = 0;
    // Output arrays (probably) don't need to be padded for alignment as only one (or three) threads per block write to them
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      uint this_count = divUp(os_int.term_type_counts[i],QMMM_BLOCK_SIZE);
      partial_out_size += this_count;
    }
    size_t gpu_size = COALESCED_DIMENSION(partial_out_size) * integral_vars.clatoms * 3 * sizeof(scalar_type);
    float mb_size = (float)gpu_size / 1048576.0f;
    cout << "QM/MM output size: " << mb_size << " MB" << endl;

    // Forces: output is partial forces
    partial_mm_forces_dev.resize(COALESCED_DIMENSION(partial_out_size), integral_vars.clatoms);

    cudaAssertNoError("QMMMIntegral::alloc_output");

    return true;
}

template <class scalar_type>
void QMMMIntegral<scalar_type>::clear( void )
{
    clatom_pos_dev.deallocate();
    clatom_chg_dev.deallocate();

    partial_mm_forces_dev.deallocate();

    cudaAssertNoError("QMMMIntegral::clear");
}

template <class scalar_type>
void QMMMIntegral<scalar_type>::calc_fock( double& Es )
{
    os_int.reload_density();

    uint partial_out_size = 0, max_partial_size = 0;
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      uint this_count = divUp(os_int.term_type_counts[i],QMMM_BLOCK_SIZE);
      if (this_count > max_partial_size) max_partial_size = this_count;
      partial_out_size += this_count;
    }

    {
      dim3 threads(os_int.dens_values.size());
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads,blockSize);
      // Zero the partial Fock matrix on the GPU
      zero_fock<double><<<gridSize,blockSize>>>(os_int.partial_fock_dev.data,os_int.dens_values.size(),1);;
    }

    //
    // The STR table for F(m,U) calculation is being accessed via texture fetches
    //
    cudaBindTextureToArray(str_tex,gammaArray);

#define qmmm_fock_parameters \
  os_int.term_type_counts[i], os_int.factor_ac_dev.data, os_int.nuc_dev.data, os_int.func_code_dev.data+offset,os_int.local_dens_dev.data+offset, \
  os_int.partial_fock_dev.data+fock_offset, os_int.dens_values.size(),clatom_pos_dev.data,clatom_chg_dev.data
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_TERM_TYPES];
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (uint i = 0; i < NUM_TERM_TYPES; i++)
    {
      uint offset = os_int.term_type_offsets[i];
      uint fock_offset = os_int.dens_offsets[i];
      dim3 threads = os_int.term_type_counts[i];
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads, blockSize);
      switch (i) {
        case 0: gpu_qmmm_fock<scalar_type,0><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 1: gpu_qmmm_fock<scalar_type,1><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 2: gpu_qmmm_fock<scalar_type,2><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 3: gpu_qmmm_fock<scalar_type,3><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 4: gpu_qmmm_fock<scalar_type,4><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
        case 5: gpu_qmmm_fock<scalar_type,5><<<gridSize,blockSize,0,stream[i]>>>( qmmm_fock_parameters ); break;
      }

      //
      // Reduce the partial Fock terms for a particular term type as soon as that kernel is done; also calculate partial energies for that type
      //
      dim3 reduceThreads = os_int.dens_counts[i];
      dim3 reduceBlockSize(QMMM_REDUCE_BLOCK_SIZE);
      dim3 reduceGridSize = divUp(reduceThreads,reduceBlockSize);
      gpu_fock_reduce<scalar_type><<<reduceGridSize,reduceBlockSize,0,stream[i]>>>( os_int.partial_fock_dev.data+fock_offset,os_int.dens_values_dev.data+fock_offset,
                                                                                    os_int.partial_energies_dev.data+os_int.energies_offsets[i],
                                                                                    os_int.dens_values.size(),max_partial_size,os_int.dens_counts[i] );
    }
    cudaDeviceSynchronize();
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }

    cudaUnbindTexture(str_tex);

    os_int.get_fock_output(Es,integral_vars.rmm_1e_output);

    cudaAssertNoError("QMMMIntegral::calc_fock");
}

template<class scalar_type>
void QMMMIntegral<scalar_type>::calc_gradient(double* qm_forces, double* mm_forces, bool do_cl, bool do_qm)
{
    os_int.reload_density();

    uint partial_out_size = 0, max_partial_size = 0;
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      uint this_count = divUp(os_int.term_type_counts[i],QMMM_BLOCK_SIZE);
      if (this_count > max_partial_size) max_partial_size = this_count;
      partial_out_size += this_count;
    }
    {
        dim3 threads(COALESCED_DIMENSION(partial_out_size),G2G::fortran_vars.atoms);
        dim3 blockSize(32,4);
        dim3 gridSize = divUp(threads,blockSize);
        zero_forces<scalar_type><<<gridSize,blockSize>>>(os_int.partial_qm_forces_dev.data,COALESCED_DIMENSION(partial_out_size),G2G::fortran_vars.atoms);
    }

    //
    // The STR table for F(m,U) calculation is being accessed via texture fetches
    //
    cudaBindTextureToArray(str_tex,gammaArray);

#define qmmm_forces_parameters \
  os_int.term_type_counts[i], os_int.factor_ac_dev.data, os_int.nuc_dev.data, os_int.dens_values_dev.data+dens_offset, os_int.func_code_dev.data+offset,os_int.local_dens_dev.data+offset, \
  partial_mm_forces_dev.data+force_offset, os_int.partial_qm_forces_dev.data+force_offset, COALESCED_DIMENSION(partial_out_size),clatom_pos_dev.data,clatom_chg_dev.data
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_TERM_TYPES];
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (uint i = 0; i < NUM_TERM_TYPES; i++)
    {
      uint offset = os_int.term_type_offsets[i];
      uint dens_offset = os_int.dens_offsets[i];
      uint force_offset = os_int.out_offsets[i];
      dim3 threads = os_int.term_type_counts[i];
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads, blockSize);
      if (do_cl) {
        if (do_qm) {
          switch (i) {
            case 0: gpu_qmmm_forces<scalar_type,0,true,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 1: gpu_qmmm_forces<scalar_type,1,true,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 2: gpu_qmmm_forces<scalar_type,2,true,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 3: gpu_qmmm_forces<scalar_type,3,true,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 4: gpu_qmmm_forces<scalar_type,4,true,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 5: gpu_qmmm_forces<scalar_type,5,true,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
          }
        } else {
          switch (i) {
            case 0: gpu_qmmm_forces<scalar_type,0,true,false><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 1: gpu_qmmm_forces<scalar_type,1,true,false><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 2: gpu_qmmm_forces<scalar_type,2,true,false><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 3: gpu_qmmm_forces<scalar_type,3,true,false><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 4: gpu_qmmm_forces<scalar_type,4,true,false><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 5: gpu_qmmm_forces<scalar_type,5,true,false><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
          }
        }
      } else {
        if (do_qm) {
          switch (i) {
            case 0: gpu_qmmm_forces<scalar_type,0,false,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 1: gpu_qmmm_forces<scalar_type,1,false,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 2: gpu_qmmm_forces<scalar_type,2,false,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 3: gpu_qmmm_forces<scalar_type,3,false,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 4: gpu_qmmm_forces<scalar_type,4,false,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
            case 5: gpu_qmmm_forces<scalar_type,5,false,true><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
          }
        } else { // This should never be launched
        }
      }
    }
    cudaDeviceSynchronize();
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }

    cudaUnbindTexture(str_tex);

    os_int.get_gradient_output(qm_forces,partial_out_size);
    if (integral_vars.clatoms > 0) get_gradient_output(mm_forces,partial_out_size);

    cudaAssertNoError("QMMMIntegral::calc_gradient");
}

template<class scalar_type>
void QMMMIntegral<scalar_type>::get_gradient_output(double* mm_forces, uint partial_out_size)
{
    G2G::HostMatrix<G2G::vec_type<scalar_type,3> > cpu_partial_mm_forces(partial_mm_forces_dev);

    //
    // Accumulate partial output
    //
    // TODO: need to think about how to accumulate individual force terms
    // Currently, we reduce on a per-block basis in the kernel, then accumulate the block results here on the host
    //
    // The energy partial results are being reduced on-device and that works very well, could probably do that for forces too
    //
    for (uint i = 0; i < integral_vars.clatoms; i++) {
      for (uint j = 0; j < partial_out_size; j++) {
        mm_forces[i + 0 * integral_vars.clatoms] += cpu_partial_mm_forces(j,i).x;
        mm_forces[i + 1 * integral_vars.clatoms] += cpu_partial_mm_forces(j,i).y;
        mm_forces[i + 2 * integral_vars.clatoms] += cpu_partial_mm_forces(j,i).z;
      }
    }

    cudaAssertNoError("QMMMIntegral::get_gradient_output");
}

#if AINT_MP && !FULL_DOUBLE
template class QMMMIntegral<float>;
#else
template class QMMMIntegral<double>;
#endif

}
