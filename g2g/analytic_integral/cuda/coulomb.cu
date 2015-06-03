#include <fstream>
#include <iostream>
#include <vector>

#include "../../common.h"
#include "../../init.h"
#include "../../cuda/cuda_extra.h"
#include "../../matrix.h"
#include "../../timer.h"
#include "../../scalar_vector_types.h"
#include "../../global_memory_pool.h"

#include "../coulomb_integral.h"
#include "../aint_init.h"
#include "../aint_common.h"

namespace G2G
{
#include "gpu_vars/g2g_gpu_variables.h"

extern double free_global_memory;
}

using std::cout;
using std::vector;
using std::endl;

namespace AINT
{

#include "gpu_vars/os_gpu_variables.h"
#include "gpu_vars/coulomb_gpu_variables.h"

#include "kernels/os_util.h"
#include "kernels/coulomb_forces.h"
#include "kernels/coulomb_energy.h"
#include "kernels/coulomb_fit.h"

// -----------------------------------------------------------------------------------------------------------------
template <class scalar_type>
bool CoulombIntegral<scalar_type>::load_aux_basis( void )
{
  std::vector<G2G::vec_type<scalar_type,2> > factor_ac_dens_cpu;
  std::vector<scalar_type> fit_dens_cpu;
  std::vector<G2G::vec_type<scalar_type,3> > nuc_dens_cpu;
  std::vector<uint> nuc_ind_dens_cpu;

  input_size = 0;

  for (uint func = 0; func < integral_vars.s_funcs_dens; func++) {
    uint nuc_ind = integral_vars.nucleii_dens(func) - 1;
    double3 nuc_pos = G2G::fortran_vars.atom_positions(nuc_ind);
    for (uint k = 0; k < integral_vars.contractions_dens(func); k++) {
      factor_ac_dens_cpu.push_back(G2G::vec_type<scalar_type,2>(integral_vars.a_values_dens(func,k),integral_vars.c_values_dens(func,k)));
      nuc_dens_cpu.push_back(G2G::vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
      nuc_ind_dens_cpu.push_back(nuc_ind);
      input_ind_cpu.push_back(input_size);
      input_size++;
    }
  }

  s_end = input_size;
  cudaMemcpyToSymbol(gpu_s_end, &s_end, sizeof(s_end), 0, cudaMemcpyHostToDevice);

  uint tmp_size = input_size;
  for (uint func = tmp_size; func < COALESCED_DIMENSION(tmp_size); func++) {
    factor_ac_dens_cpu.push_back(factor_ac_dens_cpu[tmp_size-1]);
    nuc_dens_cpu.push_back(nuc_dens_cpu[tmp_size-1]);
    nuc_ind_dens_cpu.push_back(nuc_ind_dens_cpu[tmp_size-1]);
    input_size++;
  }

  p_offset = input_size;
  cudaMemcpyToSymbol(gpu_p_offset, &p_offset, sizeof(p_offset), 0, cudaMemcpyHostToDevice);

  tmp_size = 0;
  for (uint func = integral_vars.s_funcs_dens; func < integral_vars.s_funcs_dens+integral_vars.p_funcs_dens*3; func += 3) {
    uint nuc_ind = integral_vars.nucleii_dens(func) - 1;
    double3 nuc_pos = G2G::fortran_vars.atom_positions(nuc_ind);
    for (uint k = 0; k < integral_vars.contractions_dens(func); k++) {
      for (uint f = 0; f < 3; f++) {
        factor_ac_dens_cpu.push_back(G2G::vec_type<scalar_type,2>(integral_vars.a_values_dens(func,k),integral_vars.c_values_dens(func,k)));
        nuc_dens_cpu.push_back(G2G::vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
        nuc_ind_dens_cpu.push_back(nuc_ind);
        input_ind_cpu.push_back(input_size+f);
      }
      tmp_size += 3;
      input_size += 3;
      if (tmp_size == 126) {
        for (uint f = 0; f < 2; f++) factor_ac_dens_cpu.push_back(G2G::vec_type<scalar_type,2>(integral_vars.a_values_dens(func,k),integral_vars.c_values_dens(func,k)));
        for (uint f = 0; f < 2; f++) nuc_dens_cpu.push_back(G2G::vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
        for (uint f = 0; f < 2; f++) nuc_ind_dens_cpu.push_back(nuc_ind);
        tmp_size = 0;
        input_size += 2;
      }
    }
  }

  p_end = input_size;
  cudaMemcpyToSymbol(gpu_p_end, &p_end, sizeof(p_end), 0, cudaMemcpyHostToDevice);

  tmp_size = input_size;
  for (uint func = tmp_size; func < COALESCED_DIMENSION(tmp_size); func++) {
    factor_ac_dens_cpu.push_back(factor_ac_dens_cpu[tmp_size-1]);
    nuc_dens_cpu.push_back(nuc_dens_cpu[tmp_size-1]);
    nuc_ind_dens_cpu.push_back(nuc_ind_dens_cpu[tmp_size-1]);
    input_size++;
  }

  d_offset = input_size;
  cudaMemcpyToSymbol(gpu_d_offset, &d_offset, sizeof(d_offset), 0, cudaMemcpyHostToDevice);

  tmp_size = 0;
  for (uint func = integral_vars.s_funcs_dens+integral_vars.p_funcs_dens*3; func < integral_vars.m_dens; func += 6) {
    uint nuc_ind = integral_vars.nucleii_dens(func) - 1;
    double3 nuc_pos = G2G::fortran_vars.atom_positions(nuc_ind);
    for (uint k = 0; k < integral_vars.contractions_dens(func); k++) {
      for (uint f = 0; f < 6; f++) {
        factor_ac_dens_cpu.push_back(G2G::vec_type<scalar_type,2>(integral_vars.a_values_dens(func,k),integral_vars.c_values_dens(func,k)));
        nuc_dens_cpu.push_back(G2G::vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
        nuc_ind_dens_cpu.push_back(nuc_ind);
        input_ind_cpu.push_back(input_size+f);
      }
      tmp_size += 6;
      input_size += 6;
      if (tmp_size == 126) {
        for (uint f = 0; f < 2; f++) factor_ac_dens_cpu.push_back(G2G::vec_type<scalar_type,2>(integral_vars.a_values_dens(func,k),integral_vars.c_values_dens(func,k)));
        for (uint f = 0; f < 2; f++) nuc_dens_cpu.push_back(G2G::vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
        for (uint f = 0; f < 2; f++) nuc_ind_dens_cpu.push_back(nuc_ind);
        tmp_size = 0;
        input_size += 2;
      }
    }
  }

  d_end = input_size;
  cudaMemcpyToSymbol(gpu_d_end, &d_end, sizeof(d_end), 0, cudaMemcpyHostToDevice);

  size_t gpu_size = factor_ac_dens_cpu.size() * 2 * sizeof(scalar_type);
  gpu_size += nuc_dens_cpu.size() * 3 * sizeof(scalar_type);
  gpu_size += (nuc_ind_dens_cpu.size() + input_ind_cpu.size()) * sizeof(uint);
  gpu_size += input_size * sizeof(scalar_type);
  float mb_size = (float)gpu_size / 1048576.0f;
  cout << "Coulomb aux basis input size: " << mb_size << " MB" << endl;
  if (GlobalMemoryPool::tryAlloc(gpu_size) && G2G::free_global_memory > 0.0) return false;

  factor_ac_dens_dev = factor_ac_dens_cpu;
  nuc_dens_dev = nuc_dens_cpu;
  nuc_ind_dens_dev = nuc_ind_dens_cpu;
  input_ind_dev = input_ind_cpu;

  fit_dens_dev.resize(input_size);

  cudaAssertNoError("CoulombIntegral::load_aux_basis");

  return true;
}

template <class scalar_type>
bool CoulombIntegral<scalar_type>::load_input( void )
{
    
    G2G::HostMatrix<double> Ginv_h(COALESCED_DIMENSION(integral_vars.m_dens),integral_vars.m_dens);
    for (uint i = 0; i < integral_vars.m_dens; i++)
    {
      for (uint j = 0; j < i; j++)
      {
        Ginv_h(i,j) = integral_vars.Ginv_input(i+(2*integral_vars.m_dens-(j+1))*j/2);
      }
      for (uint j = i; j < integral_vars.m_dens; j++)
      {
        Ginv_h(i,j) = integral_vars.Ginv_input(j+(2*integral_vars.m_dens-(i+1))*i/2);
      }
    }
    size_t gpu_size = Ginv_h.bytes();
    float mb_size = (float)gpu_size / 1048576.0f;
    cout << "Coulomb input size: " << mb_size << " MB" << endl;
    if (GlobalMemoryPool::tryAlloc(gpu_size) && G2G::free_global_memory > 0.0) return false;
    
    Ginv_dev = Ginv_h;

    cudaMemcpyToSymbol(gpu_out_offsets,os_int.out_offsets,NUM_TERM_TYPES*sizeof(uint),0,cudaMemcpyHostToDevice);

    cudaAssertNoError("CoulombIntegral::load_input");

    return true;
}

template <class scalar_type>
bool CoulombIntegral<scalar_type>::alloc_output( void )
{
    uint partial_out_size = 0;
    // Output arrays (probably) don't need to be padded for alignment as only one (or three) threads per block write to them
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      uint this_count = divUp(os_int.term_type_counts[i],QMMM_BLOCK_SIZE);
      partial_out_size += this_count;
    }
    size_t gpu_size = COALESCED_DIMENSION(integral_vars.m_dens) * partial_out_size * sizeof(double);
    float mb_size = (float)gpu_size / 1048576.0f;
    cout << "Coulomb output size: " << mb_size << " MB" << endl;
    if (GlobalMemoryPool::tryAlloc(gpu_size) && G2G::free_global_memory > 0.0) return false;

    rc_partial_dev.resize(COALESCED_DIMENSION(integral_vars.m_dens),partial_out_size);

    cudaAssertNoError("CoulombIntegral::alloc_output");

    return true;
}

template <class scalar_type>
void CoulombIntegral<scalar_type>::clear( void )
{
    input_ind_cpu.clear();

    if (fit_dens_dev.is_allocated() && G2G::free_global_memory > 0.0) {
      size_t gpu_size = fit_dens_dev.bytes() + factor_ac_dens_dev.bytes() + nuc_dens_dev.bytes() + nuc_ind_dens_dev.bytes() + input_ind_dev.bytes();
      GlobalMemoryPool::dealloc(gpu_size);
    }
    fit_dens_dev.deallocate();
    factor_ac_dens_dev.deallocate();
    nuc_dens_dev.deallocate();
    nuc_ind_dens_dev.deallocate();
    input_ind_dev.deallocate();
    if (Ginv_dev.is_allocated() && G2G::free_global_memory > 0.0) {
      size_t gpu_size = Ginv_dev.bytes();
      GlobalMemoryPool::dealloc(gpu_size);
    }
    Ginv_dev.deallocate();
    if (rc_partial_dev.is_allocated() && G2G::free_global_memory > 0.0) {
      size_t gpu_size = rc_partial_dev.bytes();
      GlobalMemoryPool::dealloc(gpu_size);
    }
    rc_partial_dev.deallocate();

    cudaAssertNoError("CoulombIntegral::clear");
}

template <class scalar_type>
void CoulombIntegral<scalar_type>::calc_gradient( double* qm_forces, bool cpu_fit_dens )
{

    os_int.reload_density();

    if (cpu_fit_dens) {
    G2G::HostMatrix<scalar_type> fit_dens_h(input_size);
    for (uint i = 0; i < input_size; i++) fit_dens_h(i) = 0.0;
    for (uint i = 0; i < integral_vars.m_dens; i++) {
      fit_dens_h(input_ind_cpu[i]) = integral_vars.af_input_ndens1(i);
    }
    fit_dens_dev = fit_dens_h;
    }

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

#define coulomb_forces_parameters \
  os_int.term_type_counts[i], os_int.factor_ac_dev.data, os_int.nuc_dev.data, os_int.dens_values_dev.data+dens_offset, os_int.func_code_dev.data+offset,os_int.local_dens_dev.data+offset, \
  os_int.partial_qm_forces_dev.data+force_offset, COALESCED_DIMENSION(partial_out_size),factor_ac_dens_dev.data,nuc_dens_dev.data,nuc_ind_dens_dev.data,fit_dens_dev.data, \
  s_end, p_end, d_end, p_offset, d_offset
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
      switch (i) {
        case 0: gpu_coulomb_forces<scalar_type,0><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 1: gpu_coulomb_forces<scalar_type,1><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 2: gpu_coulomb_forces<scalar_type,2><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 3: gpu_coulomb_forces<scalar_type,3><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 4: gpu_coulomb_forces<scalar_type,4><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
        case 5: gpu_coulomb_forces<scalar_type,5><<<gridSize,blockSize,0,stream[i]>>>( coulomb_forces_parameters ); break;
      }
    }
    cudaDeviceSynchronize();
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }

    cudaUnbindTexture(str_tex);

    os_int.get_gradient_output(qm_forces, partial_out_size);

    cudaAssertNoError("CoulombIntegral::calc_gradient");
}

template<class scalar_type>
void CoulombIntegral<scalar_type>::fit_aux_density( void )
{
    os_int.reload_density();
    uint partial_out_size = 0, max_partial_size = 0;
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      uint this_count = divUp(os_int.term_type_counts[i],QMMM_BLOCK_SIZE);
      if (this_count > max_partial_size) max_partial_size = this_count;
      partial_out_size += this_count;
    }

    {
      dim3 threads(input_size);
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads,blockSize);
      zero_fock<scalar_type><<<gridSize,blockSize>>>(fit_dens_dev.data,input_size,1);
    }
    {
      dim3 threads(COALESCED_DIMENSION(integral_vars.m_dens),partial_out_size);
      dim3 blockSize(32,4);
      dim3 gridSize = divUp(threads,blockSize);
      zero_fock<double><<<gridSize,blockSize>>>(rc_partial_dev.data,COALESCED_DIMENSION(integral_vars.m_dens),partial_out_size);
    }

    //
    // The STR table for F(m,U) calculation is being accessed via texture fetches
    //
    cudaBindTextureToArray(str_tex,gammaArray);

#define fit1_parameters \
  os_int.term_type_counts[i], os_int.factor_ac_dev.data, os_int.nuc_dev.data, os_int.dens_values_dev.data+dens_offset, os_int.func_code_dev.data+offset,os_int.local_dens_dev.data+offset, \
  rc_partial_dev.data+rc_offset, COALESCED_DIMENSION(integral_vars.m_dens),factor_ac_dens_dev.data,nuc_dens_dev.data, \
  s_end, p_end, d_end, p_offset, d_offset
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
      uint rc_offset = os_int.out_offsets[i] * COALESCED_DIMENSION(integral_vars.m_dens);
      dim3 threads = os_int.term_type_counts[i];
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads, blockSize);
      switch (i) {
        case 0: gpu_coulomb_fit1<scalar_type,0><<<gridSize,blockSize,0,stream[i]>>>( fit1_parameters ); break;
        case 1: gpu_coulomb_fit1<scalar_type,1><<<gridSize,blockSize,0,stream[i]>>>( fit1_parameters ); break;
        case 2: gpu_coulomb_fit1<scalar_type,2><<<gridSize,blockSize,0,stream[i]>>>( fit1_parameters ); break;
        case 3: gpu_coulomb_fit1<scalar_type,3><<<gridSize,blockSize,0,stream[i]>>>( fit1_parameters ); break;
        case 4: gpu_coulomb_fit1<scalar_type,4><<<gridSize,blockSize,0,stream[i]>>>( fit1_parameters ); break;
        case 5: gpu_coulomb_fit1<scalar_type,5><<<gridSize,blockSize,0,stream[i]>>>( fit1_parameters ); break;
      }
      dim3 reduceThreads = integral_vars.m_dens;
      dim3 reduceBlockSize(QMMM_REDUCE_BLOCK_SIZE);
      dim3 reduceGridSize = divUp(reduceThreads,reduceBlockSize);
      gpu_coulomb_rc_term_reduce<scalar_type><<<reduceGridSize,reduceBlockSize,0,stream[i]>>>( rc_partial_dev.data,COALESCED_DIMENSION(integral_vars.m_dens),
                                                                                               os_int.out_offsets[i], (i<NUM_TERM_TYPES-1?os_int.out_offsets[i+1]:partial_out_size) ,integral_vars.m_dens );
    }
    cudaDeviceSynchronize();
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }

    {
      dim3 reduceThreads = integral_vars.m_dens;
      dim3 reduceBlockSize(QMMM_REDUCE_BLOCK_SIZE);
      dim3 reduceGridSize = divUp(reduceThreads,reduceBlockSize);
      gpu_coulomb_rc_reduce<scalar_type><<<reduceGridSize,reduceBlockSize>>>( rc_partial_dev.data,COALESCED_DIMENSION(integral_vars.m_dens),integral_vars.m_dens,NUM_TERM_TYPES );
    }
    cudaDeviceSynchronize();
    {
      dim3 threads = integral_vars.m_dens;
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads,blockSize);
      gpu_coulomb_fit2<scalar_type><<<gridSize,blockSize>>>( rc_partial_dev.data,fit_dens_dev.data,Ginv_dev.data,input_ind_dev.data,integral_vars.m_dens );
    }
    cudaDeviceSynchronize();

    G2G::HostMatrix<scalar_type> fit_dens_h(fit_dens_dev);
    for (uint i = 0; i < integral_vars.m_dens; i++) {
      integral_vars.af_input_ndens1(i) = fit_dens_h(input_ind_cpu[i]);
    }

    //cudaUnbindTexture(str_tex);

    cudaAssertNoError("CoulombIntegral::fit_aux_density");
}

template<class scalar_type>
void CoulombIntegral<scalar_type>::calc_fock( double& Es )
{
    uint partial_out_size = 0, max_partial_size = 0;
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      uint this_count = divUp(os_int.term_type_counts[i],QMMM_BLOCK_SIZE);
      if (this_count > max_partial_size) max_partial_size = this_count;
      partial_out_size += this_count;
    }

    {
      // I don't know why, but doing this zeroing kernel with this grid (rather than a more logical 1D grid)
      // actually seems to speed things up...I can only guess it's due to some synchronization issue
      dim3 threads(os_int.dens_values.size(),max_partial_size);
      dim3 blockSize(32,4);//QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads,blockSize);
      //
      // Zero the partial Fock matrix on the GPU
      //
      zero_fock<double><<<gridSize,blockSize>>>(os_int.partial_fock_dev.data,os_int.dens_values.size(),1);
    }

    //
    // The STR table for F(m,U) calculation is being accessed via texture fetches
    //
    //cudaBindTextureToArray(str_tex,gammaArray);

#define coulomb_fock_parameters \
  os_int.term_type_counts[i], os_int.factor_ac_dev.data, os_int.nuc_dev.data, os_int.func_code_dev.data+offset,os_int.local_dens_dev.data+offset, \
  os_int.partial_fock_dev.data+fock_offset, os_int.dens_values.size(),factor_ac_dens_dev.data,nuc_dens_dev.data,fit_dens_dev.data, \
  s_end, p_end, d_end, p_offset, d_offset
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
        case 0: gpu_coulomb_fock<scalar_type,0><<<gridSize,blockSize,0,stream[i]>>>( coulomb_fock_parameters ); break;
        case 1: gpu_coulomb_fock<scalar_type,1><<<gridSize,blockSize,0,stream[i]>>>( coulomb_fock_parameters ); break;
        case 2: gpu_coulomb_fock<scalar_type,2><<<gridSize,blockSize,0,stream[i]>>>( coulomb_fock_parameters ); break;
        case 3: gpu_coulomb_fock<scalar_type,3><<<gridSize,blockSize,0,stream[i]>>>( coulomb_fock_parameters ); break;
        case 4: gpu_coulomb_fock<scalar_type,4><<<gridSize,blockSize,0,stream[i]>>>( coulomb_fock_parameters ); break;
        case 5: gpu_coulomb_fock<scalar_type,5><<<gridSize,blockSize,0,stream[i]>>>( coulomb_fock_parameters ); break;
      }

      //
      // Reduce the partial Fock terms for a particular term type as soon as that kernel is done; also calculate partial energies for that type
      //
      dim3 reduceThreads = os_int.dens_counts[i];
      dim3 reduceBlockSize(QMMM_REDUCE_BLOCK_SIZE);
      dim3 reduceGridSize = divUp(reduceThreads,reduceBlockSize);
      gpu_fock_reduce<scalar_type><<<reduceGridSize,reduceBlockSize,0,stream[i]>>>( os_int.partial_fock_dev.data+fock_offset, os_int.dens_values_dev.data+fock_offset,
                                                                                      os_int.partial_energies_dev.data+os_int.energies_offsets[i],
                                                                                      os_int.dens_values.size(), max_partial_size, os_int.dens_counts[i] );
    }
    cudaDeviceSynchronize();
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }

    cudaUnbindTexture(str_tex);

    os_int.get_fock_output(Es, G2G::fortran_vars.rmm_output);

    cudaAssertNoError("CoulombIntegral::calc_fock");
}

#if AINT_MP && !FULL_DOUBLE
template class CoulombIntegral<float>;
#else
template class CoulombIntegral<double>;
#endif

}
