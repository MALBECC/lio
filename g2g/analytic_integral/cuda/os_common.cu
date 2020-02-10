#include <fstream>
#include <iostream>
#include <vector>

#include "../../common.h"
#include "../../cuda/cuda_extra.h"
#include "../../init.h"
#include "../../matrix.h"
#include "../../timer.h"
#include "../../scalar_vector_types.h"
#include "../../global_memory_pool.h"

#include "../os_integral.h"
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

// On-device global data that is common to all analytic (using O-S) integral evaluation
namespace AINT
{

cudaArray* gammaArray;
__device__ __constant__ uint gpu_m;
#if !AINT_MP || FULL_DOUBLE
texture<int2, cudaTextureType2D, cudaReadModeElementType> str_tex; // Texture for STR array (used in F(m,U))
__device__ __constant__ double gpu_fac[17];
#else
texture<float, cudaTextureType2D, cudaReadModeElementType> str_tex;
__device__ __constant__ float gpu_fac[17];
#endif

__device__ __constant__ uint TERM_TYPE_GAUSSIANS[6] = { 1, 3, 9, 6, 18, 36 };
__device__ __constant__ uint gpu_atom_types[MAX_ATOMS];
__device__ __constant__ uint gpu_atom_Z[MAX_ATOMS];

template<class scalar_type>
void OSIntegral<scalar_type>::load_params(void)
{
    // Use the highest cc device available for aint stuff
    int devcount = cudaGetGPUCount();
    int devnum = -1, devmajor = -1, devminor = -1;
    cudaDeviceProp devprop;
    for (int i = 0; i < devcount; i++) {
      if (cudaGetDeviceProperties(&devprop, i) != cudaSuccess) throw std::runtime_error("Could not get device properties!");
      if (devprop.major > devmajor || (devprop.major == devmajor && devprop.minor > devminor)) {
        devnum = i;
        devmajor = devprop.major;
        devminor = devprop.minor;
      }
    }

    if (G2G::verbose > 3) cout << "  Using device " << devnum << " for analytic integral calculations." << endl;
    this->my_device = devnum;
    int previous_device; cudaGetDevice(&previous_device);
    if(cudaSetDevice(devnum) != cudaSuccess) std::cout << "Error: can't set the device " << devnum << std::endl;

    cudaMemcpyToSymbol(gpu_m, &G2G::fortran_vars.m, sizeof(G2G::fortran_vars.m), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(gpu_atom_types, G2G::fortran_vars.atom_types.data, G2G::fortran_vars.atom_types.bytes(), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(gpu_atom_Z, integral_vars.atom_Z.data, integral_vars.atom_Z.bytes(), 0, cudaMemcpyHostToDevice);
    //
    // Set up arrays needed for F(m,U) calculation in QM/MM kernels (STR and FAC) and send them to device
    // FAC is small so it's put into constant memory
    // STR is large and accessed (potentially) with a random access pattern in the first index
    // TODO: putting STR into texture for now; need to see if there's a better way to access it in the kernel
    //

    // Cast STR/FAC to appropriate type (float/double)
    G2G::HostMatrix<scalar_type> h_str(880,22), h_fac(17);
    for (uint i = 0; i < 880; i++) {
      for (uint j = 0; j < 22; j++) {
        h_str(i,j) = integral_vars.str(i,j);
      }
    }
    for (uint i = 0; i < 17; i++) {
      h_fac(i) = integral_vars.fac(i);
    }

    str_tex.normalized = false;
    str_tex.filterMode = cudaFilterModePoint;

    cudaMallocArray(&gammaArray,&str_tex.channelDesc,880,22);
    cudaMemcpyToArray(gammaArray,0,0,h_str.data,sizeof(scalar_type)*880*22,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(gpu_fac,h_fac.data,h_fac.bytes(),0,cudaMemcpyHostToDevice);

    cudaSetDevice(previous_device);
    cudaAssertNoError("OSIntegral::load_params");

}

template<class scalar_type>
void OSIntegral<scalar_type>::clear( void )
{
    func_code.clear();
    local_dens.clear();
    dens_values.clear();
    local2globaldens.clear();

    if (factor_ac_dev.is_allocated() && G2G::free_global_memory > 0.0) {
      size_t gpu_size = factor_ac_dev.bytes() + nuc_dev.bytes() + func_code_dev.bytes() + local_dens_dev.bytes() + dens_values_dev.bytes();
      GlobalMemoryPool::dealloc(gpu_size);
    }
    factor_ac_dev.deallocate();
    nuc_dev.deallocate();
    func_code_dev.deallocate();
    local_dens_dev.deallocate();
    dens_values_dev.deallocate();

    if (partial_qm_forces_dev.is_allocated() && G2G::free_global_memory > 0.0) {
      size_t gpu_size = partial_fock_dev.bytes() + partial_energies_dev.bytes() + partial_qm_forces_dev.bytes();
      GlobalMemoryPool::dealloc(gpu_size);
    }
    partial_fock_dev.deallocate();
    partial_energies_dev.deallocate();
    partial_qm_forces_dev.deallocate();

    cudaAssertNoError("OSIntegral::clear");
}

template<class scalar_type>
void OSIntegral<scalar_type>::deinit( void )
{
    clear();

    cudaFreeArray(gammaArray);

    cudaAssertNoError("OSIntegral::deinit");
}

template<class scalar_type>
bool OSIntegral<scalar_type>::load_input( void )
{
    //
    // Set up device arrays for function values and mapping function -> nuclei
    //
    // TODO: tried contracting the nuc / function value arrays to a size of total_funcs rather than m, seems to slow down the kernel
    // Doing this requires some more indexing math in the kernel, but need to test on bigger test case
    G2G::HostMatrix<G2G::vec_type<scalar_type, 2> > factor_ac_cpu(COALESCED_DIMENSION(G2G::fortran_vars.m), MAX_CONTRACTIONS);
    G2G::HostMatrixUInt nuc_cpu(G2G::fortran_vars.m, 1);

    uint localfunc = 0, func = 0;
    while (func < G2G::fortran_vars.m) {
      nuc_cpu(localfunc) = G2G::fortran_vars.nucleii(func) - 1;
      for (uint k = 0; k < G2G::fortran_vars.contractions(func); k++) {
        factor_ac_cpu(localfunc, k) = G2G::vec_type<scalar_type, 2>(G2G::fortran_vars.a_values(func, k), G2G::fortran_vars.c_values(func, k));
      }
      func++;
      localfunc++;
    }
    size_t gpu_size = factor_ac_cpu.bytes() + nuc_cpu.bytes();
    gpu_size += (func_code.size() + local_dens.size()) * sizeof(uint);
    gpu_size += dens_values.size() * sizeof(scalar_type);
//    float mb_size = (float)gpu_size / 1048576.0f;
//    cout << "O-S common input size: " << mb_size << " MB" << endl;
    if (GlobalMemoryPool::tryAlloc(gpu_size) && G2G::free_global_memory > 0.0) return false;

    factor_ac_dev = factor_ac_cpu;
    nuc_dev = nuc_cpu;

    //
    // Send input arrays (thread->primitive map and thread->density map) to the device
    //
    func_code_dev = func_code;
    local_dens_dev = local_dens; // TODO: this is probably not necessary here as reload_density gets called before all the major kernel calls
    //
    // Send reduced density matrix to the device
    //
    dens_values_dev = dens_values;

    cudaAssertNoError("OSIntegral::load_input");

    return true;
}

template<class scalar_type>
void OSIntegral<scalar_type>::reload_density( void )
{
    for (uint i = 0; i < dens_values.size(); i++)
    {
      dens_values[i] = G2G::fortran_vars.rmm_input_ndens1.data[local2globaldens[i]];
    }
    dens_values_dev = dens_values;

    cudaAssertNoError("OSIntegral::reload_density");
}

template<class scalar_type>
bool OSIntegral<scalar_type>::alloc_output( void )
{
    //
    // Allocate output arrays on device
    // Currently, each block in the kernel reduces its output, so the output arrays have length (# block)
    //
    uint partial_out_size = 0, max_partial_size = 0;
    out_offsets[0] = 0;
    // Output arrays (probably) don't need to be padded for alignment as only one (or three) threads per block write to them
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      uint this_count = divUp(term_type_counts[i],QMMM_BLOCK_SIZE);
      if (this_count > max_partial_size) max_partial_size = this_count;
      partial_out_size += this_count;
      if (i+1<NUM_TERM_TYPES) { out_offsets[i+1] = partial_out_size; }
    }
    //
    // When calculating energies, the energy gets reduced per-block; we figure out the offets/counts of different term types into the partial output energy array here
    //
    uint energies_size = 0;
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      energies_offsets[i] = energies_size;
      energies_size += divUp(dens_counts[i],QMMM_REDUCE_BLOCK_SIZE);
    }

    size_t gpu_size = COALESCED_DIMENSION(partial_out_size) * G2G::fortran_vars.atoms * 3 * sizeof(scalar_type);
    gpu_size += (dens_values.size() /* max_partial_size*/ + energies_size) * sizeof(double);
//    float mb_size = (float)gpu_size / 1048576.0f;
//    cout << "O-S common output size: " << mb_size << " MB" << endl;
    if (GlobalMemoryPool::tryAlloc(gpu_size) && G2G::free_global_memory > 0.0) return false;

    //
    // Forces: output is partial forces
    //
    partial_qm_forces_dev.resize(COALESCED_DIMENSION(partial_out_size), G2G::fortran_vars.atoms);
    //
    // Fock: ouptut is partial Fock elements
    //
    // The partial Fock matrix is partitioned by term type, so the second (partial) dimension needs to be as big as largest count of a single term type
    partial_fock_dev.resize(dens_values.size(),1);//max_partial_size);

    partial_energies_dev.resize(energies_size,1);

    cudaAssertNoError("OSIntegral::alloc_output");

    return true;
}

template<class scalar_type>
void OSIntegral<scalar_type>::get_fock_output( double& Es, G2G::FortranMatrix<double>& fock_out )
{
    uint energies_size = 0;
    for (uint i = 0; i < NUM_TERM_TYPES; i++) {
      energies_size += divUp(dens_counts[i],QMMM_REDUCE_BLOCK_SIZE);
    }

    //
    // Download reduced Fock matrix and partially reduced energies
    // The Fock matrix has been reduced to the first row of the output, so we only want that much of the device array
    //
    G2G::HostMatrix<double> cpu_fock(partial_fock_dev);
    //cudaMemcpy(cpu_fock.data,partial_fock_dev.data,cpu_fock.bytes(),cudaMemcpyDeviceToHost);
    G2G::HostMatrix<double> cpu_partial_energies(partial_energies_dev);

    //
    // Send Fock elements back to RMM(M11) and do final reduction of e-nuc energies into Es
    //
    for (uint t = 0; t < NUM_TERM_TYPES; t++) {
      for (uint i = dens_offsets[t]; i < dens_offsets[t] + dens_counts[t]; i++) {
        uint dens_ind = local2globaldens[i];
        fock_out(dens_ind) += cpu_fock(i);
      }
    }
    for (uint i = 0; i < energies_size; i++) {
      Es += cpu_partial_energies(i);
    }

    cudaAssertNoError("OSIntegral::get_fock_output");
}

template<class scalar_type>
void OSIntegral<scalar_type>::get_gradient_output(double* qm_forces, uint partial_out_size)
{
    G2G::HostMatrix<G2G::vec_type<scalar_type,3> > cpu_partial_qm_forces(partial_qm_forces_dev);

    //
    // Accumulate partial output
    //
    // TODO: need to think about how to accumulate individual force terms
    // Currently, we reduce on a per-block basis in the kernel, then accumulate the block results here on the host
    //
    // The energy partial results are being reduced on-device and that works very well, could probably do that for forces too
    //
    for (uint i = 0; i < G2G::fortran_vars.atoms; i++) {
      for (uint j = 0; j < partial_out_size; j++) {
        qm_forces[i + 0 * G2G::fortran_vars.atoms] += cpu_partial_qm_forces(j,i).x;
        qm_forces[i + 1 * G2G::fortran_vars.atoms] += cpu_partial_qm_forces(j,i).y;
        qm_forces[i + 2 * G2G::fortran_vars.atoms] += cpu_partial_qm_forces(j,i).z;
      }
    }

    cudaAssertNoError("OSIntegral::get_gradient_output");
}

#if AINT_MP && !FULL_DOUBLE
template class OSIntegral<float>;
#else
template class OSIntegral<double>;
#endif

}
