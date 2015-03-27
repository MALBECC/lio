#include <iostream>
#include <vector>

#include "../common.h"
#include "../init.h"
#include "cuda_extra.h"
#include "../matrix.h"
#include "../timer.h"
#include "../scalar_vector_types.h"

using std::cout;
using std::vector;
using std::endl;

namespace G2G {
#if FULL_DOUBLE
texture<int2, cudaTextureType2D, cudaReadModeElementType> qmmm_str_tex; // Texture for STR array (used in F(m,U))
#else
texture<float, cudaTextureType2D, cudaReadModeElementType> qmmm_str_tex;
#endif

//#include "gpu_variables.h"
#include "qmmm_coul_gpu_variables.h"
#include "kernels/qmmm_forces.h"
#include "kernels/qmmm_energy.h"
#include "kernels/coulomb_forces.h"

void gpu_set_qmmm_coul_variables(void) {
  int previous_device; cudaGetDevice(&previous_device);
  int gpu_devices = cudaGetGPUCount();
  for(int i = 0; i < gpu_devices; i++) {
    if(cudaSetDevice(i) != cudaSuccess)
      std::cout << "Error: can't set the device " << i << std::endl;
    cudaMemcpyToSymbol(gpu_m, &fortran_vars.m, sizeof(fortran_vars.m), 0, cudaMemcpyHostToDevice);

    // This is needed by d-d QM/MM forces calculations to know which orbital a thread maps to
    uint d_offset = fortran_vars.s_funcs + fortran_vars.p_funcs*3;
    cudaMemcpyToSymbol(gpu_d_offset, &d_offset, sizeof(d_offset), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(gpu_dens_gauss, &fortran_vars.gaussians_dens, sizeof(fortran_vars.gaussians_dens), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(gpu_dens_s_gauss, &fortran_vars.s_gaussians_dens, sizeof(fortran_vars.s_gaussians_dens), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(gpu_dens_p_gauss, &fortran_vars.p_gaussians_dens, sizeof(fortran_vars.p_gaussians_dens), 0, cudaMemcpyHostToDevice);
  }
  cudaSetDevice(previous_device);
  cudaAssertNoError("gpu_set_qmmm_coul_variables");
}

//
// Set up arrays needed for F(m,U) calculation in QM/MM kernels (STR and FAC) and send them to device
// FAC is small so it's put into constant memory
// STR is large and accessed (potentially) with a random access pattern in the first index
// TODO: putting STR into texture for now; need to see if there's a better way to access it in the kernel
//
template<class scalar_type>
void gpu_set_gamma_arrays()
{
  // Cast STR/FAC to appropriate type (float/double)
  HostMatrix<scalar_type> h_str(880,22), h_fac(17);
  for (uint i = 0; i < 880; i++) {
    for (uint j = 0; j < 22; j++) {
      h_str(i,j) = fortran_vars.str(i,j);
    }
  }
  for (uint i = 0; i < 17; i++) {
    h_fac(i) = fortran_vars.fac(i);
  }

  qmmm_str_tex.normalized = false;
  qmmm_str_tex.filterMode = cudaFilterModePoint;
  int previous_device; cudaGetDevice(&previous_device);
  int gpu_devices = cudaGetGPUCount();
  gammaArrays = new cudaArray*[gpu_devices];
  for(int i = 0; i < gpu_devices; i++) {
    if(cudaSetDevice(i) != cudaSuccess)
      std::cout << "Error: can't set the device " << i << std::endl;
    cudaMallocArray(&gammaArrays[i],&qmmm_str_tex.channelDesc,880,22);//GAMMA_LENGTH,6);
    cudaMemcpyToArray(gammaArrays[i],0,0,h_str.data,sizeof(scalar_type)*880*22,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(gpu_fac,h_fac.data,h_fac.bytes(),0,cudaMemcpyHostToDevice);
  }
  cudaSetDevice(previous_device);
  cudaAssertNoError("gpu_set_gamma_arrays");
}

//
// Tell the GPU how many MM atoms are being used this step
//
void gpu_set_clatoms(void)
{
  int previous_device; cudaGetDevice(&previous_device);
  int gpu_devices = cudaGetGPUCount();
  for(int i = 0; i < gpu_devices; i++) {
    if(cudaSetDevice(i) != cudaSuccess)
      std::cout << "Error: can't set the device " << i << std::endl;
  cudaMemcpyToSymbol(gpu_clatoms, &fortran_vars.clatoms, sizeof(fortran_vars.clatoms), 0, cudaMemcpyHostToDevice);
  }

  cudaSetDevice(previous_device);
  cudaAssertNoError("gpu_set_clatoms");
}

#if FULL_DOUBLE
template void gpu_set_gamma_arrays<double>( void );
#else
template void gpu_set_gamma_arrays<float>( void );
#endif

#define NUM_TERM_TYPES 6 // 6 types when using s,p,and d functions: s-s,p-s,p-p,d-s,d-p,d-d
#define MAX_TERM_TYPE 6

//
// Main QM/MM routine
// If forces = true, calculate gradients of QM/MM operator and return in qm_forces and mm_forces
// If forces = false, calculate Fock matrix elements of QM/MM operator (returned in RMM(M11) back in Fortran)
//                    and QM/MM energy of the current density (nuc-nuc in Ens, e-nuc in Es)
//
template <class scalar_type,bool forces> void g2g_qmmm(double* qm_forces, double* mm_forces, double& Ens, double& Es)
{
  int devnum = -1;
  if (cudaGetDevice(&devnum) != cudaSuccess) throw std::runtime_error("Could not get device number! (QM/MM)");
  uint i,j,ni,nj;
  uint i_orbitals, j_orbitals;
  uint nuc_i,nuc_j;
  vec_type<double,3> A,B,AmB;
  double ai,aj;
  double dsq,ksi,zeta;
  uint num_terms=0, total_num_terms = 0;
  std::vector<uint> func_code, local_dens;

  std::vector<uint> local2globaldens; // Maps reduced density/Fock -> full density/Fock
  std::vector<scalar_type> dens_values;

  uint term_type_counts[NUM_TERM_TYPES]; // Number of threads for a particular type of term (0 = s-s, 1 = p-s, etc)
  uint term_type_offsets[NUM_TERM_TYPES]; // Offsets into the input arrays for each term type
  term_type_offsets[0] = 0; // s-s starts at 0
  uint i_begin, i_end, j_begin, j_end;
  uint tmp_ind = 0;
  uint s_start = 0, p_start = fortran_vars.s_funcs, d_start = fortran_vars.s_funcs + fortran_vars.p_funcs*3, m = fortran_vars.m;
  uint dd_orb;
  if (forces) { dd_orb = 1; } // 1 for d-d means 6 threads use one of dxx, dyx, etc for func i
  else        { dd_orb = 6; }

  //                                       s-s        p-s        p-p        d-s        d-p        d-d
  uint i_begin_vals[MAX_TERM_TYPE]   = { s_start,   p_start,   p_start,   d_start,   d_start,   d_start};
  uint i_end_vals[MAX_TERM_TYPE]     = { p_start,   d_start,   d_start,   m,         m,         m      };
  uint j_begin_vals[MAX_TERM_TYPE]   = { s_start,   s_start,   p_start,   s_start,   p_start,   d_start};
  uint j_end_vals[MAX_TERM_TYPE]     = { p_start-1, p_start-1, d_start-1, p_start-1, d_start-1, m-1    };
  uint i_orbital_vals[MAX_TERM_TYPE] = { 1,         3,         3,         6,         6,         dd_orb };
  uint j_orbital_vals[MAX_TERM_TYPE] = { 1,         1,         3,         1,         3,         6      };

  uint local_dens_ind, num_dens_terms = 0, total_dens_terms = 0;
  uint dens_counts[NUM_TERM_TYPES], dens_offsets[NUM_TERM_TYPES];
  dens_offsets[0] = 0;
  uint tmp_dens_ind = 0;

  Timer nuc,check,prep,kernel,down,reduce;

  //
  // First, the energy/gradient of the nuclear-nuclear interaction between QM and MM centers
  //
  nuc.start_and_sync();
  if (forces) {
    for (i = 0; i < fortran_vars.atoms; i++) {
      double3 qm_pos = fortran_vars.atom_positions(i);
      for (j = 0; j < fortran_vars.clatoms; j++) {
        double3 mm_pos = fortran_vars.clatom_positions(j);
        double3 diff = qm_pos - mm_pos;
        double dist = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
        dist = sqrt(dist);

        double prefactor = -fortran_vars.clatom_charges(j) * (fortran_vars.atom_types(i)+1) / pow(dist,3.0);
        qm_forces[i + 0 * fortran_vars.atoms] += prefactor * diff.x;
        qm_forces[i + 1 * fortran_vars.atoms] += prefactor * diff.y;
        qm_forces[i + 2 * fortran_vars.atoms] += prefactor * diff.z;
        mm_forces[j + 0 * fortran_vars.clatoms] -= prefactor * diff.x;
        mm_forces[j + 1 * fortran_vars.clatoms] -= prefactor * diff.y;
        mm_forces[j + 2 * fortran_vars.clatoms] -= prefactor * diff.z;
      }
    }
  } else {
    for (i = 0; i < fortran_vars.atoms; i++) {
      double3 qm_pos = fortran_vars.atom_positions(i);
      for (j = 0; j < fortran_vars.clatoms; j++) {
        double3 mm_pos = fortran_vars.clatom_positions(j);
        double3 diff = qm_pos - mm_pos;
        double dist = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
        dist = sqrt(dist);

        double E = fortran_vars.clatom_charges(j) * (fortran_vars.atom_types(i)+1) / dist;
        Ens += E;
      }
    }
  }
  nuc.pause();

  //
  // Do check between all basis primitives to find those with significant overlap
  // Check the resulting Gaussian argument from two primitives to the rmax parameter; only use primitives within that cut-off
  //
  // A single thread gets mapped to a pair of significant primitives
  // We set up here arrays that tell which two functions/two primitives a thread is calculating
  // We also pick out the density matrix elements for significant functions here
  //
  check.start();
  for (uint current_term_type = 0; current_term_type < NUM_TERM_TYPES; current_term_type++) {

    term_type_counts[current_term_type] = 0;
    i_begin = i_begin_vals[current_term_type]; i_end = i_end_vals[current_term_type];
    j_begin = j_begin_vals[current_term_type]; j_end = j_end_vals[current_term_type];
    i_orbitals = i_orbital_vals[current_term_type];
    j_orbitals = j_orbital_vals[current_term_type];

    dens_counts[current_term_type] = 0;
    local_dens_ind = 0;

    // We pad the input arrays between term types, so the offsets for each term type need to be tracked
    if (current_term_type > 0) {
      tmp_ind += COALESCED_DIMENSION(term_type_counts[current_term_type-1]);
      term_type_offsets[current_term_type] = tmp_ind;
      tmp_dens_ind += COALESCED_DIMENSION(dens_counts[current_term_type-1]);
      dens_offsets[current_term_type] = tmp_dens_ind;
    }

    // function i, center A
    for (i = i_begin; i < i_end; i += i_orbitals) {
      nuc_i = fortran_vars.nucleii(i) - 1;
      A = fortran_vars.atom_positions(nuc_i);
      // function j, center B
      for (j = j_begin; j <= ((i > j_end)? j_end : i); j += j_orbitals) {
        nuc_j = fortran_vars.nucleii(j) - 1;
        B = fortran_vars.atom_positions(nuc_j);
        AmB = A - B;
        dsq = length2(AmB);
        bool use_funcs = false; // Do these two functions have any significant primitive pairs?
        // primitive ni, function i
        for (ni = 0; ni < fortran_vars.contractions(i); ni++) {
          // primitive nj, function j
          for (nj = 0; nj < fortran_vars.contractions(j); nj++) {
            ai = fortran_vars.a_values(i,ni);
            aj = fortran_vars.a_values(j,nj);
            zeta = ai + aj;
            ksi = ai * aj / zeta;
            total_num_terms += (i==j)? i_orbitals*(i_orbitals+1)/2 : i_orbitals * j_orbitals;

            if (dsq*ksi < fortran_vars.rmax) {
              use_funcs = true;
              num_terms++;
              term_type_counts[current_term_type]++;

              // Encoding which two primitives this thread will calculate a force term for in one number
              // NOTE: with a long integer, we can use this scheme up to an m of about 9000
              //       if m needs to ever be larger than that, we need to break this up into multiple arrays
              uint this_func_code = nj;                                                         // First, primitive # nj in the lowest part
              this_func_code     += ni * MAX_CONTRACTIONS;                                      // Primitive # ni after the space for nj
              this_func_code     += j  * MAX_CONTRACTIONS * MAX_CONTRACTIONS;                   // Function # j after space for primitives
              this_func_code     += i  * MAX_CONTRACTIONS * MAX_CONTRACTIONS * fortran_vars.m;  // Finally, function # i in the highest part

              func_code.push_back(this_func_code); // Which primitives the thread represents
              local_dens.push_back(local_dens_ind); // Which part of the (reduced) density matrix the thread needs

            }
          }
        }

        total_dens_terms += (i==j)? i_orbitals*(i_orbitals+1)/2 : i_orbitals * j_orbitals;
        // dens_values is a reduced density matrix that only keeps the elements of functions with significant primitive pairs
        // local2globaldens maps from a reduced density/Fock index back to the full matrix
        if (use_funcs) {
          for (uint i_orbital = 0; i_orbital < i_orbitals; i_orbital++) {
            uint j_orbital_finish = (i==j)? i_orbital+1 : j_orbitals;
            for (uint j_orbital = 0; j_orbital < j_orbital_finish; j_orbital++) {
              num_dens_terms++;
              dens_counts[current_term_type]++;

              uint dens_ind = (i+i_orbital) + (2*fortran_vars.m-((j+j_orbital)+1))*(j+j_orbital)/2;
              dens_values.push_back(fortran_vars.rmm_input_ndens1.data[dens_ind]);
              if (!forces) {
                local2globaldens.push_back(dens_ind);
              }
              local_dens_ind++;
            }
          }
        }
      }
    }
    // Pad the input arrays so the next term type has an aligned offset
    for (j = term_type_counts[current_term_type]; j < COALESCED_DIMENSION(term_type_counts[current_term_type]); j++) {
      func_code.push_back(func_code[term_type_offsets[current_term_type]]); // Use the first code from this term type
      local_dens.push_back(local_dens[term_type_offsets[current_term_type]]);
    }
    for (j = dens_counts[current_term_type]; j < COALESCED_DIMENSION(dens_counts[current_term_type]); j++) {
      dens_values.push_back(dens_values[dens_offsets[current_term_type]]);
      if (!forces) {
        local2globaldens.push_back(local2globaldens[dens_offsets[current_term_type]]);
      }
    }
  }
  check.pause();

  //std::cout << "[G2G_QMMM] Number of threads: " << num_terms << std::endl;
  //std::cout << "[G2G_QMMM] Total Gaussian pairs: " << total_num_terms << std::endl;
  //std::cout << "[G2G_QMMM] Number of significant density elements: " << num_dens_terms << std::endl;
  //std::cout << "[G2G_QMMM] Total density elements: " << total_dens_terms << std::endl;

  prep.start();

  // Pad the input so that out-of-range threads do a dummy calculation (same as the first thread), rather than branching and idling
  for (i = 0; i < QMMM_BLOCK_SIZE - (COALESCED_DIMENSION(term_type_counts[NUM_TERM_TYPES-1]) % QMMM_BLOCK_SIZE); i++) {
    func_code.push_back(func_code[term_type_offsets[NUM_TERM_TYPES-1]]);
    local_dens.push_back(local_dens[term_type_offsets[NUM_TERM_TYPES-1]]);
  }
  for (i = 0; i < QMMM_BLOCK_SIZE - (COALESCED_DIMENSION(dens_counts[NUM_TERM_TYPES-1]) % QMMM_BLOCK_SIZE); i++) {
    dens_values.push_back(dens_values[dens_offsets[NUM_TERM_TYPES-1]]);
  }

  HostMatrix<vec_type<scalar_type, 2> > factor_ac_cpu(COALESCED_DIMENSION(fortran_vars.m), MAX_CONTRACTIONS);
  HostMatrixUInt nuc_cpu(fortran_vars.m, 1);
  CudaMatrix<vec_type<scalar_type, 2> > factor_ac_gpu;
  CudaMatrixUInt nuc_gpu;

  //
  // Set up device arrays for function values and mapping function -> nuclei
  //
  // TODO: tried contracting the nuc / function value arrays to a size of total_funcs rather than m, seems to slow down the kernel
  // Doing this requires some more indexing math in the kernel, but need to test on bigger test case
  uint localfunc = 0, func = 0;
  while (func < fortran_vars.m) {
    nuc_cpu(localfunc) = fortran_vars.nucleii(func) - 1;
    for (uint k = 0; k < fortran_vars.contractions(func); k++) {
      factor_ac_cpu(localfunc, k) = vec_type<scalar_type, 2>(fortran_vars.a_values(func, k), fortran_vars.c_values(func, k));
    }
    func++;
    localfunc++;
  }
  factor_ac_gpu = factor_ac_cpu;
  nuc_gpu = nuc_cpu;

  //
  // Send the MM atom positions and charges to the device
  //
  CudaMatrix<vec_type<scalar_type, 3> > clatom_pos_gpu;
  CudaMatrix<scalar_type> clatom_chg_gpu;
  {
    HostMatrix<vec_type<scalar_type,3> > clatom_pos_cpu(fortran_vars.clatoms,1);
    HostMatrix<scalar_type> clatom_chg_cpu(fortran_vars.clatoms,1);
    for (i = 0; i < fortran_vars.clatoms; i++) {
      clatom_pos_cpu(i) = vec_type<scalar_type,3>(fortran_vars.clatom_positions(i).x,fortran_vars.clatom_positions(i).y,fortran_vars.clatom_positions(i).z);
      clatom_chg_cpu(i) = fortran_vars.clatom_charges(i);
    }
    clatom_pos_gpu = clatom_pos_cpu;
    clatom_chg_gpu = clatom_chg_cpu;
  }

  //
  // Send input arrays (thread->primitive map and thread->density map) to the device
  //
  CudaMatrixUInt dev_func_code(func_code), dev_local_dens(local_dens);
  //
  // Send reduced density matrix to the device
  //
  CudaMatrix<scalar_type> dev_dens_values(dens_values);

  //
  // Allocate output arrays on device
  // Currently, each block in the kernel reduces its output, so the output arrays have length (# block)
  //
  uint partial_out_size = 0, max_partial_size = 0;
  uint out_offsets[NUM_TERM_TYPES];
  out_offsets[0] = 0;
  // Output arrays (probably) don't need to be padded for alignment as only one (or three) threads per block write to them
  for (i = 0; i < NUM_TERM_TYPES; i++) {
    uint this_count = divUp(term_type_counts[i],QMMM_BLOCK_SIZE);
    if (this_count > max_partial_size) max_partial_size = this_count;
    partial_out_size += this_count;
    if (i+1<NUM_TERM_TYPES) { out_offsets[i+1] = partial_out_size; }
  }
  CudaMatrix<vec_type<scalar_type,3> > gpu_partial_mm_forces, gpu_partial_qm_forces;
  CudaMatrix<scalar_type> gpu_partial_fock;
  //
  // Forces: output is partial QM and MM forces
  //
  if (forces)
  {
    gpu_partial_mm_forces.resize(COALESCED_DIMENSION(partial_out_size), fortran_vars.clatoms);
    gpu_partial_qm_forces.resize(COALESCED_DIMENSION(partial_out_size), fortran_vars.atoms);
  //
  // Fock: ouptut is partial Fock elements
  //
  }
  else {
    // The partial Fock matrix is partitioned by term type, so the second (partial) dimension needs to be as big as largest count of a single term type
    gpu_partial_fock.resize(dens_values.size(),max_partial_size);
    dim3 threads(dens_values.size(),max_partial_size);
    dim3 blockSize(32,4);
    dim3 gridSize = divUp(threads,blockSize);
    //
    // Zero the partial Fock matrix on the GPU
    //
    zero_fock<scalar_type><<<gridSize,blockSize>>>(gpu_partial_fock.data,dens_values.size(),max_partial_size);
  }

  //
  // When calculating energies, the energy gets reduced per-block; we figure out the offets/counts of different term types into the partial output energy array here
  //
  uint energies_offsets[NUM_TERM_TYPES];
  uint energies_size = 0;
  CudaMatrix<scalar_type> gpu_qmmm_partial_energies;
  if (!forces) {
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      energies_offsets[i] = energies_size;
      energies_size += divUp(dens_counts[i],QMMM_REDUCE_BLOCK_SIZE);
    }
    gpu_qmmm_partial_energies.resize(energies_size,1);
  }

  //
  // The STR table for F(m,U) calculation is being accessed via texture fetches
  //
  cudaBindTextureToArray(qmmm_str_tex,gammaArrays[devnum]);
  prep.pause_and_sync();

  //
  // Forces kernel
  //
  if (forces) {

    kernel.start();
#define qmmm_forces_parameters \
  term_type_counts[i], factor_ac_gpu.data, nuc_gpu.data, dev_dens_values.data+dens_offset, dev_func_code.data+offset,dev_local_dens.data+offset, \
  gpu_partial_mm_forces.data+force_offset, gpu_partial_qm_forces.data+force_offset, COALESCED_DIMENSION(partial_out_size),clatom_pos_gpu.data,clatom_chg_gpu.data
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_TERM_TYPES];
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (i = 0; i < NUM_TERM_TYPES; i++)
    {
      uint offset = term_type_offsets[i];
      uint dens_offset = dens_offsets[i];
      uint force_offset = out_offsets[i];
      dim3 threads = term_type_counts[i];
      dim3 blockSize(QMMM_BLOCK_SIZE);
      dim3 gridSize = divUp(threads, blockSize);
      switch (i) {
        case 0: gpu_qmmm_forces<scalar_type,0><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 1: gpu_qmmm_forces<scalar_type,1><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 2: gpu_qmmm_forces<scalar_type,2><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 3: gpu_qmmm_forces<scalar_type,3><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 4: gpu_qmmm_forces<scalar_type,4><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
        case 5: gpu_qmmm_forces<scalar_type,5><<<gridSize,blockSize,0,stream[i]>>>( qmmm_forces_parameters ); break;
      }
    }
    cudaDeviceSynchronize();
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }
    kernel.pause();
  }
  //
  // Energy/Fock kernel
  //
  else {

    kernel.start();

#define qmmm_fock_parameters \
  term_type_counts[i], factor_ac_gpu.data, nuc_gpu.data, dev_func_code.data+offset,dev_local_dens.data+offset, /*num_dens_terms,*/ \
  gpu_partial_fock.data+fock_offset, dens_values.size()/*COALESCED_DIMENSION(num_dens_terms)*/,clatom_pos_gpu.data,clatom_chg_gpu.data//,fock_out_offset
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_TERM_TYPES];
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (i = 0; i < NUM_TERM_TYPES; i++)
    {
      uint offset = term_type_offsets[i];
      uint fock_offset = dens_offsets[i];
      dim3 threads = term_type_counts[i];
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
      dim3 reduceThreads = dens_counts[i];
      dim3 reduceBlockSize(QMMM_REDUCE_BLOCK_SIZE);
      dim3 reduceGridSize = divUp(reduceThreads,reduceBlockSize);
      gpu_qmmm_fock_reduce<scalar_type><<<reduceGridSize,reduceBlockSize,0,stream[i]>>>( gpu_partial_fock.data+fock_offset, dev_dens_values.data+fock_offset, gpu_qmmm_partial_energies.data+energies_offsets[i],
                                                                                         dens_values.size(), max_partial_size, dens_counts[i] );
    }
    cudaDeviceSynchronize();
    for (i = 0; i < NUM_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }
    kernel.pause();
  }

  if (forces) {
    //
    // Download the partial forces from the device
    // TODO: this could maybe be done asynchronously with the kernels; as one term type finishes we can download its forces, etc
    //
    down.start();
    HostMatrix<vec_type<scalar_type,3> > cpu_partial_mm_forces(gpu_partial_mm_forces), cpu_partial_qm_forces(gpu_partial_qm_forces);
    down.pause_and_sync();

    //
    // Accumulate partial output
    //
    // TODO: need to think about how to accumulate individual force terms
    // Currently, we reduce on a per-block basis in the kernel, then accumulate the block results here on the host
    //
    // The energy partial results are being reduced on-device and that works very well, could probably do that for forces too
    //
    reduce.start();
    for (i = 0; i < fortran_vars.atoms; i++) {
      for (j = 0; j < partial_out_size; j++) {
        qm_forces[i + 0 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).x;
        qm_forces[i + 1 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).y;
        qm_forces[i + 2 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).z;
      }
    }
    for (i = 0; i < fortran_vars.clatoms; i++) {
      for (j = 0; j < partial_out_size; j++) {
        mm_forces[i + 0 * fortran_vars.clatoms] += cpu_partial_mm_forces(j,i).x;
        mm_forces[i + 1 * fortran_vars.clatoms] += cpu_partial_mm_forces(j,i).y;
        mm_forces[i + 2 * fortran_vars.clatoms] += cpu_partial_mm_forces(j,i).z;
      }
    }
    reduce.pause();
  } else {
    //
    // Download reduced Fock matrix and partially reduced energies
    // The Fock matrix has been reduced to the first row of the output, so we only want that much of the device array
    //
    down.start();
    HostMatrix<scalar_type> cpu_fock(dens_values.size());
    cudaMemcpy(cpu_fock.data,gpu_partial_fock.data,cpu_fock.bytes(),cudaMemcpyDeviceToHost);
    HostMatrix<scalar_type> cpu_partial_energies(gpu_qmmm_partial_energies);
    down.pause_and_sync();

    //
    // Send Fock elements back to RMM(M11) and do final reduction of e-nuc energies into Es
    //
    reduce.start();
    for (uint t = 0; t < NUM_TERM_TYPES; t++) {
      for (i = dens_offsets[t]; i < dens_offsets[t] + dens_counts[t]; i++) {
        uint dens_ind = local2globaldens[i];
        fortran_vars.rmm_1e_output(dens_ind) += cpu_fock(i);
      }
    }
    for (i = 0; i < energies_size; i++) {
      Es += cpu_partial_energies(i);
    }
    reduce.pause();
  }

  //cout << "[G2G_QMMM] nuc-nuc: " << nuc << " overlap check: " << check << " kernel prep: " << prep << endl;
  //cout << "[G2G_QMMM] kernel: " << kernel << " download: " << down << " host reduction: " << reduce << endl;

  cudaUnbindTexture(qmmm_str_tex);

  cudaAssertNoError("qmmm");
}

template<class scalar_type>
void clean_gamma( void ) {
  int previous_device; cudaGetDevice(&previous_device);
  int gpu_devices = cudaGetGPUCount();
  for(int i = 0; i < gpu_devices; i++) {
    if(cudaSetDevice(i) != cudaSuccess)  std::cout << "Error: can't set the device " << i << std::endl;
    cudaFreeArray(gammaArrays[i]);
  }
  delete[] gammaArrays;
  cudaSetDevice(previous_device);
  cudaAssertNoError("clean_gamma");
}
#if FULL_DOUBLE
template void g2g_qmmm<double,true>(double* qm_forces, double* mm_forces, double& Ens, double& Es);
template void g2g_qmmm<double,false>(double* qm_forces, double* mm_forces, double& Ens, double& Es);
template void clean_gamma<double>( void );
#else
template void g2g_qmmm<float,true>(double* qm_forces, double* mm_forces, double& Ens, double& Es);
template void g2g_qmmm<float,false>(double* qm_forces, double* mm_forces, double& Ens, double& Es);
template void clean_gamma<float>( void );
#endif

#define NUM_COULOMB_TERM_TYPES 6 // 6 types when using s,p,and d functions: s-s,p-s,p-p,d-s,d-p,d-d

//
// Main QM/MM routine
// If forces = true, calculate gradients of QM/MM operator and return in qm_forces and mm_forces
// If forces = false, calculate Fock matrix elements of QM/MM operator (returned in RMM(M11) back in Fortran)
//                    and QM/MM energy of the current density (nuc-nuc in Ens, e-nuc in Es)
//
template <class scalar_type,bool forces> void g2g_coulomb(double* qm_forces, double& Es)
{
  int devnum = -1;
  if (cudaGetDevice(&devnum) != cudaSuccess) throw std::runtime_error("Could not get device number! (Coulomb)");
  uint i,j,ni,nj;
  uint i_orbitals, j_orbitals;
  uint nuc_i,nuc_j;
  vec_type<double,3> A,B,AmB;
  double ai,aj;
  double dsq,ksi,zeta;
  uint num_terms=0, total_num_terms = 0;
  std::vector<uint> func_code, local_dens;

  //std::vector<uint> local2globaldens; // Maps reduced density/Fock -> full density/Fock
  std::vector<scalar_type> dens_values;

  uint term_type_counts[NUM_COULOMB_TERM_TYPES]; // Number of threads for a particular type of term (0 = s-s, 1 = p-s, etc)
  uint term_type_offsets[NUM_COULOMB_TERM_TYPES]; // Offsets into the input arrays for each term type
  term_type_offsets[0] = 0; // s-s starts at 0
  uint i_begin, i_end, j_begin, j_end;
  uint tmp_ind = 0;
  uint s_start = 0, p_start = fortran_vars.s_funcs, d_start = fortran_vars.s_funcs + fortran_vars.p_funcs*3, m = fortran_vars.m;
  uint dd_orb;
  if (forces) { dd_orb = 6; } // 1 for d-d means 6 threads use one of dxx, dyx, etc for func i
  else        { dd_orb = 6; }

  //                                       s-s        p-s        p-p        d-s        d-p        d-d
  uint i_begin_vals[MAX_TERM_TYPE]   = { s_start,   p_start,   p_start,   d_start,   d_start,   d_start};
  uint i_end_vals[MAX_TERM_TYPE]     = { p_start,   d_start,   d_start,   m,         m,         m      };
  uint j_begin_vals[MAX_TERM_TYPE]   = { s_start,   s_start,   p_start,   s_start,   p_start,   d_start};
  uint j_end_vals[MAX_TERM_TYPE]     = { p_start-1, p_start-1, d_start-1, p_start-1, d_start-1, m-1    };
  uint i_orbital_vals[MAX_TERM_TYPE] = { 1,         3,         3,         6,         6,         dd_orb };
  uint j_orbital_vals[MAX_TERM_TYPE] = { 1,         1,         3,         1,         3,         6      };

  uint local_dens_ind, num_dens_terms = 0, total_dens_terms = 0;
  uint dens_counts[NUM_COULOMB_TERM_TYPES], dens_offsets[NUM_COULOMB_TERM_TYPES];
  dens_offsets[0] = 0;
  uint tmp_dens_ind = 0;

  Timer check,prep,kernel,down,reduce;

  //
  // Do check between all basis primitives to find those with significant overlap
  // Check the resulting Gaussian argument from two primitives to the rmax parameter; only use primitives within that cut-off
  //
  // A single thread gets mapped to a pair of significant primitives
  // We set up here arrays that tell which two functions/two primitives a thread is calculating
  // We also pick out the density matrix elements for significant functions here
  //
  check.start();
  for (uint current_term_type = 0; current_term_type < NUM_COULOMB_TERM_TYPES; current_term_type++) {

    term_type_counts[current_term_type] = 0;
    i_begin = i_begin_vals[current_term_type]; i_end = i_end_vals[current_term_type];
    j_begin = j_begin_vals[current_term_type]; j_end = j_end_vals[current_term_type];
    i_orbitals = i_orbital_vals[current_term_type];
    j_orbitals = j_orbital_vals[current_term_type];

    dens_counts[current_term_type] = 0;
    local_dens_ind = 0;

    // We pad the input arrays between term types, so the offsets for each term type need to be tracked
    if (current_term_type > 0) {
      tmp_ind += COALESCED_DIMENSION(term_type_counts[current_term_type-1]);
      term_type_offsets[current_term_type] = tmp_ind;
      tmp_dens_ind += COALESCED_DIMENSION(dens_counts[current_term_type-1]);
      dens_offsets[current_term_type] = tmp_dens_ind;
    }

    // function i, center A
    for (i = i_begin; i < i_end; i += i_orbitals) {
      nuc_i = fortran_vars.nucleii(i) - 1;
      A = fortran_vars.atom_positions(nuc_i);
      // function j, center B
      for (j = j_begin; j <= ((i > j_end)? j_end : i); j += j_orbitals) {
        nuc_j = fortran_vars.nucleii(j) - 1;
        B = fortran_vars.atom_positions(nuc_j);
        AmB = A - B;
        dsq = length2(AmB);
        bool use_funcs = false; // Do these two functions have any significant primitive pairs?
        // primitive ni, function i
        for (ni = 0; ni < fortran_vars.contractions(i); ni++) {
          // primitive nj, function j
          for (nj = 0; nj < fortran_vars.contractions(j); nj++) {
            ai = fortran_vars.a_values(i,ni);
            aj = fortran_vars.a_values(j,nj);
            zeta = ai + aj;
            ksi = ai * aj / zeta;
            total_num_terms += (i==j)? i_orbitals*(i_orbitals+1)/2 : i_orbitals * j_orbitals;

            if (dsq*ksi < fortran_vars.rmax) {
              use_funcs = true;
              num_terms++;
              term_type_counts[current_term_type]++;

              // Encoding which two primitives this thread will calculate a force term for in one number
              // NOTE: with a long integer, we can use this scheme up to an m of about 9000
              //       if m needs to ever be larger than that, we need to break this up into multiple arrays
              uint this_func_code = nj;                                                         // First, primitive # nj in the lowest part
              this_func_code     += ni * MAX_CONTRACTIONS;                                      // Primitive # ni after the space for nj
              this_func_code     += j  * MAX_CONTRACTIONS * MAX_CONTRACTIONS;                   // Function # j after space for primitives
              this_func_code     += i  * MAX_CONTRACTIONS * MAX_CONTRACTIONS * fortran_vars.m;  // Finally, function # i in the highest part

              func_code.push_back(this_func_code); // Which primitives the thread represents
              local_dens.push_back(local_dens_ind); // Which part of the (reduced) density matrix the thread needs

            }
          }
        }

        total_dens_terms += (i==j)? i_orbitals*(i_orbitals+1)/2 : i_orbitals * j_orbitals;
        // dens_values is a reduced density matrix that only keeps the elements of functions with significant primitive pairs
        // local2globaldens maps from a reduced density/Fock index back to the full matrix
        if (use_funcs) {
          for (uint i_orbital = 0; i_orbital < i_orbitals; i_orbital++) {
            uint j_orbital_finish = (i==j)? i_orbital+1 : j_orbitals;
            for (uint j_orbital = 0; j_orbital < j_orbital_finish; j_orbital++) {
              num_dens_terms++;
              dens_counts[current_term_type]++;

              uint dens_ind = (i+i_orbital) + (2*fortran_vars.m-((j+j_orbital)+1))*(j+j_orbital)/2;
              dens_values.push_back(fortran_vars.rmm_input_ndens1.data[dens_ind]);
              if (!forces) {
              //  local2globaldens.push_back(dens_ind);
              }
              local_dens_ind++;
            }
          }
        }
      }
    }
    // Pad the input arrays so the next term type has an aligned offset
    for (j = term_type_counts[current_term_type]; j < COALESCED_DIMENSION(term_type_counts[current_term_type]); j++) {
      func_code.push_back(func_code[term_type_offsets[current_term_type]]); // Use the first code from this term type
      local_dens.push_back(local_dens[term_type_offsets[current_term_type]]);
    }
    for (j = dens_counts[current_term_type]; j < COALESCED_DIMENSION(dens_counts[current_term_type]); j++) {
      dens_values.push_back(dens_values[dens_offsets[current_term_type]]);
      if (!forces) {
      //  local2globaldens.push_back(local2globaldens[dens_offsets[current_term_type]]);
      }
    }
  }
  check.pause();

  //std::cout << "[G2G_COULOMB] Number of threads: " << num_terms << std::endl;
  //std::cout << "[G2G_COULOMB] Total Gaussian pairs: " << total_num_terms << std::endl;
  //std::cout << "[G2G_COULOMB] Number of significant density elements: " << num_dens_terms << std::endl;
  //std::cout << "[G2G_COULOMB] Total density elements: " << total_dens_terms << std::endl;

  prep.start();

  // Pad the input so that out-of-range threads do a dummy calculation (same as the first thread), rather than branching and idling
  for (i = 0; i < QMMM_BLOCK_SIZE - (COALESCED_DIMENSION(term_type_counts[NUM_COULOMB_TERM_TYPES-1]) % QMMM_BLOCK_SIZE); i++) {
    func_code.push_back(func_code[term_type_offsets[NUM_COULOMB_TERM_TYPES-1]]);
    local_dens.push_back(local_dens[term_type_offsets[NUM_COULOMB_TERM_TYPES-1]]);
  }
  for (i = 0; i < QMMM_BLOCK_SIZE - (COALESCED_DIMENSION(dens_counts[NUM_COULOMB_TERM_TYPES-1]) % QMMM_BLOCK_SIZE); i++) {
    dens_values.push_back(dens_values[dens_offsets[NUM_COULOMB_TERM_TYPES-1]]);
  }

  HostMatrix<vec_type<scalar_type, 2> > factor_ac_cpu(COALESCED_DIMENSION(fortran_vars.m), MAX_CONTRACTIONS);
  HostMatrixUInt nuc_cpu(fortran_vars.m, 1);
  CudaMatrix<vec_type<scalar_type, 2> > factor_ac_gpu;
  CudaMatrixUInt nuc_gpu;

  //
  // Set up device arrays for function values and mapping function -> nuclei
  //
  // TODO: tried contracting the nuc / function value arrays to a size of total_funcs rather than m, seems to slow down the kernel
  // Doing this requires some more indexing math in the kernel, but need to test on bigger test case
  uint localfunc = 0, func = 0;
  while (func < fortran_vars.m) {
    nuc_cpu(localfunc) = fortran_vars.nucleii(func) - 1;
    for (uint k = 0; k < fortran_vars.contractions(func); k++) {
      factor_ac_cpu(localfunc, k) = vec_type<scalar_type, 2>(fortran_vars.a_values(func, k), fortran_vars.c_values(func, k));
    }
    func++;
    localfunc++;
  }
  factor_ac_gpu = factor_ac_cpu;
  nuc_gpu = nuc_cpu;

  std::vector<vec_type<scalar_type,2> > factor_ac_dens_cpu;
  std::vector<scalar_type> fit_dens_cpu;
  std::vector<vec_type<scalar_type,3> > nuc_dens_cpu;
  std::vector<uint> nuc_ind_dens_cpu;
  uint input_size = 0;
  for (func = 0; func < fortran_vars.s_funcs_dens; func++) {
    uint nuc_ind = fortran_vars.nucleii_dens(func) - 1;
    double3 nuc_pos = fortran_vars.atom_positions(nuc_ind);
    for (uint k = 0; k < fortran_vars.contractions_dens(func); k++) {
      factor_ac_dens_cpu.push_back(vec_type<scalar_type,2>(fortran_vars.a_values_dens(func,k),fortran_vars.c_values_dens(func,k)));
      nuc_dens_cpu.push_back(vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
      nuc_ind_dens_cpu.push_back(nuc_ind);
      fit_dens_cpu.push_back(fortran_vars.af_input_ndens1(func));
      input_size++;
    }
  }
  uint s_end = input_size;
  uint tmp_size = input_size;
  for (func = tmp_size; func < COALESCED_DIMENSION(tmp_size); func++) {
    factor_ac_dens_cpu.push_back(factor_ac_dens_cpu[tmp_size-1]);
    nuc_dens_cpu.push_back(nuc_dens_cpu[tmp_size-1]);
    nuc_ind_dens_cpu.push_back(nuc_ind_dens_cpu[tmp_size-1]);
    fit_dens_cpu.push_back(fit_dens_cpu[tmp_size-1]);
    input_size++;
  }
  uint p_offset = input_size;
  tmp_size = 0;
  for (func = fortran_vars.s_funcs_dens; func < fortran_vars.s_funcs_dens+fortran_vars.p_funcs_dens*3; func += 3) {
    uint nuc_ind = fortran_vars.nucleii_dens(func) - 1;
    double3 nuc_pos = fortran_vars.atom_positions(nuc_ind);
    for (uint k = 0; k < fortran_vars.contractions_dens(func); k++) {
      for (uint f = 0; f < 3; f++) {
        factor_ac_dens_cpu.push_back(vec_type<scalar_type,2>(fortran_vars.a_values_dens(func,k),fortran_vars.c_values_dens(func,k)));
        nuc_dens_cpu.push_back(vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
        nuc_ind_dens_cpu.push_back(nuc_ind);
        fit_dens_cpu.push_back(fortran_vars.af_input_ndens1(func+f));
      }
      tmp_size += 3;
      input_size += 3;
      if (tmp_size == 126) {
        for (uint f = 0; f < 2; f++) factor_ac_dens_cpu.push_back(vec_type<scalar_type,2>(fortran_vars.a_values_dens(func,k),fortran_vars.c_values_dens(func,k)));
        for (uint f = 0; f < 2; f++) nuc_dens_cpu.push_back(vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
        for (uint f = 0; f < 2; f++) nuc_ind_dens_cpu.push_back(nuc_ind);
        for (uint f = 0; f < 2; f++) fit_dens_cpu.push_back(fortran_vars.af_input_ndens1(func+f));
        tmp_size = 0;
        input_size += 2;
      }
    }
  }
  uint p_end = input_size;
  tmp_size = input_size;
  for (func = tmp_size; func < COALESCED_DIMENSION(tmp_size); func++) {
    factor_ac_dens_cpu.push_back(factor_ac_dens_cpu[tmp_size-1]);
    nuc_dens_cpu.push_back(nuc_dens_cpu[tmp_size-1]);
    nuc_ind_dens_cpu.push_back(nuc_ind_dens_cpu[tmp_size-1]);
    fit_dens_cpu.push_back(fit_dens_cpu[tmp_size-1]);
    input_size++;
  }
  uint d_offset = input_size;
  tmp_size = 0;
  for (func = fortran_vars.s_funcs_dens+fortran_vars.p_funcs_dens*3; func < fortran_vars.m_dens; func += 6) {
    uint nuc_ind = fortran_vars.nucleii_dens(func) - 1;
    double3 nuc_pos = fortran_vars.atom_positions(nuc_ind);
    for (uint k = 0; k < fortran_vars.contractions_dens(func); k++) {
      for (uint f = 0; f < 6; f++) {
        factor_ac_dens_cpu.push_back(vec_type<scalar_type,2>(fortran_vars.a_values_dens(func,k),fortran_vars.c_values_dens(func,k)));
        nuc_dens_cpu.push_back(vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
        nuc_ind_dens_cpu.push_back(nuc_ind);
        fit_dens_cpu.push_back(fortran_vars.af_input_ndens1(func+f));
      }
      tmp_size += 6;
      input_size += 6;
      if (tmp_size == 126) {
        for (uint f = 0; f < 2; f++) factor_ac_dens_cpu.push_back(vec_type<scalar_type,2>(fortran_vars.a_values_dens(func,k),fortran_vars.c_values_dens(func,k)));
        for (uint f = 0; f < 2; f++) nuc_dens_cpu.push_back(vec_type<scalar_type, 3>(nuc_pos.x, nuc_pos.y, nuc_pos.z));
        for (uint f = 0; f < 2; f++) nuc_ind_dens_cpu.push_back(nuc_ind);
        for (uint f = 0; f < 2; f++) fit_dens_cpu.push_back(fortran_vars.af_input_ndens1(func+f));
        tmp_size = 0;
        input_size += 2;
      }
    }
  }
  uint d_end = input_size;
  CudaMatrix<vec_type<scalar_type,2> > factor_ac_dens_gpu(factor_ac_dens_cpu);
  CudaMatrix<vec_type<scalar_type, 3> > nuc_dens_gpu(nuc_dens_cpu);
  CudaMatrix<scalar_type> fit_dens_gpu(fit_dens_cpu);
  CudaMatrix<uint> nuc_ind_dens_gpu(nuc_ind_dens_cpu);

  //
  // Send input arrays (thread->primitive map and thread->density map) to the device
  //
  CudaMatrixUInt dev_func_code(func_code), dev_local_dens(local_dens);
  //
  // Send reduced density matrix to the device
  //
  CudaMatrix<scalar_type> dev_dens_values(dens_values);

  //
  // Allocate output arrays on device
  // Currently, each block in the kernel reduces its output, so the output arrays have length (# block)
  //
  uint partial_out_size = 0;//, max_partial_size = 0;
  uint out_offsets[NUM_COULOMB_TERM_TYPES];
  // Output arrays (probably) don't need to be padded for alignment as only one (or three) threads per block write to them
  for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
    uint this_count = divUp(term_type_counts[i],QMMM_BLOCK_SIZE);
    //if (this_count > max_partial_size) max_partial_size = this_count;
    out_offsets[i] = partial_out_size;
    partial_out_size += this_count;
  }
  CudaMatrix<vec_type<scalar_type,3> > gpu_partial_qm_forces;
  //CudaMatrix<scalar_type> gpu_partial_fock;
  //
  // Forces: output is partial QM and MM forces
  //
  if (forces)
  {
    gpu_partial_qm_forces.resize(COALESCED_DIMENSION(partial_out_size), fortran_vars.atoms);
    dim3 threads(COALESCED_DIMENSION(partial_out_size),fortran_vars.atoms);
    dim3 blockSize(32,4);
    dim3 gridSize = divUp(threads,blockSize);
    zero_forces<scalar_type><<<gridSize,blockSize>>>(gpu_partial_qm_forces.data,COALESCED_DIMENSION(partial_out_size),fortran_vars.atoms);
    //cudaMemset(gpu_partial_qm_forces.data, 0, COALESCED_DIMENSION(partial_out_size) * fortran_vars.atoms * sizeof(vec_type<scalar_type,3>));

  //
  // Fock: ouptut is partial Fock elements
  //
  }
  /*else {
    // The partial Fock matrix is partitioned by term type, so the second (partial) dimension needs to be as big as largest count of a single term type
    gpu_partial_fock.resize(dens_values.size(),max_partial_size);
    dim3 threads(dens_values.size(),max_partial_size);
    dim3 blockSize(32,4);
    dim3 gridSize = divUp(threads,blockSize);
    //
    // Zero the partial Fock matrix on the GPU
    //
    zero_fock<scalar_type><<<gridSize,blockSize>>>(gpu_partial_fock.data,dens_values.size(),max_partial_size);
  }*/

  //
  // When calculating energies, the energy gets reduced per-block; we figure out the offets/counts of different term types into the partial output energy array here
  //
  /*uint energies_offsets[NUM_COULOMB_TERM_TYPES];
  uint energies_size = 0;
  CudaMatrix<scalar_type> gpu_qmmm_partial_energies;
  if (!forces) {
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      energies_offsets[i] = energies_size;
      energies_size += divUp(dens_counts[i],QMMM_REDUCE_BLOCK_SIZE);
    }
    gpu_qmmm_partial_energies.resize(energies_size,1);
  }*/

  //
  // The STR table for F(m,U) calculation is being accessed via texture fetches
  //
  cudaBindTextureToArray(qmmm_str_tex,gammaArrays[devnum]);
  prep.pause_and_sync();

  //
  // Forces kernel
  //
  if (forces) {

    kernel.start();
#define coulomb_forces_parameters \
  term_type_counts[i], factor_ac_gpu.data, nuc_gpu.data, dev_dens_values.data+dens_offset, dev_func_code.data+offset,dev_local_dens.data+offset, \
  gpu_partial_qm_forces.data+force_offset, COALESCED_DIMENSION(partial_out_size),factor_ac_dens_gpu.data,nuc_dens_gpu.data,nuc_ind_dens_gpu.data,fit_dens_gpu.data, \
  s_end, p_end, d_end, p_offset, d_offset
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_COULOMB_TERM_TYPES];
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++)
    {
      uint offset = term_type_offsets[i];
      uint dens_offset = dens_offsets[i];
      uint force_offset = out_offsets[i];
      dim3 threads = term_type_counts[i];
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
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }
    kernel.pause();
  }
  //
  // Energy/Fock kernel
  //
  /*else {

    kernel.start();

#define qmmm_fock_parameters \
  term_type_counts[i], factor_ac_gpu.data, nuc_gpu.data, dev_func_code.data+offset,dev_local_dens.data+offset, \
  gpu_partial_fock.data+fock_offset, dens_values.size(),clatom_pos_gpu.data,clatom_chg_gpu.data//,fock_out_offset
    // Each term type is calculated asynchronously
    cudaStream_t stream[NUM_COULOMB_TERM_TYPES];
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      cudaStreamCreate(&stream[i]);
    }
    //
    // Begin launching kernels (one for each type of term, 0 = s-s, 1 = p-s, etc)
    //
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++)
    {
      uint offset = term_type_offsets[i];
      uint fock_offset = dens_offsets[i];
      dim3 threads = term_type_counts[i];
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
      dim3 reduceThreads = dens_counts[i];
      dim3 reduceBlockSize(QMMM_REDUCE_BLOCK_SIZE);
      dim3 reduceGridSize = divUp(reduceThreads,reduceBlockSize);
      gpu_qmmm_fock_reduce<scalar_type><<<reduceGridSize,reduceBlockSize,0,stream[i]>>>( gpu_partial_fock.data+fock_offset, dev_dens_values.data+fock_offset, gpu_qmmm_partial_energies.data+energies_offsets[i],
                                                                                         dens_values.size(), max_partial_size, dens_counts[i] );
    }
    cudaDeviceSynchronize();
    for (i = 0; i < NUM_COULOMB_TERM_TYPES; i++) {
      cudaStreamDestroy(stream[i]);
    }
    kernel.pause();
  }*/

  if (forces) {
    //
    // Download the partial forces from the device
    // TODO: this could maybe be done asynchronously with the kernels; as one term type finishes we can download its forces, etc
    //
    down.start();
    HostMatrix<vec_type<scalar_type,3> > cpu_partial_qm_forces(gpu_partial_qm_forces);
    down.pause_and_sync();

    //
    // Accumulate partial output
    //
    // TODO: need to think about how to accumulate individual force terms
    // Currently, we reduce on a per-block basis in the kernel, then accumulate the block results here on the host
    //
    // The energy partial results are being reduced on-device and that works very well, could probably do that for forces too
    //
    reduce.start();
    for (i = 0; i < fortran_vars.atoms; i++) {
      for (j = 0; j < partial_out_size; j++) {
        qm_forces[i + 0 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).x;
        qm_forces[i + 1 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).y;
        qm_forces[i + 2 * fortran_vars.atoms] += cpu_partial_qm_forces(j,i).z;
      }
    }
    reduce.pause();
  } /*else {
    //
    // Download reduced Fock matrix and partially reduced energies
    // The Fock matrix has been reduced to the first row of the output, so we only want that much of the device array
    //
    down.start();
    HostMatrix<scalar_type> cpu_fock(dens_values.size());
    cudaMemcpy(cpu_fock.data,gpu_partial_fock.data,cpu_fock.bytes(),cudaMemcpyDeviceToHost);
    HostMatrix<scalar_type> cpu_partial_energies(gpu_qmmm_partial_energies);
    down.pause_and_sync();

    //
    // Send Fock elements back to RMM(M11) and do final reduction of e-nuc energies into Es
    //
    reduce.start();
    for (uint t = 0; t < NUM_COULOMB_TERM_TYPES; t++) {
      for (i = dens_offsets[t]; i < dens_offsets[t] + dens_counts[t]; i++) {
        uint dens_ind = local2globaldens[i];
        fortran_vars.rmm_1e_output(dens_ind) += cpu_fock(i);
      }
    }
    for (i = 0; i < energies_size; i++) {
      Es += cpu_partial_energies(i);
    }
    reduce.pause();
  }*/

  //cout << "[G2G_COULOMB] overlap check: " << check << " kernel prep: " << prep << endl;
  //cout << "[G2G_COULOMB] kernel: " << kernel << " download: " << down << " host reduction: " << reduce << endl;

  cudaUnbindTexture(qmmm_str_tex);

  cudaAssertNoError("qmmm");
}

#if FULL_DOUBLE
template void g2g_coulomb<double,true>(double* qm_forces, double& Es);
template void g2g_coulomb<double,false>(double* qm_forces, double& Es);
#else
template void g2g_coulomb<float,true>(double* qm_forces, double& Es);
template void g2g_coulomb<float,false>(double* qm_forces, double& Es);
#endif

}
