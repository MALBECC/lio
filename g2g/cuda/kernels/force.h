template <class scalar_type>
__global__ void gpu_compute_forces(uint points, scalar_type* force_factors,
                                   vec_type<scalar_type, 4>* density_deriv,
                                   vec_type<scalar_type, 4>* forces,
                                   uint nucleii_count) {
  uint atom = index_x(blockDim, blockIdx, threadIdx);
  bool valid_thread = (atom < nucleii_count);

  vec_type<scalar_type, 4> atom_force(0.0f, 0.0f, 0.0f, 0.0f);

  __shared__ scalar_type factor_sh[FORCE_BLOCK_SIZE];

  for (uint point_base = 0; point_base < points;
       point_base += FORCE_BLOCK_SIZE) {
    __syncthreads();
    if (point_base + threadIdx.x < points)
      factor_sh[threadIdx.x] = force_factors[point_base + threadIdx.x];
    __syncthreads();

    if (valid_thread) {
      for (uint point_sub = 0;
           point_sub < FORCE_BLOCK_SIZE && (point_base + point_sub < points);
           point_sub++) {
        vec_type<scalar_type, 4> density_deriv_local(
            density_deriv[COALESCED_DIMENSION(points) * atom +
                          (point_base + point_sub)]);  // uncoalesced
        atom_force += density_deriv_local * factor_sh[point_sub];
      }
    }
  }

  if (valid_thread) {
    forces[atom] = atom_force;
  }
}

///==========================================================================================================================
template <class scalar_type>
__global__ void gpu_compute_forces_open(
    uint points, scalar_type* force_factors_a, scalar_type* force_factors_b,
    vec_type<scalar_type, 4>* density_deriv_a,
    vec_type<scalar_type, 4>* density_deriv_b,
    vec_type<scalar_type, 4>* forces_a, vec_type<scalar_type, 4>* forces_b,
    uint nucleii_count) {
  uint atom = index_x(blockDim, blockIdx, threadIdx);
  bool valid_thread = (atom < nucleii_count);

  vec_type<scalar_type, 4> atom_force_a(0.0f, 0.0f, 0.0f, 0.0f);
  vec_type<scalar_type, 4> atom_force_b(0.0f, 0.0f, 0.0f, 0.0f);

  __shared__ scalar_type factor_a_sh[FORCE_BLOCK_SIZE];
  __shared__ scalar_type factor_b_sh[FORCE_BLOCK_SIZE];

  for (uint point_base = 0; point_base < points;
       point_base += FORCE_BLOCK_SIZE) {
    __syncthreads();
    if (point_base + threadIdx.x < points) {
      factor_a_sh[threadIdx.x] = force_factors_a[point_base + threadIdx.x];
      factor_b_sh[threadIdx.x] = force_factors_b[point_base + threadIdx.x];
    }
    __syncthreads();

    if (valid_thread) {
      for (uint point_sub = 0;
           point_sub < FORCE_BLOCK_SIZE && (point_base + point_sub < points);
           point_sub++) {
        vec_type<scalar_type, 4> density_deriv_a_local(
            density_deriv_a[COALESCED_DIMENSION(points) * atom +
                            (point_base + point_sub)]);  // uncoalesced
        vec_type<scalar_type, 4> density_deriv_b_local(
            density_deriv_b[COALESCED_DIMENSION(points) * atom +
                            (point_base + point_sub)]);  // uncoalesced
        atom_force_a += density_deriv_a_local * factor_a_sh[point_sub];
        atom_force_b += density_deriv_b_local * factor_b_sh[point_sub];
      }
    }
  }

  if (valid_thread) {
    forces_a[atom] = atom_force_a;
    forces_b[atom] = atom_force_b;
  }
}
