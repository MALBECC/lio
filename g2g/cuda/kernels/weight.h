// TODO: precomputar/precargar las cosas por atomo una sola vez para todos los
// puntos, pasar esto a nucleii_count, etc, coalescing

template <class scalar_type>
__global__ void gpu_compute_weights(uint points,
                                    vec_type<scalar_type, 4>* point_positions,
                                    vec_type<scalar_type, 4>* atom_position_rm,
                                    scalar_type* weights, uint* nucleii,
                                    uint nucleii_count) {
  uint point = index_x(blockDim, blockIdx, threadIdx);

  __shared__ vec_type<scalar_type, 3> atom_position_sh[MAX_ATOMS];
  __shared__ uint nucleii_sh[MAX_ATOMS];
  __shared__ scalar_type rm_sh[MAX_ATOMS];

  for (uint i = 0; i < gpu_atoms; i += WEIGHT_BLOCK_SIZE) {
    if (i + threadIdx.x < gpu_atoms) {
      vec_type<scalar_type, 4> atom_position_rm_local(
          atom_position_rm[i + threadIdx.x]);
      atom_position_sh[i + threadIdx.x] =
          vec_type<scalar_type, 3>(atom_position_rm_local);
      rm_sh[i + threadIdx.x] = atom_position_rm_local.w;
    }
  }

  for (uint i = 0; i < nucleii_count; i += WEIGHT_BLOCK_SIZE) {
    if (i + threadIdx.x < nucleii_count) {
      nucleii_sh[i + threadIdx.x] = nucleii[i + threadIdx.x];
    }
  }

  __syncthreads();

  if (point >= points) return;

  vec_type<scalar_type, 4> point_position4 = point_positions[point];
  vec_type<scalar_type, 3> point_position = vec_type<scalar_type, 3>(
      point_position4.x, point_position4.y, point_position4.z);
  uint atom_of_point = (uint)floor(point_position4.w);

  scalar_type P_total = 0.0f;
  scalar_type P_atom = 0.0f;
  bool found_atom = false;

  for (uint nuc_i = 0; nuc_i < nucleii_count; nuc_i++) {
    uint atom_i = nucleii_sh[nuc_i];
    scalar_type P_curr = 1.0f;

    for (uint nuc_j = 0; nuc_j < nucleii_count; nuc_j++) {
      uint atom_j = nucleii_sh[nuc_j];
      if (atom_i == atom_j) continue;

      scalar_type u =
          (::distance(point_position, atom_position_sh[atom_i]) -
           ::distance(point_position, atom_position_sh[atom_j])) /
          ::distance(atom_position_sh[atom_i], atom_position_sh[atom_j]);

      scalar_type x;
      x = rm_sh[atom_i] / rm_sh[atom_j];
      x = (x - 1.0f) / (x + 1.0f);
      u += (x / (x * x - 1.0f)) * (1.0f - u * u);

#pragma unroll 3
      for (uint i = 0; i < 3; i++) u = 1.5f * u - 0.5f * (u * u * u);

      u = 0.5f * (1.0f - u);

      P_curr *= u;
    }

    if (atom_i == atom_of_point) {
      found_atom = true;
      P_atom = P_curr;
    }
    P_total += P_curr;
  }

  if (!found_atom) {
    P_atom = 1.0f;
    uint atom_i = atom_of_point;

    for (uint nuc_j = 0; nuc_j < nucleii_count; nuc_j++) {
      uint atom_j = nucleii_sh[nuc_j];

      scalar_type u =
          (::distance(point_position, atom_position_sh[atom_i]) -
           ::distance(point_position, atom_position_sh[atom_j])) /
          ::distance(atom_position_sh[atom_i], atom_position_sh[atom_j]);

      scalar_type x;
      x = rm_sh[atom_i] / rm_sh[atom_j];
      x = (x - 1.0f) / (x + 1.0f);
      u += (x / (x * x - 1.0f)) * (1.0f - u * u);

#pragma unroll 3
      for (uint i = 0; i < 3; i++) u = 1.5f * u - 0.5f * (u * u * u);

      u = 0.5f * (1.0f - u);

      P_atom *= u;
    }
  }

  if (P_total == 0.0f)
    weights[point] = 0.0f;
  else
    weights[point] = P_atom / P_total;
}
