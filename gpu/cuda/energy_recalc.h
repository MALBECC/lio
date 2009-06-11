
/* TODO: juntar esto con la version en functions.h */
__device__ void compute_function_(float3 point_position, uint contractions,
  float* factor_a, float* factor_c, uint nuc, float& fj, uint4 functions, uint idx)
{
	float3 atom_nuc_position = gpu_atom_positions[nuc]; // TODO: ver si al usar memoria compartida para esto, pago menos precio por todos los misses
	float3 v = point_position - atom_nuc_position;
	float dist = length2(v);

	float t = 0.0f;

	for (uint contraction = 0; contraction < contractions; contraction++) {
		t += expf(-(factor_a[contraction] * dist)) * factor_c[contraction];
	}

  if (idx < functions.x) { fj = t; }
  else if (idx < (functions.x + functions.y * 3)) {
    uint p_idx = (idx - functions.x) % 3;   // TODO: arreglar este modulo y el de abajo, son lentos
    switch(p_idx) {
      case 0: fj = v.x * t; break;
      case 1: fj = v.y * t; break;
      case 2: fj = v.z * t; break;
    }
  }
  else {
    uint d_idx = (idx - functions.x - functions.y * 3) % 6;
    switch(d_idx) {
      case 0: fj = t * v.x * v.x * gpu_normalization_factor; break;
      case 1: fj = t * v.y * v.x;                            break;
      case 2: fj = t * v.y * v.y * gpu_normalization_factor; break;
      case 3: fj = t * v.z * v.x;                            break;
      case 4: fj = t * v.z * v.y;                            break;
      case 5: fj = t * v.z * v.z * gpu_normalization_factor; break;
    }
  }
}

template<bool compute_energy>
__global__ void gpu_compute_density(float* output, float4* point_weight_positions, uint points, float* rdm,
  uint4 functions, uint2* nuc_contractions, float2* factor_ac)
{
  const uint point = index_x(blockDim, blockIdx, threadIdx);
  const uint m = functions.w;

  float partial_density = 0.0f;

  __shared__ float rdm_sh[NCO_BATCH_SIZE];
  __shared__ uint nuc_sh[DENSITY_BLOCK_SIZE];
  __shared__ uint contractions_sh[DENSITY_BLOCK_SIZE];
  __shared__ float factor_a_sh[DENSITY_BLOCK_SIZE][MAX_CONTRACTIONS];
  __shared__ float factor_c_sh[DENSITY_BLOCK_SIZE][MAX_CONTRACTIONS];
  __shared__ float w_sh[NCO_BATCH_SIZE][DENSITY_BLOCK_SIZE + 1];   // TODO: poner en 0

  bool valid_thread = (point < points);

  float4 point_weight_position;
  if (valid_thread) point_weight_position = point_weight_positions[point];
  float point_weight = point_weight_position.w;
  float3 point_position = to_float3(point_weight_position);

  for (uint i = 0; i < gpu_nco; i += NCO_BATCH_SIZE) {
    #pragma unroll 16
    for (uint ii = 0; ii < NCO_BATCH_SIZE; ii++) w_sh[ii][threadIdx.x] = 0.0f;

    for (uint j = 0; j < m; j += DENSITY_BLOCK_SIZE) {
      if (j + threadIdx.x < m) {
        uint2 nuc_contraction = nuc_contractions[j + threadIdx.x];
        nuc_sh[threadIdx.x] = nuc_contraction.x;
        contractions_sh[threadIdx.x] = nuc_contraction.y;
        for (uint contraction = 0; contraction < contractions_sh[threadIdx.x]; contraction++) {
          float2 factor_ac_local = factor_ac[COALESCED_DIMENSION(m) * contraction + (j + threadIdx.x)];
          factor_a_sh[threadIdx.x][contraction] = factor_ac_local.x;
          factor_c_sh[threadIdx.x][contraction] = factor_ac_local.y;
        }
      }

      __syncthreads();

      for (uint jj = 0; jj < DENSITY_BLOCK_SIZE; jj++) {
        if (j + jj < m) {
          float fj;
          compute_function_(point_position, contractions_sh[jj], factor_a_sh[jj], factor_c_sh[jj], nuc_sh[jj], fj, functions, j + jj);

          __syncthreads();

          if (i + threadIdx.x < NCO_BATCH_SIZE) {
            if  (i + threadIdx.x < gpu_nco) rdm_sh[threadIdx.x] = rdm[COALESCED_DIMENSION(gpu_nco) * (j + jj) + (i + threadIdx.x)];
            else rdm_sh[threadIdx.x] = 0.0f;
          }

          __syncthreads();

          #pragma unroll 16
          for (uint ii = 0; ii < NCO_BATCH_SIZE; ii++) w_sh[ii][threadIdx.x] += fj * rdm_sh[ii];
        }
      }
    }

    #pragma unroll 16
    for (uint ii = 0; ii < NCO_BATCH_SIZE; ii++) partial_density += w_sh[ii][threadIdx.x] * w_sh[ii][threadIdx.x];
  }

  partial_density *= 2.0f;

  /***** compute energy / factor *****/
  float exc_corr, y2a;
  if (compute_energy) {
    gpu_pot<true, false>(partial_density, exc_corr, y2a);
    if (valid_thread) output[point] = (partial_density * point_weight) * exc_corr;
	}
	else {
		gpu_pot<false, true>(partial_density, exc_corr, y2a);
    if (valid_thread) output[point] = point_weight * y2a;
  }
}
