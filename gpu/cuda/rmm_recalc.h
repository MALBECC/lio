/* TODO: juntar esto con la version en functions.h */
__device__ __host__ void compute_function_rmm(float3 point_position, 
  float* factor_a, float* factor_c, uint nuc, float& fj, uint4 functions, uint idx)
{
	float3 atom_nuc_position = gpu_atom_positions[nuc]; // TODO: ver si al usar memoria compartida para esto, pago menos precio por todos los misses
	float3 v = point_position - atom_nuc_position;
	float dist = length2(v);

	float t = 0.0f;

	for (uint contraction = 0; contraction < MAX_CONTRACTIONS; contraction++) {
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

__global__ void gpu_update_rmm(float* factors, uint points, float4* point_weight_positions, uint4 functions, 
  float* rmm, uint* nuc, float2* factor_ac)
{
	uint3 pos = index(blockDim, blockIdx, threadIdx);

	uint i = pos.x; // columna
	uint j = pos.y; // fila
  uint m = functions.w;

  // calculate this rmm
	float rmm_local = 0.0f;

	__shared__ float factor_sh[RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY];
  __shared__ float3 point_positions_sh[RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY];

  __shared__ float fi_sh[RMM_BLOCK_SIZE_XY];
  __shared__ float fj_sh[RMM_BLOCK_SIZE_XY];

  __shared__ uint nuc_i_sh[RMM_BLOCK_SIZE_XY];
  __shared__ uint nuc_j_sh[RMM_BLOCK_SIZE_XY];

  __shared__ float factor_a_i_sh[RMM_BLOCK_SIZE_XY][MAX_CONTRACTIONS];  /* TODO: esto creo que tiene bank conflicts */
  __shared__ float factor_a_j_sh[RMM_BLOCK_SIZE_XY][MAX_CONTRACTIONS];
  __shared__ float factor_c_i_sh[RMM_BLOCK_SIZE_XY][MAX_CONTRACTIONS];
  __shared__ float factor_c_j_sh[RMM_BLOCK_SIZE_XY][MAX_CONTRACTIONS];

  /* en los bloques que no se computa nada (fuera del triangulo inferior izquierdo), directamente hago return */
  /* [ first_fi > first_fj ] */
  if (blockIdx.x * blockDim.x > blockIdx.y * blockDim.y) return;
  
  /* Fi */
  if (threadIdx.y == 0) {
    uint idx = blockIdx.x * blockDim.x;  // first Fi
    if (idx + threadIdx.x < m) {
      nuc_i_sh[threadIdx.x] = nuc[idx + threadIdx.x];
      for (uint contraction = 0; contraction < MAX_CONTRACTIONS; contraction++) {
        float2 factor_ac_local = factor_ac[COALESCED_DIMENSION(m) * contraction + (idx + threadIdx.x)];
        factor_a_i_sh[threadIdx.x][contraction] = factor_ac_local.x;
        factor_c_i_sh[threadIdx.x][contraction] = factor_ac_local.y;
      }
    }
  }

  /* Fj */
  if (blockIdx.x * blockDim.x != blockIdx.y * blockDim.y) {
    if (threadIdx.y == 1) {
      uint idx = blockIdx.y * blockDim.y;  // first Fj
      if (idx + threadIdx.x < m) {
        nuc_j_sh[threadIdx.x] = nuc[idx + threadIdx.x];
        for (uint contraction = 0; contraction < MAX_CONTRACTIONS; contraction++) {
          float2 factor_ac_local = factor_ac[COALESCED_DIMENSION(m) * contraction + (idx + threadIdx.x)];
          factor_a_j_sh[threadIdx.x][contraction] = factor_ac_local.x;
          factor_c_j_sh[threadIdx.x][contraction] = factor_ac_local.y;
        }
      }
    }
  }

	for (uint point_base = 0; point_base < points; point_base += (RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY)) {
		uint abs_threadIdx = threadIdx.y * blockDim.x + threadIdx.x;  // absolute threadId inside block

		if (point_base + abs_threadIdx < points) {
			factor_sh[abs_threadIdx] = factors[point_base + abs_threadIdx];
      float4 point_data = point_weight_positions[point_base + abs_threadIdx];
      point_positions_sh[abs_threadIdx] = to_float3(point_data);
    }

    //#pragma unroll 16
		for (uint point_sub = 0; point_sub < (RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY); point_sub++) {
      if (point_base + point_sub < points) {

        __syncthreads();

        /* compute Fi */
        if (threadIdx.y == 0) {
          uint idx = blockIdx.x * blockDim.x;  // first Fi
          compute_function_rmm(point_positions_sh[point_sub], factor_a_i_sh[threadIdx.x], factor_c_i_sh[threadIdx.x],
            nuc_i_sh[threadIdx.x], fi_sh[threadIdx.x], functions, idx + threadIdx.x);
        }

        if (blockIdx.x * blockDim.x != blockIdx.y * blockDim.y) {
          /* copmute Fj */
          if (threadIdx.y == 1) {
            uint idx = blockIdx.y * blockDim.y;  // first Fj
            compute_function_rmm(point_positions_sh[point_sub], factor_a_j_sh[threadIdx.x], factor_c_j_sh[threadIdx.x],
              nuc_j_sh[threadIdx.x], fj_sh[threadIdx.x], functions, idx + threadIdx.x);
          }
        }

        __syncthreads();

        /* NOTE: on blocks with some valid threads, the computation contributes 0 to rmm_local */
        if (blockIdx.x * blockDim.x == blockIdx.y * blockDim.y)
          rmm_local += fi_sh[threadIdx.x] * fi_sh[threadIdx.y] * factor_sh[point_sub];
        else
          rmm_local += fi_sh[threadIdx.x] * fj_sh[threadIdx.y] * factor_sh[point_sub];
      }
    }
	}

	/* quiero triangulo inferior solamente TODO: sacar esto */
  if (i<= j && i < m && j < m) rmm[COALESCED_DIMENSION(m) * j + i] = rmm_local;
}
