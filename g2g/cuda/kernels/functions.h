
// -*- mode: c -*-

/* TODO: coalescear contractions y demas */
template<bool do_forces>
static __device__ __host__ void compute_function(uint m, uint idx, float3 point_position, uint contractions,
  float* factor_a_sh, float* factor_c_sh, uint nuc, float& t, float& tg, float3& v)
{
	float3 atom_nuc_position = gpu_atom_positions[nuc]; // TODO: ver si al usar memoria compartida para esto, pago menos precio por todos los misses
	v = point_position - atom_nuc_position;
	float dist = length2(v);

	t = 0.0f;
	if (do_forces) tg = 0.0f;

	for (uint contraction = 0; contraction < contractions; contraction++) {
		float t0 = expf(-(factor_a_sh[contraction] * dist)) * factor_c_sh[contraction];
		t += t0;
		if (do_forces) tg += t0 * factor_a_sh[contraction];
	}
}

/**
 * gpu_compute_functions
 */
template<bool do_forces>
__global__ void gpu_compute_functions(float4* point_positions, uint points, uint* contractions, float2* factor_ac,
																			uint* nuc, float* function_values, float4* gradient_values, uint4 functions)
{
	dim3 pos = index(blockDim, blockIdx, threadIdx);
	uint point = pos.x;

	/**** Load Point Information ****/
  bool valid_thread = (point < points);
  float3 point_position;
  if (valid_thread) {
    float4 point_position4 = point_positions[point];
    point_position = to_float3(point_position4);
  }

	/** Compute functions ***/
	float t, tg;
	float3 v;

  __shared__ uint nuc_sh[FUNCTIONS_BLOCK_SIZE];
  __shared__ uint contractions_sh[FUNCTIONS_BLOCK_SIZE];
  __shared__ float factor_a_sh[FUNCTIONS_BLOCK_SIZE+1][MAX_CONTRACTIONS];
  __shared__ float factor_c_sh[FUNCTIONS_BLOCK_SIZE+1][MAX_CONTRACTIONS];

  for (uint i = 0; i < functions.w; i += FUNCTIONS_BLOCK_SIZE) {
    if (i + threadIdx.x < functions.w) {
      nuc_sh[threadIdx.x] = nuc[i + threadIdx.x];
      contractions_sh[threadIdx.x] = contractions[i + threadIdx.x];
      for (uint contraction = 0; contraction < contractions_sh[threadIdx.x]; contraction++) {
        float2 factor_ac_local = factor_ac[COALESCED_DIMENSION(functions.w) * contraction + (i + threadIdx.x)];
        factor_a_sh[threadIdx.x][contraction] = factor_ac_local.x;
        factor_c_sh[threadIdx.x][contraction] = factor_ac_local.y;
      }
    }

    __syncthreads();

    // TODO: se podrian evitar los modulos
    if (valid_thread) {
      for (uint ii = 0; ii < FUNCTIONS_BLOCK_SIZE && (i + ii < functions.w); ii++) {
        compute_function<do_forces>(functions.w, ii, point_position, contractions_sh[ii], factor_a_sh[ii], factor_c_sh[ii], nuc_sh[ii], t, tg, v);
        uint idx = COALESCED_DIMENSION(points) * (i + ii) + point;

        if (i + ii < functions.x) {
          function_values[idx] = t;
        }
        else if (i + ii < (functions.x + functions.y * 3)) {
          uint p_idx = ((i + ii) - functions.x) % 3;
          switch(p_idx) {
            case 0: function_values[idx] = v.x * t; break;
            case 1: function_values[idx] = v.y * t; break;
            case 2: function_values[idx] = v.z * t; break;
          }
        }
        else {
          uint d_idx = ((i + ii) - functions.x - functions.y * 3) % 6;
          switch(d_idx) {
            case 0: function_values[idx] = t * v.x * v.x * gpu_normalization_factor; break;
            case 1: function_values[idx] = t * v.y * v.x;                            break;
            case 2: function_values[idx] = t * v.y * v.y * gpu_normalization_factor; break;
            case 3: function_values[idx] = t * v.z * v.x;                            break;
            case 4: function_values[idx] = t * v.z * v.y;                            break;
            case 5: function_values[idx] = t * v.z * v.z * gpu_normalization_factor; break;
          }
        }
      }
    }

    __syncthreads();
  }

  #if 0
	// s functions
	for (uint i = 0; i < functions.x; i++, base_idx++) {
		compute_function<do_forces>(i, point_position, contractions, factor_ac, nuc, spd, t, tg, v);

    uint idx = COALESCED_DIMENSION(points) * base_idx + point;
		function_values[idx] = t;

	}

	// p functions
	for (uint i = 0; i < functions.y; i++, base_idx+=3) {
		compute_function<do_forces>(functions.x + i, point_position, contractions, factor_ac, nuc, spd, t, tg, v);

    uint3 idx;
    idx.x = COALESCED_DIMENSION(points) * (base_idx + 0) + point;
    idx.y = COALESCED_DIMENSION(points) * (base_idx + 1) + point;
    idx.z = COALESCED_DIMENSION(points) * (base_idx + 2) + point;

		function_values[idx.x] = v.x * t;
		function_values[idx.y] = v.y * t;
		function_values[idx.z] = v.z * t;

		if (do_forces) {
      float3 c = v * (2.0f * tg);
			gradient_values[idx.x] = to_float4(v * c.x - make_float3(t, 0, 0));
			gradient_values[idx.y] = to_float4(v * c.y - make_float3(0, t, 0));
			gradient_values[idx.z] = to_float4(v * c.z - make_float3(0, 0, t));
		}
	}

	// d functions
	for (uint i = 0; i < functions.z; i++, base_idx+=6) {
		compute_function<do_forces>(functions.x + functions.y + i, point_position, contractions, factor_ac, nuc, spd, t, tg, v);

    uint3 idx1, idx2;
    idx1.x = COALESCED_DIMENSION(points) * (base_idx + 0) + point;
    idx1.y = COALESCED_DIMENSION(points) * (base_idx + 1) + point;
    idx1.z = COALESCED_DIMENSION(points) * (base_idx + 2) + point;
    idx2.x = COALESCED_DIMENSION(points) * (base_idx + 3) + point;
    idx2.y = COALESCED_DIMENSION(points) * (base_idx + 4) + point;
    idx2.z = COALESCED_DIMENSION(points) * (base_idx + 5) + point;

		function_values[idx1.x] = t * v.x * v.x * gpu_normalization_factor;
		function_values[idx1.y] = t * v.y * v.x;
		function_values[idx1.z] = t * v.y * v.y * gpu_normalization_factor;
		function_values[idx2.x] = t * v.z * v.x;
		function_values[idx2.y] = t * v.z * v.y;
		function_values[idx2.z] = t * v.z * v.z * gpu_normalization_factor;

		if (do_forces) {
			gradient_values[idx1.x] = to_float4(v * 2.0f * tg * v.x * v.x * gpu_normalization_factor - make_float3(2 * t * v.x * gpu_normalization_factor, 0, 0));
			gradient_values[idx1.y] = to_float4(v * 2.0f * tg * v.y * v.x                            - make_float3(t * v.y, t * v.x, 0));
			gradient_values[idx1.z] = to_float4(v * 2.0f * tg * v.y * v.y * gpu_normalization_factor - make_float3(0, 2 * t * v.y * gpu_normalization_factor, 0));
			gradient_values[idx2.x] = to_float4(v * 2.0f * tg * v.z * v.x                            - make_float3(t * v.z, 0, t * v.x));
			gradient_values[idx2.y] = to_float4(v * 2.0f * tg * v.z * v.y                            - make_float3(0, t * v.z, t * v.y));
			gradient_values[idx2.z] = to_float4(v * 2.0f * tg * v.z * v.z * gpu_normalization_factor - make_float3(0, 0, 2 * t * v.z * gpu_normalization_factor));
		}
	}
  #endif
}
