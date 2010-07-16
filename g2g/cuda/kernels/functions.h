
// -*- mode: c -*-

/* TODO: coalescear contractions y demas */
template<bool do_forces, bool do_gga>
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
		if (do_forces || do_gga) tg += t0 * factor_a_sh[contraction];
	}
}

/**
 * gpu_compute_functions
 */
template<bool do_forces, bool do_gga>
__global__ void gpu_compute_functions(float4* point_positions, uint points, uint* contractions, float2* factor_ac,
																			uint* nuc, float* function_values, float4* gradient_values, float4* hessian_values, uint4 functions)
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
        compute_function<do_forces, do_gga>(functions.w, ii, point_position, contractions_sh[ii], factor_a_sh[ii], factor_c_sh[ii], nuc_sh[ii], t, tg, v);
        uint idx = COALESCED_DIMENSION(points) * (i + ii) + point;

        if (i + ii < functions.x) {
          function_values[idx] = t;
          if (do_forces) gradient_values[idx] = to_float4(v * (-2.0f * tg));
        }
        else if (i + ii < (functions.x + functions.y * 3)) {
          uint p_idx = ((i + ii) - functions.x) % 3;
          switch(p_idx) {
            case 0:
              function_values[idx] = v.x * t;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(t, 0.0f, 0.0f) - v * 2.0f * tg * v.x);
            break;
            case 1:
              function_values[idx] = v.y * t;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(0.0f, t, 0.0f) - v * 2.0f * tg * v.y);
            break;
            case 2:
              function_values[idx] = v.z * t;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(0.0f, 0.0f, t) - v * 2.0f * tg * v.z);
            break;
          }
        }
        else {
          uint d_idx = ((i + ii) - functions.x - functions.y * 3) % 6;
          switch(d_idx) {
            case 0:
              function_values[idx] = t * v.x * v.x * gpu_normalization_factor;
              if (do_forces || do_gga) gradient_values[idx] = to_float4((make_float3(2.0f * v.x, 0.0f      , 0.0f      ) * t - 2.0f * tg * v * v.x * v.x) * gpu_normalization_factor);
            break;
            case 1:
              function_values[idx] = t * v.y * v.x;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(v.y        , v.x       , 0.0f      ) * t - 2.0f * tg * v * v.y * v.x);
            break;
            case 2:
              function_values[idx] = t * v.y * v.y * gpu_normalization_factor;
              if (do_forces || do_gga) gradient_values[idx] = to_float4((make_float3(0.0f      , 2.0f * v.y, 0.0f      ) * t - 2.0f * tg * v * v.y * v.y) * gpu_normalization_factor);
            break;
            case 3:
              function_values[idx] = t * v.z * v.x;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(v.z        , 0.0f      , v.x       ) * t - 2.0f * tg * v * v.z * v.x);
            break;
            case 4:
              function_values[idx] = t * v.z * v.y;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(0.0f       , v.z       , v.y       ) * t - 2.0f * tg * v * v.z * v.y);
            break;
            case 5:
              function_values[idx] = t * v.z * v.z * gpu_normalization_factor;
              if (do_forces || do_gga) gradient_values[idx] = to_float4((make_float3(0.0f      , 0.0f      , 2.0f * v.z) * t - 2.0f * tg * v * v.z * v.z) * gpu_normalization_factor);
            break;
          }
        }
      }
    }

    __syncthreads();
  }
}
