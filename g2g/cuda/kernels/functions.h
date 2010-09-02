
// -*- mode: c -*-

/* TODO: coalescear contractions y demas */
template<bool do_forces, bool do_gga>
static __device__ __host__ void compute_function(uint m, uint idx, float3 point_position, uint contractions,
  float* factor_a_sh, float* factor_c_sh, uint nuc, float& t, float& tg, float& th, float3& v)
{
	float3 atom_nuc_position = gpu_atom_positions[nuc]; // TODO: ver si al usar memoria compartida para esto, pago menos precio por todos los misses
	v = point_position - atom_nuc_position;
	float dist = length2(v);

	t = 0.0f;
	if (do_forces || do_gga) tg = 0.0f;
  if (do_gga) th = 0.0f;

	for (uint contraction = 0; contraction < contractions; contraction++) {
		float t0 = expf(-(factor_a_sh[contraction] * dist)) * factor_c_sh[contraction];
		t += t0;
		if (do_forces || do_gga) tg += t0 * factor_a_sh[contraction];
    if (do_gga) th += t0 * (factor_a_sh[contraction] * factor_a_sh[contraction]);
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
	float t, tg, th;
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
        compute_function<do_forces, do_gga>(functions.w, ii, point_position, contractions_sh[ii], factor_a_sh[ii], factor_c_sh[ii], nuc_sh[ii], t, tg, th, v);
        uint idx = COALESCED_DIMENSION(points) * (i + ii) + point;
        uint hidx1, hidx2;
        if (do_gga) {
          hidx1 = COALESCED_DIMENSION(points) * (2 * (i + ii) + 0) + point;
          hidx2 = COALESCED_DIMENSION(points) * (2 * (i + ii) + 1) + point;
        }

				float3 vxxy, vyzz;
				if (do_gga) { vxxy = make_float3(v.x, v.x, v.y); vyzz = make_float3(v.y, v.z, v.z); }

        if (i + ii < functions.x) {
          function_values[idx] = t;
          if (do_forces || do_gga) gradient_values[idx] = to_float4(v * (-2.0f * tg));
          if (do_gga) {
            hessian_values[hidx1] = to_float4((v * v) * 4.0f * th - 2.0f * tg); // Fxx, Fxy, Fxz
            hessian_values[hidx2] = to_float4(vxxy * vyzz * 4.0f * th); // Fxy, Fxz, Fyz
          }
        }
        else if (i + ii < (functions.x + functions.y * 3)) {
          uint p_idx = ((i + ii) - functions.x) % 3;
          switch(p_idx) {
            case 0:
              function_values[idx] = v.x * t;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(t, 0.0f, 0.0f) - v * 2.0f * tg * v.x);
              if (do_gga) {
                hessian_values[hidx1] = to_float4(4.0f * th * v.x * (v * v) - tg * v.x * make_float3(6.0f, 2.0f, 2.0f));
                hessian_values[hidx2] = to_float4(4.0f * th * v.x * (vxxy * vyzz) - 2.0f * tg * make_float3(v.y, v.z, 0.0f));
              }
            break;
            case 1:
              function_values[idx] = v.y * t;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(0.0f, t, 0.0f) - v * 2.0f * tg * v.y);
              if (do_gga) {
                hessian_values[hidx1] = to_float4(4.0f * th * v.y * (v * v) - tg * v.y * make_float3(2.0f, 6.0f, 2.0f));
                hessian_values[hidx2] = to_float4(4.0f * th * v.y * (vxxy * vyzz) - 2.0f * tg * make_float3(v.x, 0.0f, v.z));
              }
            break;
            case 2:
              function_values[idx] = v.z * t;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(0.0f, 0.0f, t) - v * 2.0f * tg * v.z);
              if (do_gga) {
                hessian_values[hidx1] = to_float4(4.0f * th * v.z * (v * v) - tg * v.z * make_float3(2.0f, 2.0f, 6.0f));
                hessian_values[hidx2] = to_float4(4.0f * th * v.z * (vxxy * vyzz) - 2.0f * tg * make_float3(0.0f, v.x, v.y));
              }
            break;
          }
        }
        else {
          uint d_idx = ((i + ii) - functions.x - functions.y * 3) % 6;
          switch(d_idx) {
            case 0:
              function_values[idx] = t * v.x * v.x * gpu_normalization_factor;
              if (do_forces || do_gga) gradient_values[idx] = to_float4((make_float3(2.0f * v.x, 0.0f      , 0.0f      ) * t - 2.0f * tg * v * v.x * v.x) * gpu_normalization_factor);
              if (do_gga) {
                hessian_values[hidx1] = to_float4((4.0f * th * (v.x * v.x) * (v * v)       - tg * (v.x * v.x)   * make_float3(10.0f, 2.0f, 2.0f) + make_float3(2.0f * t, 0.0f    , 0.0f)) * gpu_normalization_factor);
                hessian_values[hidx2] = to_float4((4.0f * th * (v.x * v.x) * (vxxy * vyzz) - tg * (vxxy * vyzz) * make_float3(4.0f,  4.0f, 0.0f)                               ) * gpu_normalization_factor);
              }
            break;
            case 1:
              function_values[idx] = t * v.y * v.x;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(v.y        , v.x       , 0.0f      ) * t - 2.0f * tg * v * v.y * v.x);
              if (do_gga) {
                hessian_values[hidx1] = to_float4((4.0f * th * (v.x * v.y) * (v * v)       - tg * (v.x * v.y)   * make_float3(6.0f,  6.0f, 2.0f)                               ));
                hessian_values[hidx2] = to_float4((4.0f * th * (v.x * v.y) * (vxxy * vyzz) - tg * make_float3(2.0f * (v.x * v.x + v.y * v.y), 2.0f * v.y * v.z, 2.0f * v.x * v.z) + make_float3(t     , 0.0f    , 0.0f)));
              }
            break;
            case 2:
              function_values[idx] = t * v.y * v.y * gpu_normalization_factor;
              if (do_forces || do_gga) gradient_values[idx] = to_float4((make_float3(0.0f      , 2.0f * v.y, 0.0f      ) * t - 2.0f * tg * v * v.y * v.y) * gpu_normalization_factor);
              if (do_gga) {
                hessian_values[hidx1] = to_float4((4.0f * th * (v.y * v.y) * (v * v)       - tg * (v.y * v.y)   * make_float3(2.0f, 10.0f, 2.0f) + make_float3(0.0f    , 2.0f * t, 0.0f)) * gpu_normalization_factor);
                hessian_values[hidx2] = to_float4((4.0f * th * (v.y * v.y) * (vxxy * vyzz) - tg * (vxxy * vyzz) * make_float3(4.0f,  0.0f, 4.0f)                                ) * gpu_normalization_factor);
              }
            break;
            case 3:
              function_values[idx] = t * v.z * v.x;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(v.z        , 0.0f      , v.x       ) * t - 2.0f * tg * v * v.z * v.x);
              if (do_gga) {
                hessian_values[hidx1] = to_float4((4.0f * th * (v.x * v.z) * (v * v)       - tg * (v.x * v.z)   * make_float3(6.0f,  2.0f, 6.0f)                                ));
                hessian_values[hidx2] = to_float4((4.0f * th * (v.x * v.z) * (vxxy * vyzz) - tg * make_float3(2.0f * v.y * v.z, 2.0f * (v.x * v.x + v.z * v.z), 2.0f * v.x * v.y) + make_float3(0.0f,      t,     0.0f)));
              }
            break;
            case 4:
              function_values[idx] = t * v.z * v.y;
              if (do_forces || do_gga) gradient_values[idx] = to_float4(make_float3(0.0f       , v.z       , v.y       ) * t - 2.0f * tg * v * v.z * v.y);
              if (do_gga) {
                hessian_values[hidx1] = to_float4((4.0f * th * (v.y * v.z) * (v * v)       - tg * (v.y * v.z)   * make_float3(2.0f,  6.0f, 6.0f)                                ));
                hessian_values[hidx2] = to_float4((4.0f * th * (v.y * v.z) * (vxxy * vyzz) - tg * make_float3(2.0f * v.x * v.z, 2.0f * v.x * v.y, 2.0f * (v.y * v.y + v.z * v.z)) + make_float3(0.0f,      0.0f,     t)));
              }
            break;
            case 5:
              function_values[idx] = t * v.z * v.z * gpu_normalization_factor;
              if (do_forces || do_gga) gradient_values[idx] = to_float4((make_float3(0.0f      , 0.0f      , 2.0f * v.z) * t - 2.0f * tg * v * v.z * v.z) * gpu_normalization_factor);
              if (do_gga) {
                hessian_values[hidx1] = to_float4((4.0f * th * (v.z * v.z) * (v * v)       - tg * (v.z * v.z)   * make_float3(2.0f,  2.0f, 10.0f) + make_float3(0.0f,      0.0f, 2.0f * t)) * gpu_normalization_factor);
                hessian_values[hidx2] = to_float4((4.0f * th * (v.z * v.z) * (vxxy * vyzz) - tg * (vxxy * vyzz) * make_float3(0.0f,  4.0f, 4.0f)                                 ) * gpu_normalization_factor);
              }
            break;
          }
        }
      }
    }

    __syncthreads();
  }
}
