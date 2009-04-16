// -*- mode: c -*-

/* TODO: coalescear contractions y demas */
template<bool do_forces>
__device__ __host__ void compute_function(uint idx, float3 point_position, uint* contractions, float2* factor_ac, uint* nuc,
																uint spd, float& t, float& tg, float3& v)
{
  _EMU(printf("idx: %i nuc: %i cont: %i\n", idx, nuc[idx], contractions[idx]));
	float3 atom_nuc_position = gpu_atom_positions[nuc[idx]];
	v = point_position - atom_nuc_position;
	float dist = length2(v);
	uint func_contractions = contractions[idx];

	t = 0.0f;
	if (do_forces) tg = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[contraction * spd + idx];
		float t0 = expf(-(curr_factor_ac.x * dist)) * curr_factor_ac.y;
		t += t0;
		if (do_forces) tg += t0 * curr_factor_ac.x;
	}
}

/**
 * gpu_compute_functions
 */
template<bool do_forces>
__global__ void gpu_compute_functions(float3* point_positions, uint points, uint* contractions, float2* factor_ac,
																			uint* nuc, float* function_values, float4* gradient_values, uint4 functions, uint spd)
{
	dim3 pos = index(blockDim, blockIdx, threadIdx);
	uint point = pos.x;
	
	/**** Load Point Information ****/
	if (point >= points) return; 	// esto se puede evitar computando basura	
	float3 point_position = point_positions[point];

	/** Compute functions ***/
  uint base_idx = 0;

	float t, tg;
	float3 v;

	// s functions
	for (uint i = 0; i < functions.x; i++, base_idx++) {
		compute_function<do_forces>(i, point_position, contractions, factor_ac, nuc, spd, t, tg, v);

    uint idx = COALESCED_DIMENSION(points) * base_idx + point;
		function_values[idx] = t;
		if (do_forces) { gradient_values[idx] = to_float4(v * (2.0f * tg)); }
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
}
