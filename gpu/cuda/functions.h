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
		float rexp = curr_factor_ac.x * dist;
		//if (rexp > 30.0f) continue;
		float t0 = expf(-rexp) * curr_factor_ac.y;
		t += t0;
		if (do_forces) tg += t0 * curr_factor_ac.x;
		//_EMU(printf("func(%i): %f, a: %f, c: %f (%i/%i contractions) spd: %i\n", idx, t, curr_factor_ac.x, curr_factor_ac.y, contraction, func_contractions, spd));
	}
}

/**
 * gpu_compute_functions
 */
template<bool do_forces>
__global__ void gpu_compute_functions(float3* point_positions, uint points, uint* contractions, float2* factor_ac,
																			uint* nuc, float* function_values, float3* gradient_values, uint4 functions, uint spd)
{
	dim3 pos = index(blockDim, blockIdx, threadIdx);
	uint point = pos.x;
	
	/**** Load Point Information ****/
	if (point >= points) return; 	// esto se puede evitar computando basura	
	float3 point_position = point_positions[point];

	/** Compute functions ***/
	//uint base = point * COALESCED_DIMENSION(functions.w);
  uint func_base = 0;
  uint grad_base = point * functions.w;

	float t, tg;
	float3 v;
	
	// s functions
	for (uint i = 0; i < functions.x; i++, func_base++, grad_base++) {
		compute_function<do_forces>(i, point_position, contractions, factor_ac, nuc, spd, t, tg, v);

    uint base = COALESCED_DIMENSION(points) * func_base + point;
		function_values[base] = t;
		if (do_forces) gradient_values[grad_base] = v * (2.0f * tg);
	}
	
	// p functions
	for (uint i = 0; i < functions.y; i++, func_base+=3, grad_base+=3) {
		compute_function<do_forces>(functions.x + i, point_position, contractions, factor_ac, nuc, spd, t, tg, v);

    
		function_values[COALESCED_DIMENSION(points) * (func_base + 0) + point] = v.x * t;
		function_values[COALESCED_DIMENSION(points) * (func_base + 1) + point] = v.y * t;
		function_values[COALESCED_DIMENSION(points) * (func_base + 2) + point] = v.z * t;

		if (do_forces) {
			gradient_values[grad_base + 0] = v * 2.0f * v.x * tg - make_float3(t, 0, 0);
			gradient_values[grad_base + 1] = v * 2.0f * v.y * tg - make_float3(0, t, 0);
			gradient_values[grad_base + 2] = v * 2.0f * v.z * tg - make_float3(0, 0, t);
		}
	}
	
	// d functions
	for (uint i = 0; i < functions.z; i++, func_base+=6, grad_base+=6) {
		compute_function<do_forces>(functions.x + functions.y + i, point_position, contractions, factor_ac, nuc, spd, t, tg, v);
		
		float tx = t * v.x;
		float ty = t * v.y;
		float tz = t * v.z;

		function_values[COALESCED_DIMENSION(points) * (func_base + 0) + point] = tx * v.x * gpu_normalization_factor;
		function_values[COALESCED_DIMENSION(points) * (func_base + 1) + point] = ty * v.x;
		function_values[COALESCED_DIMENSION(points) * (func_base + 2) + point] = ty * v.y * gpu_normalization_factor;
		function_values[COALESCED_DIMENSION(points) * (func_base + 3) + point] = tz * v.x;
		function_values[COALESCED_DIMENSION(points) * (func_base + 4) + point] = tz * v.y;
		function_values[COALESCED_DIMENSION(points) * (func_base + 5) + point] = tz * v.z * gpu_normalization_factor;

		if (do_forces) {
			float tgx = tg * v.x;
			float tgy = tg * v.y;
			float tgz = tg * v.z;

			gradient_values[grad_base + 0] = v * 2.0f * tgx * v.x * gpu_normalization_factor - make_float3(2 * tx * gpu_normalization_factor, 0, 0);
			gradient_values[grad_base + 1] = v * 2.0f * tgy * v.x                            - make_float3(ty, tx, 0);
			gradient_values[grad_base + 2] = v * 2.0f * tgy * v.y * gpu_normalization_factor - make_float3(0, 2 * ty * gpu_normalization_factor, 0);
			gradient_values[grad_base + 3] = v * 2.0f * tgz * v.x                            - make_float3(tz, 0, tx);
			gradient_values[grad_base + 4] = v * 2.0f * tgz * v.y                            - make_float3(0, tz, ty);
			gradient_values[grad_base + 5] = v * 2.0f * tgz * v.z * gpu_normalization_factor - make_float3(0, 0, 2 * tz * gpu_normalization_factor);
		}
	}
}
