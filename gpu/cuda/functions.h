// -*- mode: c -*-
__device__ float compute_function(uint idx, float3 point_position, uint* contractions, float2* factor_ac, float3 atom_nuc_position,
                                  uint spd)
{
	float dist = distance2(point_position, atom_nuc_position);
	uint func_contractions = contractions[idx];

	float t = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[contraction * spd + idx];
		float rexp = curr_factor_ac.x * dist;
		if (rexp > 30.0f) continue; // TODO: esto se puede sacar?
		t += expf(-rexp) * curr_factor_ac.y;
		//_EMU(printf("func(%i): %f, a: %f, c: %f (%i/%i contractions) spd: %i\n", idx, t, curr_factor_ac.x, curr_factor_ac.y, contraction, func_contractions, spd));
	}
	return t;
}

/**
 * gpu_compute_functions
 * 	grid: points 
 *  TODO: agregar otra dimension para parelelizar las funciones
 *  TODO: paralelizar por puntos (pueden cargar cosas comunes)
 */
__global__ void gpu_compute_functions(float3* point_positions, uint points, uint* contractions, float2* factor_ac, uint* nuc,
																			float* function_values, uint4 functions, uint spd)
{
	dim3 pos = index(blockDim, blockIdx, threadIdx);
	uint point = pos.x;
	
	/**** Load Point Information ****/
	if (point >= points) return; 	// esto se puede evitar computando basura	
	float3 point_position = point_positions[point];

	/** Compute functions ***/
	uint base = point * functions.w;
	
	// s functions
	for (uint i = 0; i < functions.x; i++) {
		float3 atom_nuc_position = gpu_atom_positions[nuc[i]];
		function_values[base + i] = compute_function(i, point_position, contractions, factor_ac, atom_nuc_position,spd);
	}
	
	// p functions
	base += functions.x;
	for (uint i = 0; i < functions.y; i++) {
		float3 atom_nuc_position = gpu_atom_positions[nuc[functions.x + i]];
		float t = compute_function(functions.x + i, point_position, contractions, factor_ac, atom_nuc_position,spd);
		for (uint j = 0; j < 3; j++) {
			function_values[base + i * 3 + j] = (float3_elem(point_position, j) - float3_elem(atom_nuc_position, j)) * t;
		}
	}
	
	// d functions
	base += functions.y * 3;
	for (uint i = 0; i < functions.z; i++) {
		float3 atom_nuc_position = gpu_atom_positions[nuc[functions.x + functions.y + i]];
		float t = compute_function(functions.y + functions.x + i, point_position, contractions, factor_ac, atom_nuc_position,spd);
	
		float3 v = point_position - atom_nuc_position;

		uint jk = 0;
		for (uint j = 0; j < 3; j++) {
			float t1 = float3_elem(v, j);
			for (uint k = 0; k <= j; k++) {
				float t2 = (j == k ? float3_elem(v, k) * gpu_normalization_factor : float3_elem(v, k));
				function_values[base + i * 6 + jk] = t * t1 * t2;
				jk++;
			}
		}
	}
}
