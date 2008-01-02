// -*- mode: c -*-

__device__ void calc_function_s(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																const float3* atom_positions, const float* factor_a, const float* factor_c, uint func_index, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float dist = distance2(point_position, atom_positions[atom_nuc]);

	uint func_contractions = contractions[func_index];
	*func_value = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float rexp = factor_a[func_index * MAX_CONTRACTIONS + contraction] * dist;
		if (rexp > 30.0f) continue;
		*func_value += expf(-rexp) * factor_c[func_index * MAX_CONTRACTIONS + contraction];
	}	
}


// ****** P *******
__device__ void calc_function_p(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																const float3* atom_positions, const float* factor_a, const float* factor_c, uint func_index, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float dist = distance2(point_position, atom_positions[atom_nuc]);

	uint func_contractions = contractions[func_index];
	for (uint i = 0; i < 3; i++) func_value[i] = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float rexp = factor_a[func_index * MAX_CONTRACTIONS + contraction] * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * factor_c[func_index * MAX_CONTRACTIONS + contraction];
		float3 v = (point_position - atom_positions[atom_nuc]) * t;
		
		for (uint i = 0; i < 3; i++) func_value[i] += elem(v,i);
	}
}

// ************** D **************
__device__ void calc_function_d(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																const float3* atom_positions, const float* factor_a, const float* factor_c, uint func_index,
																float normalization_factor, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float dist = distance2(point_position, atom_positions[atom_nuc]);

	uint func_contractions = contractions[func_index];
	for (uint i = 0; i < 6; i++) func_value[i] = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float rexp = factor_a[func_index * MAX_CONTRACTIONS + contraction] * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * factor_c[func_index * MAX_CONTRACTIONS + contraction];

		float3 v = point_position - atom_positions[atom_nuc];

		uint ij = 0;
		for (uint i = 0; i < 3; i++) {
			float t1 = elem(v, i);
			for (uint j = 0; j <= i; j++) {
				float t2 = (i == j ? elem(v, j) * normalization_factor : elem(v, j));
				func_value[ij] += t * t1 * t2;
				ij++;
			}
		}
	}
}

