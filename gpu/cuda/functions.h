// -*- mode: c -*-

// ********* S ***********
__device__ void calc_function_s(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																const float3* atom_positions, const float3* atom_positions_shared,
																const float2* factor_ac, uint func_index, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float3 atom_nuc_position = get_atom_position(atom_nuc, atom_positions_shared, atom_positions);
	float dist = distance2(point_position, atom_nuc_position);

	uint func_contractions = contractions[func_index];
	*func_value = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[func_index * MAX_CONTRACTIONS + contraction];
		float rexp = curr_factor_ac.x * dist;
		if (rexp > 30.0f) continue;
		*func_value += expf(-rexp) * curr_factor_ac.y;
	}	
}


// ****** P *******
__device__ void calc_function_p(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																const float3* atom_positions, const float3* atom_positions_shared,
																const float2* factor_ac, uint func_index, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float3 atom_nuc_position = get_atom_position(atom_nuc, atom_positions_shared, atom_positions);
	float dist = distance2(point_position, atom_nuc_position);

	uint func_contractions = contractions[func_index];
	for (uint i = 0; i < 3; i++) func_value[i] = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[func_index * MAX_CONTRACTIONS + contraction];
		float rexp = curr_factor_ac.x * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * curr_factor_ac.y;
		float3 v = (point_position - atom_nuc_position) * t;
		
		for (uint i = 0; i < 3; i++) func_value[i] += elem(v,i);
	}
}

// ************** D **************
__device__ void calc_function_d(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																const float3* atom_positions, const float3* atom_positions_shared,
																const float2* factor_ac, uint func_index, float normalization_factor, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float3 atom_nuc_position = get_atom_position(atom_nuc, atom_positions_shared, atom_positions);	
	float dist = distance2(point_position, atom_nuc_position);

	uint func_contractions = contractions[func_index];
	for (uint i = 0; i < 6; i++) func_value[i] = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[func_index * MAX_CONTRACTIONS + contraction];
		float rexp = curr_factor_ac.x * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * curr_factor_ac.y;

		float3 v = point_position - atom_nuc_position;

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

