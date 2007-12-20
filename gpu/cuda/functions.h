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
__device__ void calc_single_function_p(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																			 const float3* atom_positions, const float* factor_a, const float* factor_c, uint func_index, uint subfunc, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float dist = distance2(point_position, atom_positions[atom_nuc]);

	uint func_contractions = contractions[func_index];
	*func_value = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float rexp = factor_a[func_index * MAX_CONTRACTIONS + contraction] * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * factor_c[func_index * MAX_CONTRACTIONS + contraction];
		
		float v = (elem(point_position,subfunc) - elem(atom_positions[atom_nuc],subfunc)) * t;
		
		*func_value += v;
	}
}

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
__device__ void calc_single_function_d(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
																			 const float3* atom_positions, const float* factor_a, const float* factor_c, uint func_index, uint subfunc,
																			 float normalization_factor, float* func_value)
{
	uint atom_nuc = nuc[func_index];
	float dist = distance2(point_position, atom_positions[atom_nuc]);

	uint func_contractions = contractions[func_index];
	*func_value = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float rexp = factor_a[func_index * MAX_CONTRACTIONS + contraction] * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * factor_c[func_index * MAX_CONTRACTIONS + contraction];

		float3 v = point_position - atom_positions[atom_nuc];
		
		float t1 = 1.0f, t2 = 1.0f;
		switch(subfunc) {
			case 0:
				t1 = elem(v, 0);
				t2 = elem(v, 0) * normalization_factor;
			break;
			case 1:
				t1 = elem(v, 1);
				t2 = elem(v, 0);			
			break;
			case 2:
				t1 = elem(v, 1);
				t2 = elem(v, 1) * normalization_factor;			
			break;
			case 3:
				t1 = elem(v, 2);
				t2 = elem(v, 0);			
			break;
			case 4:
				t1 = elem(v, 2);
				t2 = elem(v, 1);			
			break;
			case 5:
				t1 = elem(v, 2);
				t2 = elem(v, 2) * normalization_factor;			
			break;
		}
		
		*func_value += t * t1 * t2;
	}
}

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

#if 0
__device__ void calc_function(const uint3& num_funcs, const uint* nuc, const uint* contractions, const float3& point_position,
															const float3* atom_positions, const float* factor_a, const float* factor_c, uint big_func_index,
															bool normalize, float& func_value)
{
	uint funcs_s = num_funcs.x;
	uint funcs_p = num_funcs.y;
	
	uint funcs_sp = funcs_s + funcs_p * 3;

	if (big_func_index < funcs_s) {
		calc_function_s(num_funcs, nuc, contractions, point_position, atom_positions, factor_a, factor_c, big_func_index, &func_value);
	}
	else if (big_func_index < funcs_sp) {
		uint func_index = (big_func_index - funcs_s) / 3;
		float complete_func_value[3];
		calc_function_p(num_funcs, nuc, contractions, point_position, atom_positions, factor_a, factor_c, func_index, complete_func_value);
		func_value = complete_func_value[(big_func_index - funcs_s) % 3];
	}
	else {
		float normalization_factor = (normalize ? rsqrtf(3.0f) : 1.0f);
		uint func_index = (big_func_index - funcs_s - funcs_p * 3) / 6;
		float complete_func_value[6];
		calc_function_d(num_funcs, nuc, contractions, point_position, atom_positions, factor_a, factor_c, func_index,
										normalization_factor, complete_func_value);
		func_value = complete_func_value[(big_func_index - funcs_s - funcs_p * 3) % 6];
	}
}
#endif
