// -*- mode: c -*-

// ********* S ***********
__device__ void calc_function_s(const uint3& num_funcs, const uint2* nuc_contractions, const float3& point_position,
																const float2* factor_ac, uint func_index, float* func_value)
{
	uint atom_nuc = nuc_contractions[func_index].x;
	float3 atom_nuc_position = gpu_atom_positions[atom_nuc];
	float dist = distance2(point_position, atom_nuc_position);

	uint func_contractions = nuc_contractions[func_index].y;
	float func_value_local = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[func_index * MAX_CONTRACTIONS + contraction];
		float rexp = curr_factor_ac.x * dist;
		if (rexp > 30.0f) continue;
		func_value_local += expf(-rexp) * curr_factor_ac.y;
	}
	
	*func_value = func_value_local;
}

__device__ void calc_function_s(const uint3& num_funcs, const uint2* nuc_contractions, const float3& point_position,
															  const float2* factor_ac, uint func_index, float* func_value, float3* Fg)
{
	uint atom_nuc = nuc_contractions[func_index].x;
	float3 atom_nuc_position = gpu_atom_positions[atom_nuc];
	float dist = distance2(point_position, atom_nuc_position);

	uint func_contractions = nuc_contractions[func_index].y;
	float func_value_local = 0.0f;
	float Fgx = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[func_index * MAX_CONTRACTIONS + contraction];
		float rexp = curr_factor_ac.x * dist;
		if (rexp > 30.0f) continue;
		
		float term = expf(-rexp) * curr_factor_ac.y;
		func_value_local += term;
    Fgx += term * curr_factor_ac.x;
	}

	*Fg = (point_position - atom_nuc_position) * 2.0f * Fgx;
	*func_value = func_value_local;
}

// ****** P *******
__device__ void calc_function_p(const uint3& num_funcs, const uint2* nuc_contractions, const float3& point_position,
																const float2* factor_ac, uint func_index, float* func_value)
{
	uint atom_nuc = nuc_contractions[func_index].x;
	float3 atom_nuc_position = gpu_atom_positions[atom_nuc];
	float dist = distance2(point_position, atom_nuc_position);

	uint func_contractions = nuc_contractions[func_index].y;
	
	float func_value_local[3];
	for (uint i = 0; i < 3; i++) func_value_local[i] = 0.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[func_index * MAX_CONTRACTIONS + contraction];
		float rexp = curr_factor_ac.x * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * curr_factor_ac.y;
		float3 v = (point_position - atom_nuc_position) * t;
		
		for (uint i = 0; i < 3; i++) func_value_local[i] += elem(v,i);
	}
	
	for (uint i = 0; i < 3; i++) func_value[i] = func_value_local[i];
}

__device__ void calc_function_p(const uint3& num_funcs, const uint2* nuc_contractions, const float3& point_position,
																const float2* factor_ac, uint func_index, float* func_value, float3* Fg)
{
	uint atom_nuc = nuc_contractions[func_index].x;
	float3 atom_nuc_position = gpu_atom_positions[atom_nuc];
	float dist = distance2(point_position, atom_nuc_position);

	uint func_contractions = nuc_contractions[func_index].y;
	
	float func_value_local[3];
	float3 local_Fg[3];
	for (uint i = 0; i < 3; i++) { func_value_local[i] = 0.0f; local_Fg[i] = make_float3(0.0f,0.0f,0.0f); }
	
	float3 txyz = (point_position - atom_nuc_position) * 2.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[func_index * MAX_CONTRACTIONS + contraction];
		float rexp = curr_factor_ac.x * dist;
		if (rexp > 30.0f) continue;

		float t = expf(-rexp) * curr_factor_ac.y;
		float3 term = (point_position - atom_nuc_position) * t;
		
		func_value_local[0] += term.x;
		local_Fg[0] = local_Fg[0] + txyz * term.x * curr_factor_ac.x - make_float3(t, 0.0f,0.0f);
		
		func_value_local[1] += term.y;
		local_Fg[1] = local_Fg[1] + txyz * term.y * curr_factor_ac.x - make_float3(0.0f, t, 0.0f);

		func_value_local[2] += term.z;
		local_Fg[2] = local_Fg[2] + txyz * term.z * curr_factor_ac.x - make_float3(0.0f, 0.0f, t);
	}
	
	for (uint i = 0; i < 3; i++) { func_value[i] = func_value_local[i]; Fg[i] = local_Fg[i]; }
}

// ************** D **************
__device__ void calc_function_d(const uint3& num_funcs, const uint2* nuc_contractions, const float3& point_position,
																const float2* factor_ac, uint func_index, float normalization_factor, float* func_value)
{
	uint atom_nuc = nuc_contractions[func_index].x;
	float3 atom_nuc_position = gpu_atom_positions[atom_nuc];
	float dist = distance2(point_position, atom_nuc_position);

	uint func_contractions = nuc_contractions[func_index].y;
	
	float func_value_local[6];
	for (uint i = 0; i < 6; i++) func_value_local[i] = 0.0f;

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
				func_value_local[ij] += t * t1 * t2;
				ij++;
			}
		}
	}
	
	for (uint i = 0; i < 6; i++) func_value[i] = func_value_local[i];
}

__device__ void calc_function_d(const uint3& num_funcs, const uint2* nuc_contractions, const float3& point_position,
																const float2* factor_ac, uint func_index, float normalization_factor, float* func_value, float3* Fg)
{
	uint atom_nuc = nuc_contractions[func_index].x;
	float3 atom_nuc_position = gpu_atom_positions[atom_nuc];
	float dist = distance2(point_position, atom_nuc_position);

	uint func_contractions = nuc_contractions[func_index].y;
	
	float func_value_local[6];
	float3 local_Fg[6];
	for (uint i = 0; i < 6; i++) { func_value_local[i] = 0.0f; local_Fg[i] = make_float3(0.0f,0.0f,0.0f); }

	float3 t3 = (point_position - atom_nuc_position) * 2.0f;

	for (uint contraction = 0; contraction < func_contractions; contraction++) {
		float2 curr_factor_ac = factor_ac[func_index * MAX_CONTRACTIONS + contraction];
		float rexp = curr_factor_ac.x * dist;
		if (rexp > 30.0f) continue;

		float t0 = expf(-rexp) * curr_factor_ac.y;

		float3 v = point_position - atom_nuc_position;

		uint ij = 0;
		for (uint i = 0; i < 3; i++) {
			float t1 = elem(v, i);
			for (uint j = 0; j <= i; j++) {
				float t = (i == j ? normalization_factor * t0 : t0);
				float t2 = elem(v,j);
				float term = t * t1 * t2;
				func_value_local[ij] += term;
				local_Fg[ij] = local_Fg[ij] + t3 * term * curr_factor_ac.x;
				
				if (i == 0) local_Fg[ij].x = local_Fg[ij].x - t * t2;
				else if (i == 1) local_Fg[ij].y = local_Fg[ij].y - t * t2;
				else if (i == 2) local_Fg[ij].z = local_Fg[ij].z - t * t2;

				if (j == 0) local_Fg[ij].x = local_Fg[ij].x - t * t1;
				else if (j == 1) local_Fg[ij].y = local_Fg[ij].y - t * t1;
				else if (j == 2) local_Fg[ij].z = local_Fg[ij].z - t * t1;

				ij++;
			}
		}
	}
	
	for (uint i = 0; i < 6; i++) { func_value[i] = func_value_local[i]; Fg[i] = local_Fg[i]; }
}

