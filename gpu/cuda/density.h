/*--------------------------------------------------- Local Density Functionals -------------------------------------------------------------*/
__device__ void local_density_kernel(float& density, uint3 num_funcs, const uint* nuc, const uint* contractions, float3 point_position,
																		 const float3* atom_positions, const float3* atom_positions_shared, bool normalize, const float2* factor_ac,
																		 const float* rmm, uint nco, uint big_index, float* F, uint Ndens)
{
	const uint& funcs_s = num_funcs.x;
	const uint& funcs_p = num_funcs.y;
	const uint& funcs_d = num_funcs.z;
	const uint& m = funcs_s + funcs_p * 3 + funcs_d * 6;
	uint func_real = 0;
	
	density = 0.0f;

	/* s functions */
	for (uint func = 0; func < funcs_s; func++, func_real++)
		calc_function_s(num_funcs, nuc, contractions, point_position, atom_positions, atom_positions_shared, factor_ac, func, &F[func_real]);
	
	/* p functions */
	for (uint func = funcs_s; func <  funcs_s + funcs_p; func++, func_real+=3)
		calc_function_p(num_funcs, nuc, contractions, point_position, atom_positions, atom_positions_shared, factor_ac, func, &F[func_real]);
	
	/* d functions */
	float normalization_factor = (normalize ? rsqrtf(3.0f) : 1.0f);
	
	for (uint func = (funcs_s + funcs_p); func < (funcs_s + funcs_p + funcs_d); func++, func_real+=6)
		calc_function_d(num_funcs, nuc, contractions, point_position, atom_positions, atom_positions_shared, factor_ac, func, normalization_factor, &F[func_real]);
	
	/* density */
	if (Ndens == 1) {
		uint k = 0;
		
		for (uint i = 0; i < m; i++) {
			float w = 0.0f;
			for (uint j = i; j < m; j++) {
				w += rmm[k] * F[j];
				k++;				
			}
			
			density += F[i] * w;
		}
	}
	else {
		for (uint i = 0; i < nco; i++) {
			float w = 0.0f;
			for (uint j = 0; j < m; j++) {
				w += rmm[i * m + j] * F[j];
			}
			density += w * w;
		}
		density *= 2.0f;
	}
}

/*---------------------------------------------- Density Functionals with Force Calculation --------------------------------------------------------------*/
__device__ void density_deriv_kernel(float& density, uint3 num_funcs, const uint* nuc, const uint* contractions, float3 point_position,
																		 const float3* atom_positions, const float3* atom_positions_shared, bool normalize, const float2* factor_ac,
																		 const float* rmm, uint nco, uint big_index, float* F, uint Ndens, float3* dd, float3* Fg, float3* w3, uint atoms_n)
{
	const uint& funcs_s = num_funcs.x;
	const uint& funcs_p = num_funcs.y;
	const uint& funcs_d = num_funcs.z;
	const uint& m = funcs_s + funcs_p * 3 + funcs_d * 6;
	uint func_real = 0;
	
	density = 0.0f;
	
	/* s functions */
	for (uint func = 0; func < funcs_s; func++, func_real++)
		calc_function_s(num_funcs, nuc, contractions, point_position, atom_positions, atom_positions_shared, factor_ac, func, &F[func_real], &Fg[func_real]);
	
	/* p functions */
	for (uint func = funcs_s; func <  funcs_s + funcs_p; func++, func_real+=3)
		calc_function_p(num_funcs, nuc, contractions, point_position, atom_positions, atom_positions_shared, factor_ac, func, &F[func_real], &Fg[func_real]);
	
	/* d functions */
	float normalization_factor = (normalize ? rsqrtf(3.0f) : 1.0f);
	
	for (uint func = (funcs_s + funcs_p); func < (funcs_s + funcs_p + funcs_d); func++, func_real+=6)
		calc_function_d(num_funcs, nuc, contractions, point_position, atom_positions, atom_positions_shared, factor_ac, func, normalization_factor, &F[func_real], &Fg[func_real]);
	
	for (uint j = 0; j < atoms_n; j++) {
		dd[j] = make_float3(0.0f,0.0f,0.0f);
	}
	

	// density (Ndens is never 1 here)
	uint k = 0; 
	for (uint i = 0; i < nco; i++) {
		float w = 0.0f;
		for (uint atom_i = 0; atom_i < atoms_n; atom_i++) { w3[atom_i] = make_float3(0.0f,0.0f,0.0f); }
														
		func_real = 0;
		for (uint func = 0; func < funcs_s; func++, func_real++, k++) {
			float this_rmm = rmm[k];
			uint atom_nuc = nuc[func];
			w += this_rmm * F[func_real];			
			w3[atom_nuc] = w3[atom_nuc] + Fg[func_real] * this_rmm;
		}		
		for (uint func = funcs_s; func < funcs_s + funcs_p; func++) {
			uint atom_nuc = nuc[func];
			for (uint subfunc = 0; subfunc < 3; subfunc++, func_real++, k++) {
				float this_rmm = rmm[k];
				w += this_rmm * F[func_real];				
				w3[atom_nuc] = w3[atom_nuc] + Fg[func_real] * this_rmm;
			}
		}		
		for (uint func = (funcs_s + funcs_p); func < (funcs_s + funcs_p + funcs_d); func++) {
			uint atom_nuc = nuc[func];
			for (uint subfunc = 0; subfunc < 6; subfunc++, func_real++, k++) {
				float this_rmm = rmm[k];
				w += this_rmm * F[func_real];				
				w3[atom_nuc] = w3[atom_nuc] + Fg[func_real] * this_rmm;
			}
		}
		// TODO: ver si el ciclo se puede hacer por cada atomo
		
		for (uint j = 0; j < atoms_n; j++) {
			dd[j] = dd[j] + w3[j] * w;
		}

		density += w * w;
	}
	density *= 2.0f;
}

