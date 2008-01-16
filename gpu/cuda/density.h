/*--------------------------------------------------- Local Density Functionals -------------------------------------------------------------*/
__device__ void local_density_kernel(float& density, uint3 num_funcs, const uint* nuc, const uint* contractions, float3 point_position,
																		 const float3* atom_positions, const float3* atom_positions_shared, bool normalize, const float2* factor_ac,
																		 /*const float* factor_c, */const float* rmm, uint nco, uint big_index, float* F, uint Ndens)
{
	const uint& funcs_s = num_funcs.x;
	const uint& funcs_p = num_funcs.y;
	const uint& funcs_d = num_funcs.z;
	const uint& m = funcs_s + funcs_p * 3 + funcs_d * 6;
	uint func_real = 0;
	
	density = 0.0f;	

	/* s functions */
	for (uint func = 0; func < funcs_s; func++, func_real++)
		calc_function_s(num_funcs, nuc, contractions, point_position, atom_positions, atom_positions_shared, factor_ac, /*factor_c,*/ func, &F[func_real]);
	
	/* p functions */
	for (uint func = funcs_s; func <  funcs_s + funcs_p; func++, func_real+=3)
		calc_function_p(num_funcs, nuc, contractions, point_position, atom_positions, atom_positions_shared, factor_ac, /*factor_c,*/ func, &F[func_real]);
	
	/* d functions */
	float normalization_factor = (normalize ? rsqrtf(3.0f) : 1.0f);
	
	for (uint func = (funcs_s + funcs_p); func < (funcs_s + funcs_p + funcs_d); func++, func_real+=6)
		calc_function_d(num_funcs, nuc, contractions, point_position, atom_positions, atom_positions_shared, factor_ac, /*factor_c, */func, normalization_factor, &F[func_real]);
	
	/*#ifdef _DEBUG
	for (uint i = 0; i < m; i++)
		_EMU(printf("func %i %.12e\n", big_index, F[i]));
	#endif*/

	/* density */
	if (Ndens == 1) {
		uint k = 0;
		
		for (uint i = 0; i < m; i++) {
			float Fi = F[i];
			if (Fi == 0.0f) {
				k += m - i;
				continue;
			}
			
			float w = 0.0f;
			for (uint j = i; j < m; j++) {
				w += rmm[k] * F[j];
				k++;				
			}
			
			density += Fi * w;
		}
	}
	else {
		uint k = 0;
		for (uint i = 0; i < nco; i++) {
			float w = 0.0f;
			for (uint j = 0; j < m; j++) {
				w += rmm[k] * F[j];
				k++;
			}
			density += w * w;
		}
		density *= 2.0f;		
	}
}


#if 0
/*---------------------------------------------- Density Functionals with Force Calculation --------------------------------------------------------------*/
__device__ void density_deriv_kernel(float& density, uint3 num_funcs, const uint* nuc, const uint* nuc_inverse, const uint* contractions,
																		 float3 point_position, const float3* atom_positions, bool normalize, const float* factor_a,
																		 const float* factor_c, const float* rmm, uint nco, uint big_index, float4* F, uint Ndens,
																		 float* forces)
{
	const uint& funcs_s = num_funcs.x;
	const uint& funcs_p = num_funcs.y;
	const uint& funcs_d = num_funcs.z;
	const uint& m = funcs_s + funcs_p * 3 + funcs_d * 6;
	uint func_real = 0;
	
	density = 0.0f;

	/* s functions */
	for (uint func = 0; func < funcs_s; func++, func_real++)
		calc_function_s(num_funcs, nuc, contractions, point_position, atom_positions, factor_a, factor_c, func, &F[func_real]);
	
	/* p functions */
	for (uint func = funcs_s; func <  funcs_s + funcs_p; func++, func_real+=3)
		calc_function_p(num_funcs, nuc, contractions, point_position, atom_positions, factor_a, factor_c, func, &F[func_real]);
	
	/* d functions */
	float normalization_factor = (normalize ? rsqrtf(3.0f) : 1.0f);
	
	for (uint func = (funcs_s + funcs_p); func < (funcs_s + funcs_p + funcs_d); func++, func_real+=6)
		calc_function_d(num_funcs, nuc, contractions, point_position, atom_positions, factor_a, factor_c, func, normalization_factor, &F[func_real]);
	
	/*#ifdef _DEBUG
	for (uint i = 0; i < m; i++)
		_EMU(printf("func %i %.12e\n", big_index, F[i]));
	#endif*/
		

	// On each float4 of F: f.w is function value; F.x,F.y,F.z is the gradient
		
	/* density */
	if (Ndens == 1) {
		uint k = 0;
		
		for (uint i = 0; i < m; i++) {
			float4 Fi = F[i];
			
			if (Fi.w == 0.0f) {
				k += m - i;
				continue;
			}
			
			float w = 0.0f;
			for (uint j = i; j < m; j++) {
				w += rmm[k] * F[j].w;
				k++;				
			}
			
			density += Fi.w * w;						
		}
		
		/* TODO: agregar lo que agrego abajo (para calcular fuerzas cuando Ndens == 1) */
	}
	else {
		uint k = 0;
		for (uint i = 0; i < nco; i++) {
			float w = 0.0f;
			for (uint j = 0; j < m; j++) {
				w += rmm[k] * F[j];
				k++;
			}
			density += w * w;
		}
		density *= 2.0f;		
	}
}


#endif
