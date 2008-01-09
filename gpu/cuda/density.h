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
			if (F[i] == 0.0f) {
				k += m - i;
				continue;
			}
			
			float w = 0.0f;
			for (uint j = i; j < m; j++) {
				w += rmm[k] * F[j];
				k++;				
			}
			
			density += F[i] * w;
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
/*---------------------------------------------- Non-Local Density Functionals --------------------------------------------------------------*/
__device__ void nonlocal_density_kernel(float& density, uint3 num_funcs, const uint* nuc, const uint* contractions, float3 point_position,
																				const float3* atom_positions, bool normalize, const float* factor_a,
																				const float* factor_c, const float* rmm, uint nco, uint big_index, float* F, uint Ndens,
																				float Dx, float Dy, float Dz,float Dxx, float Dxy, float Dxz, float Dyy, float Dzz, float Dyz)
{
	const uint& funcs_s = num_funcs.x;
	const uint& funcs_p = num_funcs.y;
	const uint& funcs_d = num_funcs.z;
	const uint& m = funcs_s + funcs_p * 3 + funcs_d * 6;
	uint func_real = 0;
	
	density = 0.0f;
	Dx = Dy = Dz = Dxx = Dxy = Dxz = Dyy = Dzz = Dyz = 0.0f;	

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

	/* density */
	if (Ndens == 1) {
		uint k = 0;
		
		for (uint i = 0; i < m; i++) {
			if (F[i] == 0.0f) {
				k += m - i;
				continue;
			}
			
			float w = 0.0f;			
			for (uint j = i; j < m; j++) {
				w += rmm[k] * F[j];
				k++;				
			}
			
			density += F[i] * w;
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


#endif
