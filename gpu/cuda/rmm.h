/*
 * Funcion llamada para cada (i,j) en RMM, para calcular RMM(i,j) -> un thread por cada punto
 */
template <const uint* const curr_layers, uint grid_n>
__global__ void calc_new_rmm(const float3* atom_positions, const uint* types, const float3* point_positions,
														 const float* wang, const uint atoms_n, uint nco, uint3 num_funcs,
														 const uint* nuc, const uint* contractions, bool normalize, const float* factor_a, const float* factor_c,
														 const float* rmm, float* rmm_output, const float* factors, const float* all_functions)
{
	const uint m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6;

	uint3 abs_idx3d = index(blockDim, blockIdx, threadIdx);

#if 0
	uint abs_idx = abs_idx3d.x * divUp(m, 2) + abs_idx3d.y;

	uint i = (uint)floor((-((-(m + 1.0f) + 0.5f) + sqrtf(((m + 1.0f) - 0.5f) * ((m + 1) - 0.5f) - 2.0f * abs_idx))));
	uint j = abs_idx - ((m + 1) * i - i * (i + 1) / 2);
	_EMU(printf("rmm_idx: %i, i: %i, j: %i m: %i f: %.12e\n", abs_idx, i, j, m, -(m + 1.0f) + 0.5f + sqrtf(((m + 1) - 0.5f) * ((m + 1) - 0.5f) - 2 * abs_idx)));	
	//_EMU(assert(j < i));
	
	if (abs_idx > (m * (m + 1)) / 2) return;
	
	uint rmm_idx = abs_idx;	
#else
	uint i = abs_idx3d.x;
	uint j = abs_idx3d.y;
	if (i >= m || j >= m || i > j) { return; }
	
	//uint point_atom_i = threadIdx.x;
	
	uint rmm_idx = (i * m - (i * (i - 1)) / 2) + (j - i);
	_EMU(printf("rmm i: %i j: %i rmm_idx: %i\n", i, j, rmm_idx));
#endif
	
	dim3 energySize(atoms_n, MAX_LAYERS, grid_n);
	
	float Fi = 0.0f, Fj = 0.0f;
	
	// calculate this rmm
	float rmm_local = 0.0f;

	for (uint atom_i = 0; atom_i < atoms_n; atom_i++) {
		uint atom_i_type = types[atom_i];
		uint atom_i_layers = curr_layers[atom_i_type];

		float3 atom_i_position = atom_positions[atom_i];
		float rm = rm_factor[atom_i_type];
		
		float tmp0 = (PI / (atom_i_layers + 1.0f));
		
		_EMU(printf("layers: %i\n", atom_i_layers));

		for (uint layer_atom_i = 0; layer_atom_i < atom_i_layers; layer_atom_i++) {

			float tmp1 = tmp0 * (layer_atom_i + 1.0f);
			float x = cosf(tmp1);
			float r1 = rm * (1.0f + x) / (1.0f - x);

			for (uint point_atom_i = 0; point_atom_i < grid_n; point_atom_i++) {
				uint factor_idx = index_from3d(energySize, dim3(atom_i, layer_atom_i, point_atom_i));
				float factor = factors[factor_idx];

				float3 point_position = atom_i_position + point_positions[point_atom_i] * r1;

				Fi = all_functions[factor_idx * m + i];
				Fj = (i == j ? Fi : all_functions[factor_idx * m + j]);

				rmm_local += factor * Fi * Fj;
			}			
		}
	}

	rmm_output[rmm_idx] = rmm_local;

	_EMU(printf("rmm value(%i) %.12e\n", rmm_idx, rmm_output[rmm_idx]));		
}
