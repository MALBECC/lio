
/*
 * Funcion llamada para cada (i,j) en RMM, para calcular RMM(i,j) -> un thread por cada punto
 */
template <uint grid_n, const uint* const curr_layers>
__global__ void calc_new_rmm(const float3* atom_positions, const uint* types, const float3* point_positions,
														 const uint atoms_n, uint nco, uint3 num_funcs,
														 const uint* nuc, const uint* contractions, bool normalize, const float2* factor_ac,
														 const float* rmm, float* rmm_output, const float* factors, const float* all_functions)
{
	const uint m = num_funcs.x + num_funcs.y * 3 + num_funcs.z * 6;

	uint3 abs_idx3d = index(blockDim, blockIdx, threadIdx);

#if 0
	uint abs_idx = abs_idx3d.x * divUp(m, 2) + abs_idx3d.y;
	uint rmm_idx = abs_idx;

	uint i = (uint)floor((-((-(m + 1.0f) + 0.5f) + sqrtf(((m + 1.0f) - 0.5f) * ((m + 1) - 0.5f) - 2.0f * abs_idx))));
	uint j = abs_idx - ((m + 1) * i - i * (i + 1) / 2);	
#else
	uint i = abs_idx3d.x; // columna
	uint j = abs_idx3d.y; // fila
	
	uint rmm_idx = (i * m - (i * (i - 1)) / 2) + (j - i);
#endif
	
	bool valid_thread = true;
	if (i >= m || j >= m || i > j) valid_thread = false;	// quiero triangulo inferior solamente
	
	_EMU(printf("rmm i: %i j: %i rmm_idx: %i\n", i, j, rmm_idx));
	
	dim3 energySize(atoms_n, MAX_LAYERS, grid_n);
	
	// calculate this rmm
	float rmm_local = 0.0f;
	
	__shared__ float factor_local;
	__shared__ float functions_i_local[RMM_BLOCK_SIZE_X];
	__shared__ float functions_j_local[RMM_BLOCK_SIZE_Y];

	for (uint atom_i = 0; atom_i < atoms_n; atom_i++) {
		uint atom_i_type = types[atom_i];
		uint atom_i_layers = curr_layers[atom_i_type];

		_EMU(printf("layers: %i\n", atom_i_layers));

		for (uint layer_atom_i = 0; layer_atom_i < atom_i_layers; layer_atom_i++) {
			for (uint point_atom_i = 0; point_atom_i < grid_n; point_atom_i++) {
				uint factor_idx = index_from3d(energySize, dim3(atom_i, layer_atom_i, point_atom_i));
				
				__syncthreads();				 // por si el escritor se adelanta a los lectores				

				/* cache into local memory */
				if (threadIdx.x == 0 && threadIdx.y == 0) {
					factor_local = factors[factor_idx];
					_EMU(printf("load %.12e %.12e %i %i %i\n", factor_local, factors[factor_idx], factor_idx, threadIdx.x, threadIdx.y));
				}
				
				if (threadIdx.y == 0) { functions_i_local[threadIdx.x] = all_functions[factor_idx * m + i]; }
				if (threadIdx.x == 0) { functions_j_local[threadIdx.y] = all_functions[factor_idx * m + j]; }

				float factor = 0.0f;
				float Fi = 0.0f, Fj = 0.0f;
				
			  /* fill local variables from local cache */
				__syncthreads();	// por si los lectores se adelantan al escritor				
				if (valid_thread) {
					factor = factor_local;
					Fi = functions_i_local[threadIdx.x];
					Fj = (i == j ? Fi : functions_j_local[threadIdx.y]);
					
					_EMU(printf("read %.12e %.12e %i %i %i\n", factor, factors[factor_idx], factor_idx, threadIdx.x, threadIdx.y));
					_EMU(printf("read %.12e %.12e %.12e %.12e %i %i %i\n", all_functions[factor_idx * m + i], all_functions[factor_idx * m + j],
											functions_i_local[threadIdx.x], functions_j_local[threadIdx.y],
											factor_idx, threadIdx.x, threadIdx.y));					
				}
				
				/* compute */
				if (valid_thread) { rmm_local += factor * Fi * Fj; }				
			}			
		}
	}

	if (valid_thread) {
		rmm_output[rmm_idx] = rmm_local;
		_EMU(printf("rmm value(%i) %.12e\n", rmm_idx, rmm_output[rmm_idx]));
	}
}
