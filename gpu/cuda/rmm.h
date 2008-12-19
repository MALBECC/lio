
/*
 * Funcion llamada para cada (i,j) en RMM, para calcular RMM(i,j) -> un thread por cada punto
 */

#define USE_SHARED 1

// TODO: esto desperdicia la mitad de los threads

__global__ void gpu_update_rmm(float* factors, uint points, float* rmm, float* function_values, uint m)
{
	uint3 pos = index(blockDim, blockIdx, threadIdx);

	uint i = pos.x; // columna
	uint j = pos.y; // fila
	
	bool valid_thread = true;
	if (i >= m || j >= m || i > j) valid_thread = false;	// quiero triangulo inferior solamente TODO: sacar esto

	uint rmm_idx = (i * m - (i * (i - 1)) / 2) + (j - i);
	
	// calculate this rmm
	float rmm_local = 0.0f;
	
	__shared__ float functions_i_local[RMM_BLOCK_SIZE_X];
	__shared__ float functions_j_local[RMM_BLOCK_SIZE_Y];

#if USE_SHARED
	__shared__ float factor_local[RMM_BLOCK_SIZE_X * RMM_BLOCK_SIZE_Y];

	for (uint point_base = 0; point_base < points; point_base += (RMM_BLOCK_SIZE_X * RMM_BLOCK_SIZE_Y)) {
		uint abs_threadIdx = threadIdx.y * blockDim.x + threadIdx.x;

		if (point_base + abs_threadIdx < points)
			factor_local[abs_threadIdx] = factors[point_base + abs_threadIdx];

		for (uint point_sub = 0; point_sub < (RMM_BLOCK_SIZE_X * RMM_BLOCK_SIZE_Y) && (point_base + point_sub < points); point_sub++) {
			uint point = (point_base + point_sub);
			if (threadIdx.y == 0 && i < m) functions_i_local[threadIdx.x] = function_values[point * m + i];
			if (threadIdx.x == 0 && j < m) functions_j_local[threadIdx.y] = function_values[point * m + j];

			__syncthreads();

			float Fi = 0.0f, Fj = 0.0f;
				
		 	/* fill local variables from local cache */
			if (valid_thread) {
				Fi = functions_i_local[threadIdx.x];
				Fj = (i == j ? Fi : functions_j_local[threadIdx.y]);
						
				rmm_local += factor_local[point_sub] * Fi * Fj; 					
			}

			__syncthreads();
		}
	}

#else	
	__shared__ float factor_local;

	for (uint point = 0; point < points; point++) {				
		__syncthreads();				 // por si el escritor se adelanta a los lectores				

		/* cache into local memory */
		if (threadIdx.x == 0 && threadIdx.y == 0) factor_local = factors[point];				
		if (threadIdx.y == 0 && i < m) functions_i_local[threadIdx.x] = function_values[point * m + i];
		if (threadIdx.x == 0 && j < m) functions_j_local[threadIdx.y] = function_values[point * m + j];


		float factor = 0.0f;
		float Fi = 0.0f, Fj = 0.0f;
				
		 /* fill local variables from local cache */

		__syncthreads();	// por si los lectores se adelantan al escritor				

		if (valid_thread) {
			factor = factor_local;
			Fi = functions_i_local[threadIdx.x];
			Fj = (i == j ? Fi : functions_j_local[threadIdx.y]);
					
			/*_EMU(printf("read %.12e %.12e %i %i %i\n", factor, factors[factor_idx], factor_idx, threadIdx.x, threadIdx.y));
				_EMU(printf("read %.12e %.12e %.12e %.12e %i %i %i\n", all_functions[factor_idx * m + i], all_functions[factor_idx * m + j],
							functions_i_local[threadIdx.x], functions_j_local[threadIdx.y],
							factor_idx, threadIdx.x, threadIdx.y));*/ 
			//if (all_forces && j == 0 && i < atoms_n) force = force + dd[factor_idx * m + i] * factor_local;

			rmm_local += factor * Fi * Fj; 					
		}
	}
#endif

	if (valid_thread) rmm[rmm_idx] = rmm_local;
}
