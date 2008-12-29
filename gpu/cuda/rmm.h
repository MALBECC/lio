
/*
 * Funcion llamada para cada (i,j) en RMM, para calcular RMM(i,j) -> un thread por cada punto
 */

// TODO: esto desperdicia la mitad de los threads -> quizas se puede armar una grilla sin los bloques que no hagan nada
// TODO: medir bien cuanto esta coalesceando y cuanto no

__global__ void gpu_update_rmm(float* factors, uint points, float* rmm, float* function_values, uint m)
{
	uint3 pos = index(blockDim, blockIdx, threadIdx);

	uint i = pos.x; // columna
	uint j = pos.y; // fila
	
	bool valid_thread = (i < m && j < m && i <= j); // quiero triangulo inferior solamente TODO: sacar esto

	uint rmm_idx = (i * m - (i * (i - 1)) / 2) + (j - i);

  // calculate this rmm
	float rmm_local = 0.0f;

  __shared__ float functions_i_local[RMM_BLOCK_SIZE_XY][RMM_BLOCK_SIZE_XY];	// Fi[point][i]
	__shared__ float functions_j_local[RMM_BLOCK_SIZE_XY][RMM_BLOCK_SIZE_XY]; // Fj[point][j]
	__shared__ float factor_local[RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY];			// factor[point]

	for (uint point_base = 0; point_base < points; point_base += (RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY)) {
		uint abs_threadIdx = threadIdx.y * blockDim.x + threadIdx.x;
		
		__syncthreads();

		/* all threads load a point */
		if (point_base + abs_threadIdx < points)
			factor_local[abs_threadIdx] = factors[point_base + abs_threadIdx];
		
		__syncthreads();

		for (uint point_sub = 0; point_sub < (RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY) && (point_base + point_sub < points); point_sub++) {
			uint point = (point_base + point_sub);
			
			/* every RMM_BLOCK_SIZE_XY iterations, Fi and Fj get filled with RMM_BLOCK_SIZE_XY functions, for RMM_BLOCK_SIZE_XY different points */
      if (point_sub % RMM_BLOCK_SIZE_XY == 0) {
				
				__syncthreads();
				
				if (point + threadIdx.y < points) {
					/* every row of the block loads functions for a different point, and every column loads a different function */
	        if (i < m) functions_i_local[threadIdx.y][threadIdx.x] = function_values[(point + threadIdx.y) * COALESCED_DIMENSION(m) + i];
  	      if ((blockDim.y * blockIdx.y + threadIdx.x) < m) functions_j_local[threadIdx.y][threadIdx.x] = function_values[(point + threadIdx.y) * COALESCED_DIMENSION(m) + (blockIdx.y * blockDim.y + threadIdx.x)];
				}
				
  			__syncthreads();				
      }

			float Fi = 0.0f, Fj = 0.0f;		
				
		 	/* fill local variables from local cache */
			if (valid_thread) {
				Fi = functions_i_local[point_sub % RMM_BLOCK_SIZE_XY][threadIdx.x];
				Fj = (i == j ? Fi : functions_j_local[point_sub % RMM_BLOCK_SIZE_XY][threadIdx.y]);
				rmm_local += factor_local[point_sub] * Fi * Fj; 					
			}
		}
	}

  if (valid_thread) rmm[rmm_idx] = rmm_local;
}
