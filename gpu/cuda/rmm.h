
/*
 * Funcion llamada para cada (i,j) en RMM, para calcular RMM(i,j) -> un thread por cada punto
 */

// TODO: esto desperdicia la mitad de los threads -> quizas se puede armar una grilla sin los bloques que no hagan nada

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
      uint point_mod = (point_sub % RMM_BLOCK_SIZE_XY);
			
			/* every RMM_BLOCK_SIZE_XY iterations, Fi and Fj get filled with RMM_BLOCK_SIZE_XY functions, for RMM_BLOCK_SIZE_XY different points */
      if (point_mod == 0) {
				
				__syncthreads();
				
				if (point + threadIdx.x < points) {
          if ((blockIdx.x * blockDim.x + threadIdx.y) < m) {
            functions_i_local[threadIdx.x][threadIdx.y] = function_values[COALESCED_DIMENSION(points) * (blockIdx.x * blockDim.x + threadIdx.y) + (point + threadIdx.x)];
            functions_i_local[threadIdx.x][threadIdx.y] *= factor_local[point_sub + threadIdx.x];
          }
          else functions_i_local[threadIdx.x][threadIdx.y] = 0.0f;

          if ((blockIdx.y * blockDim.y + threadIdx.y) < m) functions_j_local[threadIdx.x][threadIdx.y] = function_values[COALESCED_DIMENSION(points) * (blockIdx.y * blockDim.y + threadIdx.y) + (point + threadIdx.x)];
          else functions_j_local[threadIdx.x][threadIdx.y] = 0.0f;
				}
        else {
          functions_i_local[threadIdx.x][threadIdx.y] = 0.0f;
          functions_j_local[threadIdx.x][threadIdx.y] = 0.0f;
        }
				
  			__syncthreads();				
      }

			float Fi = 0.0f, Fj = 0.0f;		
				
		 	/* fill local variables from local cache */
      /* NOTE: this condition avoids computation on blocks where no thread is valid, on blocks with some valid threads, the computation
       * is still performed but contributes 0 to rmm_local (this avoids instruction serialization) */
      if (blockIdx.x * blockDim.x <= blockIdx.y * blockDim.y) {
				Fi = functions_i_local[point_mod][threadIdx.x];
        Fj = functions_j_local[point_mod][threadIdx.y];
        rmm_local += Fi * Fj;
			}
		}
	}

  // TODO: coalescear (escribir usando dos indices)
  //if (valid_thread) rmm[rmm_idx] = rmm_local;
  if (valid_thread) rmm[COALESCED_DIMENSION(m) * j + i] = rmm_local;
}
