
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

  // calculate this rmm
	float rmm_local = 0.0f;

  __shared__ float functions_i_local[RMM_BLOCK_SIZE_XY+1][RMM_BLOCK_SIZE_XY+1];	// Fi[point][i]
	__shared__ float functions_j_local[RMM_BLOCK_SIZE_XY+1][RMM_BLOCK_SIZE_XY+1]; // Fj[point][j]
	__shared__ float factor_local[RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY];			// factor[point]    // TODO: esto seguramente tiene bank conflicts

	for (uint point_base = 0; point_base < points; point_base += (RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY)) {
		uint abs_threadIdx = threadIdx.y * blockDim.x + threadIdx.x;  // absolute threadId inside block
		
		__syncthreads();

		/* all threads load a point */
		if (point_base + abs_threadIdx < points)
			factor_local[abs_threadIdx] = factors[point_base + abs_threadIdx];
		
		__syncthreads();

    #pragma unroll 16
		for (uint point_sub = 0; point_sub < (RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY); point_sub++) {
      if (point_base + point_sub < points) {
        uint point = (point_base + point_sub);
        uint point_mod = (point_sub % RMM_BLOCK_SIZE_XY);

        /* every RMM_BLOCK_SIZE_X iterations, Fi and Fj get filled with RMM_BLOCK_SIZE_Y functions, for RMM_BLOCK_SIZE_X different points */
        if (point_mod == 0) {

          __syncthreads();

          if (point + threadIdx.x < points) {
            uint first_fi = blockIdx.x * blockDim.x;
            uint first_fj = blockIdx.y * blockDim.y;

            if ((first_fi + threadIdx.y) < m) {
              functions_i_local[threadIdx.x][threadIdx.y] = function_values[COALESCED_DIMENSION(points) * (first_fi + threadIdx.y) + (point + threadIdx.x)];
              functions_i_local[threadIdx.x][threadIdx.y] *= factor_local[point_sub + threadIdx.x];
            }
            else functions_i_local[threadIdx.x][threadIdx.y] = 0.0f;

            if ((first_fj + threadIdx.y) < m) functions_j_local[threadIdx.x][threadIdx.y] = function_values[COALESCED_DIMENSION(points) * (first_fj + threadIdx.y) + (point + threadIdx.x)];
            else functions_j_local[threadIdx.x][threadIdx.y] = 0.0f;
          }
          else {
            functions_i_local[threadIdx.x][threadIdx.y] = 0.0f;
            functions_j_local[threadIdx.x][threadIdx.y] = 0.0f;
          }

          __syncthreads();
        }

        /* fill local variables from local cache */
        /* NOTE: this condition avoids computation on blocks where no thread is valid; on blocks with some valid threads, the computation
         * is still performed but contributes 0 to rmm_local (this avoids instruction serialization) */
        if (blockIdx.x * blockDim.x <= blockIdx.y * blockDim.y) {
          rmm_local += functions_i_local[point_mod][threadIdx.x] * functions_j_local[point_mod][threadIdx.y];
        }
      }
    }
	}

  if (valid_thread) rmm[COALESCED_DIMENSION(m) * j + i] = rmm_local;
}
