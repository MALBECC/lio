
/*
 * Funcion llamada para cada (i,j) en RMM, para calcular RMM(i,j) -> un thread por cada punto
 */

// TODO: esto desperdicia la mitad de los threads -> quizas se puede armar una grilla sin los bloques que no hagan nada

template<class scalar_type>
__global__ void gpu_update_rmm(scalar_type* factors, uint points, scalar_type* rmm, scalar_type* function_values, uint m)
{
    // skip computation on blocks outside of the valid triangle
    if (blockIdx.x * blockDim.x > blockIdx.y * blockDim.y) return;

    uint3 pos = index(blockDim, blockIdx, threadIdx);

    uint i = pos.x; // columna
    uint j = pos.y; // fila
    uint first_fi = blockIdx.x * blockDim.x;
    uint first_fj = blockIdx.y * blockDim.y;

    bool valid_thread = (i < m && j < m && i <= j); // quiero triangulo inferior solamente TODO: sacar esto

    // calculate this rmm
    scalar_type rmm_local = 0.0f;

    __shared__ scalar_type functions_i_local[RMM_BLOCK_SIZE_XY][RMM_BLOCK_SIZE_XY+1];	// Fi[point][i]
    __shared__ scalar_type functions_j_local[RMM_BLOCK_SIZE_XY][RMM_BLOCK_SIZE_XY+1]; // Fj[point][j]
    __shared__ scalar_type factor_local[RMM_BLOCK_SIZE_XY * RMM_BLOCK_SIZE_XY];			// factor[point]    // TODO: esto seguramente tiene bank conflicts

    uint inc = RMM_BLOCK_SIZE_XY*RMM_BLOCK_SIZE_XY;
    for (uint point_base = 0; point_base < points; point_base += inc) {
        uint abs_threadIdx = threadIdx.y * blockDim.x + threadIdx.x;  // absolute threadId inside block

        __syncthreads();

        /* all threads load a point */
        if (point_base + abs_threadIdx < points)
            factor_local[abs_threadIdx] = factors[point_base + abs_threadIdx];

        __syncthreads();

        uint last_point = point_base + inc;
        for (uint point = point_base; point < last_point; point += RMM_BLOCK_SIZE_XY) {
            if (point < points) {
                /* every RMM_BLOCK_SIZE_X iterations, Fi and Fj get filled with RMM_BLOCK_SIZE_Y functions, for RMM_BLOCK_SIZE_X different points */
                __syncthreads();

                if (point + threadIdx.x < points) {

                    if ((first_fi + threadIdx.y) < m) {
                        functions_i_local[threadIdx.x][threadIdx.y] = function_values[COALESCED_DIMENSION(points) * (first_fi + threadIdx.y) + (point + threadIdx.x)];
                        functions_i_local[threadIdx.x][threadIdx.y] *= factor_local[(point-point_base) + threadIdx.x];
                    }
                    /* on blocks with some invalid threads, the computation contributes 0 to rmm_local (this avoids instruction serialization) */
                    else functions_i_local[threadIdx.x][threadIdx.y] = 0.0f;

                    if ((first_fj + threadIdx.y) < m) functions_j_local[threadIdx.x][threadIdx.y] = function_values[COALESCED_DIMENSION(points) * (first_fj + threadIdx.y) + (point + threadIdx.x)];
                    else functions_j_local[threadIdx.x][threadIdx.y] = 0.0f;
                }
                else {
                    functions_i_local[threadIdx.x][threadIdx.y] = 0.0f;
                    functions_j_local[threadIdx.x][threadIdx.y] = 0.0f;
                }

                __syncthreads();
#pragma unroll 16
                for (uint point_sub = 0; point_sub < RMM_BLOCK_SIZE_XY; point_sub++) {
                    rmm_local += functions_i_local[point_sub][threadIdx.x] * functions_j_local[point_sub][threadIdx.y];
                }
            }
        }
    }

    if (valid_thread) rmm[COALESCED_DIMENSION(m) * j + i] = rmm_local;
}
