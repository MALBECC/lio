// TODO: si se juntara con energy.h (teniendo un if templatizado tipo do_forces, que hace que solo se desde donde se debe (j >= i)) se leeria RMM una sola vez

template<class scalar_type>
__global__ void gpu_compute_density_derivs(scalar_type* function_values, vec_type<scalar_type,4>* gradient_values,
                                           uint* nuc, vec_type<scalar_type,4>* density_deriv, uint points, uint m, uint nuc_count)
{
  uint point = index_x(blockDim, blockIdx, threadIdx);
  bool valid_thread = (point < points);

  __shared__ scalar_type rdm_sh[DENSITY_DERIV_BATCH_SIZE];
  __shared__ uint nuc_sh[DENSITY_DERIV_BATCH_SIZE2];

  for (uint bi = 0; bi < m; bi += DENSITY_DERIV_BATCH_SIZE2) {
    __syncthreads();
    if (threadIdx.x < DENSITY_DERIV_BATCH_SIZE2) {
      if (bi + threadIdx.x < m) nuc_sh[threadIdx.x] = nuc[bi + threadIdx.x];
    }
    __syncthreads();

    for (uint i = 0; i < DENSITY_DERIV_BATCH_SIZE2 && (i + bi) < m; i++) {
      scalar_type w = 0.0f;

      for (uint bj = 0; bj < m; bj += DENSITY_DERIV_BATCH_SIZE) {
        __syncthreads();
        if (threadIdx.x < DENSITY_DERIV_BATCH_SIZE) {
            //fetch es una macro para tex2D definida en energy.h
          if (bj + threadIdx.x < m)
            rdm_sh[threadIdx.x] = fetch(rmm_input_gpu_tex, (float)(bi+i), (float)(bj+threadIdx.x));
          else
            rdm_sh[threadIdx.x] = 0.0f;
        }
        __syncthreads();

        if (valid_thread) {
          for (uint j = 0; j < DENSITY_DERIV_BATCH_SIZE && (bj + j) < m; j++) {
            scalar_type fj = function_values[COALESCED_DIMENSION(points) * (bj + j) + point];
            w += rdm_sh[j] * fj * ((bi + i) == (bj + j) ? 2 : 1);
          }
        }
      }

      if (valid_thread) {
        vec_type<scalar_type,4> Fgi = gradient_values[COALESCED_DIMENSION(points) * (bi + i) + point];
        uint nuci = nuc_sh[i];
        density_deriv[COALESCED_DIMENSION(points) * nuci + point] -= Fgi * w;
        // TODO: esto accede demasiado a memoria, quizas se pueda acumular en sh_mem
      }
    }
  }
}

//===================================================================================================================
template<class scalar_type>
__global__ void gpu_compute_density_derivs_open(scalar_type* function_values, vec_type<scalar_type,4>* gradient_values,
                                           uint* nuc, vec_type<scalar_type,4>* density_deriv_a, vec_type<scalar_type,4>* density_deriv_b,
                                           uint points, uint m, uint nuc_count)
{
  uint point = index_x(blockDim, blockIdx, threadIdx);
  bool valid_thread = (point < points);

  __shared__ scalar_type rdm_a_sh[DENSITY_DERIV_BATCH_SIZE];
  __shared__ scalar_type rdm_b_sh[DENSITY_DERIV_BATCH_SIZE];
  __shared__ uint nuc_sh[DENSITY_DERIV_BATCH_SIZE2];

  for (uint bi = 0; bi < m; bi += DENSITY_DERIV_BATCH_SIZE2) {
    __syncthreads();
    if (threadIdx.x < DENSITY_DERIV_BATCH_SIZE2) {
      if (bi + threadIdx.x < m) nuc_sh[threadIdx.x] = nuc[bi + threadIdx.x];
    }
    __syncthreads();

    for (uint i = 0; i < DENSITY_DERIV_BATCH_SIZE2 && (i + bi) < m; i++) {
      vec_type<scalar_type,4> Fgi;
      if (valid_thread) Fgi = gradient_values[COALESCED_DIMENSION(points) * (bi + i) + point];
      scalar_type w_a = 0.0f;
      scalar_type w_b = 0.0f;

      for (uint bj = 0; bj < m; bj += DENSITY_DERIV_BATCH_SIZE) {
        __syncthreads();
        if (threadIdx.x < DENSITY_DERIV_BATCH_SIZE) {
            //fetch es una macro para tex2D definida en energy.h
            // scalar_type rmd_local = fetch(rmm_input_gpu_tex, (float)(bi+i), (float)(bj+threadIdx.x));
          if (bj + threadIdx.x < m) {
            rdm_a_sh[threadIdx.x] = fetch(rmm_input_gpu_tex, (float)(bi+i), (float)(bj+threadIdx.x));
            rdm_b_sh[threadIdx.x] = fetch(rmm_input_gpu_tex2, (float)(bi+i), (float)(bj+threadIdx.x));
          }    //rdm[COALESCED_DIMENSION(m) * (bi + i) + (bj + threadIdx.x)];
          else {
            rdm_a_sh[threadIdx.x] = 0.0f;
            rdm_b_sh[threadIdx.x] = 0.0f;
          }
	      }
        __syncthreads();

        if (valid_thread) {
          for (uint j = 0; j < DENSITY_DERIV_BATCH_SIZE && (bj + j) < m; j++) {
            scalar_type fj = function_values[COALESCED_DIMENSION(points) * (bj + j) + point];
            w_a += rdm_a_sh[j] * fj * ((bi + i) == (bj + j) ? 2 : 1);
            w_b += rdm_b_sh[j] * fj * ((bi + i) == (bj + j) ? 2 : 1);
          }
        }
      }

      if (valid_thread) {
        uint nuci = nuc_sh[i];
        // TODO: esto accede demasiado a memoria, quizas se pueda acumular en sh_mem
        density_deriv_a[COALESCED_DIMENSION(points) * nuci + point] -= Fgi * w_a;
        density_deriv_b[COALESCED_DIMENSION(points) * nuci + point] -= Fgi * w_b;
      }
    }
  }
}

