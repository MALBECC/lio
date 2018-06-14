#if FULL_DOUBLE
static __inline__ __device__ double fetch_double(texture<int2, 2> t, float x,
                                                 float y) {
  int2 v = tex2D(t, x, y);
  return __hiloint2double(v.y, v.x);
}
#define fetch(t, x, y) fetch_double(t, x, y)
#else
#define fetch(t, x, y) tex2D(t, x, y)
#endif

template <class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_compute_density(
    scalar_type* const energy, scalar_type* const factor,
    const scalar_type* const point_weights, uint points,
    const scalar_type* function_values,
    const vec_type<scalar_type, 4>* gradient_values,
    const vec_type<scalar_type, 4>* hessian_values, uint m,
    scalar_type* out_partial_density, vec_type<scalar_type, 4>* out_dxyz,
    vec_type<scalar_type, 4>* out_dd1, vec_type<scalar_type, 4>* out_dd2) {
  uint point = blockIdx.x;

  uint i = threadIdx.x + blockIdx.y * 2 * DENSITY_BLOCK_SIZE;

  uint i2 = i + DENSITY_BLOCK_SIZE;

  uint min_i = blockIdx.y * 2 * DENSITY_BLOCK_SIZE + DENSITY_BLOCK_SIZE;

  bool valid_thread = (i < m);
  bool valid_thread2 = (i2 < m);

  scalar_type w = 0.0f;
  scalar_type w2 = 0.0f;
  vec_type<scalar_type, 3> w3, ww1, ww2;
  vec_type<scalar_type, 3> w32, ww12, ww22;
  if (!lda) {
    w3 = ww1 = ww2 = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    w32 = ww12 = ww22 = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
  }

  int position = threadIdx.x;

  __shared__ scalar_type fj_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fh1j_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fh2j_sh[DENSITY_BLOCK_SIZE];

  // Si nos vamos a pasar del bloque con el segundo puntero, hacemos que haga la
  // misma cuenta
  if (min_i > m) {
    min_i = min_i - DENSITY_BLOCK_SIZE;
  }

  for (int bj = 0; bj <= min_i; bj += DENSITY_BLOCK_SIZE) {
    __syncthreads();
    if (bj + position < m) {
      fj_sh[position] = function_values[(m)*point + (bj + position)];
      if (!lda) {
        fgj_sh[position] = vec_type<scalar_type, 3>(
            gradient_values[(m)*point + (bj + position)]);

        fh1j_sh[position] = vec_type<scalar_type, 3>(
            hessian_values[(m)*2 * point + (2 * (bj + position) + 0)]);
        fh2j_sh[position] = vec_type<scalar_type, 3>(
            hessian_values[(m)*2 * point + (2 * (bj + position) + 1)]);
      }
    }

    __syncthreads();
    scalar_type fjreg;
    vec_type<scalar_type, 3> fgjreg;
    vec_type<scalar_type, 3> fh1jreg;
    vec_type<scalar_type, 3> fh2jreg;

    if (valid_thread) {
      for (int j = 0; j < DENSITY_BLOCK_SIZE; j++) {
        fjreg = fj_sh[j];

        if (!lda) {
          fgjreg = fgj_sh[j];
          fh1jreg = fh1j_sh[j];
          fh2jreg = fh2j_sh[j];
        }
        // fetch es una macro para tex2D

        if ((bj + j) <= i) {
          scalar_type rdm_this_thread =
              fetch(rmm_input_gpu_tex, (float)(bj + j), (float)i);
          w += rdm_this_thread * fjreg;

          if (!lda) {
            w3 += fgjreg * rdm_this_thread;
            ww1 += fh1jreg * rdm_this_thread;
            ww2 += fh2jreg * rdm_this_thread;
          }
        }

        if (valid_thread2 && ((bj + j) <= i2)) {
          scalar_type rdm_this_thread2 =
              fetch(rmm_input_gpu_tex, (float)(bj + j), (float)i2);
          w2 += rdm_this_thread2 * fjreg;

          if (!lda) {
            w32 += fgjreg * rdm_this_thread2;
            ww12 += fh1jreg * rdm_this_thread2;
            ww22 += fh2jreg * rdm_this_thread2;
          }
        }
      }
    }
  }

  scalar_type partial_density(0.0f);
  vec_type<scalar_type, 3> dxyz, dd1, dd2;
  dxyz = dd1 = dd2 = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

  if (valid_thread) {
    scalar_type Fi = function_values[(m)*point + i];
    vec_type<scalar_type, 3> Fgi, Fhi1, Fhi2;
    if (!lda) {
      Fgi = vec_type<scalar_type, 3>(gradient_values[(m)*point + i]);

      Fhi1 =
          vec_type<scalar_type, 3>(hessian_values[(m)*2 * point + (2 * i + 0)]);
      Fhi2 =
          vec_type<scalar_type, 3>(hessian_values[(m)*2 * point + (2 * i + 1)]);
    }

    partial_density = Fi * w;
    if (!lda) {
      dxyz = Fgi * w + w3 * Fi;
      dd1 = Fgi * w3 * 2.0f + Fhi1 * w + ww1 * Fi;

      vec_type<scalar_type, 3> FgXXY(Fgi.x, Fgi.x, Fgi.y);
      vec_type<scalar_type, 3> w3YZZ(w3.y, w3.z, w3.z);
      vec_type<scalar_type, 3> FgiYZZ(Fgi.y, Fgi.z, Fgi.z);
      vec_type<scalar_type, 3> w3XXY(w3.x, w3.x, w3.y);

      dd2 = FgXXY * w3YZZ + FgiYZZ * w3XXY + Fhi2 * w + ww2 * Fi;
    }

    if (valid_thread2) {
      scalar_type Fi2 = function_values[(m)*point + i2];
      vec_type<scalar_type, 3> Fgi2, Fhi12, Fhi22;
      if (!lda) {
        Fgi2 = vec_type<scalar_type, 3>(gradient_values[(m)*point + i2]);

        Fhi12 = vec_type<scalar_type, 3>(
            hessian_values[(m)*2 * point + (2 * i2 + 0)]);
        Fhi22 = vec_type<scalar_type, 3>(
            hessian_values[(m)*2 * point + (2 * i2 + 1)]);
      }

      partial_density += Fi2 * w2;
      if (!lda) {
        dxyz += Fgi2 * w2 + w32 * Fi2;
        dd1 += Fgi2 * w32 * 2.0f + Fhi12 * w2 + ww12 * Fi2;

        vec_type<scalar_type, 3> FgXXY(Fgi2.x, Fgi2.x, Fgi2.y);
        vec_type<scalar_type, 3> w3YZZ(w32.y, w32.z, w32.z);
        vec_type<scalar_type, 3> FgiYZZ(Fgi2.y, Fgi2.z, Fgi2.z);
        vec_type<scalar_type, 3> w3XXY(w32.x, w32.x, w32.y);

        dd2 += FgXXY * w3YZZ + FgiYZZ * w3XXY + Fhi22 * w2 + ww22 * Fi2;
      }
    }
  }

  __syncthreads();
  // Estamos reutilizando la memoria shared por block para hacer el acumulado
  // por block.
  // No hace falta poner en cero porque si no es valid_thread, ya estan en cero
  fj_sh[position] = partial_density;
  fgj_sh[position] = dxyz;
  fh1j_sh[position] = dd1;
  fh2j_sh[position] = dd2;

  __syncthreads();

  for (int j = 2; j <= DENSITY_BLOCK_SIZE; j = j * 2) {
    int index = position + DENSITY_BLOCK_SIZE / j;
    if (position < DENSITY_BLOCK_SIZE / j) {
      fj_sh[position] += fj_sh[index];
      fgj_sh[position] += fgj_sh[index];
      fh1j_sh[position] += fh1j_sh[index];
      fh2j_sh[position] += fh2j_sh[index];
    }
  }

  if (threadIdx.x == 0) {
    const int myPoint = blockIdx.y * points + blockIdx.x;
    out_partial_density[myPoint] = fj_sh[position];
    // printf("%.4e ",out_partial_density);
    out_dxyz[myPoint] = vec_type<scalar_type, 4>(fgj_sh[position]);
    out_dd1[myPoint] = vec_type<scalar_type, 4>(fh1j_sh[position]);
    out_dd2[myPoint] = vec_type<scalar_type, 4>(fh2j_sh[position]);
  }
}
