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
__global__ void GS_compute_partial(
    uint points, const scalar_type* function_values, uint m,
    const vec_type<scalar_type, 4>* gradient_values,
    scalar_type* out_partial_density, vec_type<scalar_type, 4>* out_dxyz) {
  uint point = blockIdx.x;
  uint i = threadIdx.x + blockIdx.y * 2 * DENSITY_BLOCK_SIZE;
  uint i2 = i + DENSITY_BLOCK_SIZE;
  uint min_i = blockIdx.y * 2 * DENSITY_BLOCK_SIZE + DENSITY_BLOCK_SIZE;
  bool valid_thread = (i < m);
  bool valid_thread2 = (i2 < m);

  // GS density variables
  scalar_type w, w2;
  w = w2 = 0.0f;
  vec_type<scalar_type, 3> w3, w32;
  w3 = w32 = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

  int position = threadIdx.x;

  __shared__ scalar_type fj_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh[DENSITY_BLOCK_SIZE];

  // Si nos vamos a pasar del bloque con el segundo puntero, hacemos que haga la
  // misma cuenta
  if (min_i > m) {
    min_i = min_i - DENSITY_BLOCK_SIZE;
  }

  for (int bj = 0; bj <= min_i; bj += DENSITY_BLOCK_SIZE) {
    __syncthreads();
    if (bj + position < m) {
      fj_sh[position] = function_values[(m)*point + (bj + position)];
      fgj_sh[position] = vec_type<scalar_type, 3>(
          gradient_values[(m)*point + (bj + position)]);
    }
    __syncthreads();
    scalar_type fjreg;
    vec_type<scalar_type, 3> fgjreg;

    if (valid_thread) {
      for (int j = 0; j < DENSITY_BLOCK_SIZE; j++) {
        fjreg = fj_sh[j];
        fgjreg = fgj_sh[j];
        if ((bj + j) <= i) {
          scalar_type rdm_this_thread;

          // GS density
          rdm_this_thread = fetch(rmm_gpu_tex, (float)(bj + j), (float)i);
          w += rdm_this_thread * fjreg;
          w3 += fgjreg * rdm_this_thread;
        }

        if (valid_thread2 && ((bj + j) <= i2)) {
          scalar_type rdm_this_thread2;

          // GS density
          rdm_this_thread2 = fetch(rmm_gpu_tex, (float)(bj + j), (float)i2);
          w2 += rdm_this_thread2 * fjreg;
          w32 += fgjreg * rdm_this_thread2;
        }
      }  // loop j
    }    // valid tread
  }      // loop bj

  // GS density
  scalar_type partial_density(0.0f);
  vec_type<scalar_type, 3> dxyz = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

  if (valid_thread) {
    scalar_type Fi = function_values[(m)*point + i];
    vec_type<scalar_type, 3> Fgi;
    Fgi = vec_type<scalar_type, 3>(gradient_values[(m)*point + i]);
    // GS density
    partial_density = Fi * w;
    dxyz = Fgi * w + w3 * Fi;

    if (valid_thread2) {
      scalar_type Fi2 = function_values[(m)*point + i2];
      vec_type<scalar_type, 3> Fgi2, Fhi12, Fhi22;
      Fgi2 = vec_type<scalar_type, 3>(gradient_values[(m)*point + i2]);
      // GS density
      partial_density += Fi2 * w2;
      dxyz += Fgi2 * w2 + w32 * Fi2;
    }  // valid thread 2
  }    // valid thread

  __syncthreads();
  // GS density
  fj_sh[position] = partial_density;
  fgj_sh[position] = dxyz;
  __syncthreads();

  for (int j = 2; j <= DENSITY_BLOCK_SIZE; j = j * 2) {
    int index = position + DENSITY_BLOCK_SIZE / j;
    if (position < DENSITY_BLOCK_SIZE / j) {
      // GS density
      fj_sh[position] += fj_sh[index];
      fgj_sh[position] += fgj_sh[index];
    }
  }

  if (threadIdx.x == 0) {
    const int myPoint = blockIdx.y * points + blockIdx.x;
    // GS density
    out_partial_density[myPoint] = fj_sh[position];
    out_dxyz[myPoint] = vec_type<scalar_type, 4>(fgj_sh[position]);
  }
}
