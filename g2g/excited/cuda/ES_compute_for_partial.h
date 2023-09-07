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
__global__ void ES_compute_for_partial(const scalar_type* const point_weights,uint points,
                                const scalar_type* function_values, uint m,
                                const vec_type<scalar_type, 4>* gradient_values,
                                scalar_type* out_partial_tred, vec_type<scalar_type, 4>* out_tredxyz,
                                scalar_type* out_partial_diff, vec_type<scalar_type, 4>* out_diffxyz)
{
  uint point = blockIdx.x;
  uint i = threadIdx.x + blockIdx.y * 2 * DENSITY_BLOCK_SIZE;
  uint i2 = i + DENSITY_BLOCK_SIZE;
  uint min_i = blockIdx.y * 2 * DENSITY_BLOCK_SIZE + DENSITY_BLOCK_SIZE;
  bool valid_thread = (i < m);
  bool valid_thread2 = (i2 < m);

// Transition density variables
  scalar_type z, z2; z = z2 = 0.0f;
  vec_type<scalar_type, 3> z3, z32; z3 = z32 = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

// Difference density variables
  scalar_type x, x2; x = x2 = 0.0f;
  vec_type<scalar_type, 3> x3, x32; x3 = x32 = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

  int position = threadIdx.x;

  __shared__ scalar_type fj_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh[DENSITY_BLOCK_SIZE];

  __shared__ scalar_type fj_sh_tred[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh_tred[DENSITY_BLOCK_SIZE];

  __shared__ scalar_type fj_sh_diff[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh_diff[DENSITY_BLOCK_SIZE];

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

          // Transition density
          rdm_this_thread = fetch(tred_gpu_for_tex, (float)(bj + j), (float)i);
          z  += rdm_this_thread * fjreg;
          z3 += fgjreg * rdm_this_thread;

          // Difference density
          rdm_this_thread = fetch(diff_gpu_for_tex, (float)(bj + j), (float)i);
          x  += rdm_this_thread * fjreg;
          x3 += fgjreg * rdm_this_thread;
        }

        if (valid_thread2 && ((bj + j) <= i2)) {
          scalar_type rdm_this_thread2;

          // Transition density
          rdm_this_thread2 = fetch(tred_gpu_for_tex, (float)(bj + j), (float)i2);
          z2  += rdm_this_thread2 * fjreg;
          z32 += fgjreg * rdm_this_thread2;

          // Difference density
          rdm_this_thread2 = fetch(diff_gpu_for_tex, (float)(bj + j), (float)i2);
          x2  += rdm_this_thread2 * fjreg;
          x32 += fgjreg * rdm_this_thread2;
        }
      } // loop j
    } // valid tread
  } // loop bj

  // Transition density
  scalar_type partial_tred(0.0f);
  vec_type<scalar_type, 3> tredxyz = vec_type<scalar_type, 3>(0.0f,0.0f,0.0f);

  // Difference density
  scalar_type partial_diff(0.0f);
  vec_type<scalar_type, 3> diffxyz = vec_type<scalar_type, 3>(0.0f,0.0f,0.0f);


  if (valid_thread) {
    scalar_type Fi = function_values[(m)*point + i];
    vec_type<scalar_type, 3> Fgi;
    Fgi = vec_type<scalar_type, 3>(gradient_values[(m)*point + i]);
    // Transition density
    partial_tred = Fi * z;
    tredxyz = Fgi * z + z3 * Fi;
    // Difference density
    partial_diff = Fi * x;
    diffxyz = Fgi * x + x3 * Fi;

    if (valid_thread2) {
      scalar_type Fi2 = function_values[(m)*point + i2];
      vec_type<scalar_type, 3> Fgi2, Fhi12, Fhi22;
      Fgi2 = vec_type<scalar_type, 3>(gradient_values[(m)*point + i2]);
      // Transition density
      partial_tred += Fi2 * z2;
      tredxyz += Fgi2 * z2 + z32 * Fi2;
      // Difference density
      partial_diff += Fi2 * x2;
      diffxyz += Fgi2 * x2 + x32 * Fi2;
    } // valid thread 2
  } // valid thread

  __syncthreads();
// Transition density
  fj_sh_tred[position] = partial_tred;
  fgj_sh_tred[position] = tredxyz;
// Difference density
  fj_sh_diff[position] = partial_diff;
  fgj_sh_diff[position] = diffxyz;
  __syncthreads();

  for (int j = 2; j <= DENSITY_BLOCK_SIZE; j = j * 2) {
    int index = position + DENSITY_BLOCK_SIZE / j;
    if (position < DENSITY_BLOCK_SIZE / j) {
      // Transition density
      fj_sh_tred[position] += fj_sh_tred[index];
      fgj_sh_tred[position] += fgj_sh_tred[index];

      // Difference density
      fj_sh_diff[position] += fj_sh_diff[index];
      fgj_sh_diff[position] += fgj_sh_diff[index];
    }
  }

  if (threadIdx.x == 0) {
    const int myPoint = blockIdx.y * points + blockIdx.x;
    // Transition density
    out_partial_tred[myPoint] = fj_sh_tred[position];
    out_tredxyz[myPoint] = vec_type<scalar_type, 4>(fgj_sh_tred[position]);
    // Difference density
    out_partial_diff[myPoint] = fj_sh_diff[position];
    out_diffxyz[myPoint] = vec_type<scalar_type, 4>(fgj_sh_diff[position]);
  }
}

template <class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void ES_compute_for_partial(const scalar_type* const point_weights,uint points,
                                const scalar_type* function_values, uint m,
                                const vec_type<scalar_type, 4>* gradient_values,
                                scalar_type* out_partial_tred, vec_type<scalar_type, 4>* out_tredxyz,
                                scalar_type* out_partial_diff, vec_type<scalar_type, 4>* out_diffxyz,
                                scalar_type* out_partial_density, vec_type<scalar_type, 4>* out_dxyz)
{

  uint point = blockIdx.x;
  uint i = threadIdx.x + blockIdx.y * 2 * DENSITY_BLOCK_SIZE;
  uint i2 = i + DENSITY_BLOCK_SIZE;
  uint min_i = blockIdx.y * 2 * DENSITY_BLOCK_SIZE + DENSITY_BLOCK_SIZE;
  bool valid_thread = (i < m);
  bool valid_thread2 = (i2 < m);

// GS density variables
  scalar_type w, w2; w = w2 = 0.0f;
  vec_type<scalar_type, 3> w3, w32; w3 = w32 = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

// Transition density variables
  scalar_type z, z2; z = z2 = 0.0f;
  vec_type<scalar_type, 3> z3, z32; z3 = z32 = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

// Difference density variables
  scalar_type x, x2; x = x2 = 0.0f;
  vec_type<scalar_type, 3> x3, x32; x3 = x32 = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

  int position = threadIdx.x;

// GS density
  __shared__ scalar_type fj_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh[DENSITY_BLOCK_SIZE];

// Transition density
  __shared__ scalar_type fj_sh_tred[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh_tred[DENSITY_BLOCK_SIZE];

// Difference density
  __shared__ scalar_type fj_sh_diff[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh_diff[DENSITY_BLOCK_SIZE];

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
          rdm_this_thread = fetch(rmm_gpu_for_tex, (float)(bj + j), (float)i);
          w  += rdm_this_thread * fjreg;
          w3 += fgjreg * rdm_this_thread;

          // Transition density
          rdm_this_thread = fetch(tred_gpu_for_tex, (float)(bj + j), (float)i);
          z  += rdm_this_thread * fjreg;
          z3 += fgjreg * rdm_this_thread;

          // Difference density
          rdm_this_thread = fetch(diff_gpu_for_tex, (float)(bj + j), (float)i);
          x  += rdm_this_thread * fjreg;
          x3 += fgjreg * rdm_this_thread;
        }

        if (valid_thread2 && ((bj + j) <= i2)) {
          scalar_type rdm_this_thread2;

          // GS density
          rdm_this_thread2 = fetch(rmm_gpu_for_tex, (float)(bj + j), (float)i2);
          w2  += rdm_this_thread2 * fjreg;
          w32 += fgjreg * rdm_this_thread2;

          // Transition density
          rdm_this_thread2 = fetch(tred_gpu_for_tex, (float)(bj + j), (float)i2);
          z2  += rdm_this_thread2 * fjreg;
          z32 += fgjreg * rdm_this_thread2;

          // Difference density
          rdm_this_thread2 = fetch(diff_gpu_for_tex, (float)(bj + j), (float)i2);
          x2  += rdm_this_thread2 * fjreg;
          x32 += fgjreg * rdm_this_thread2;
        }
      } // loop j
    } // valid tread
  } // loop bj

  // GS density
  scalar_type partial_density(0.0f);
  vec_type<scalar_type, 3> dxyz = vec_type<scalar_type, 3>(0.0f,0.0f,0.0f);

  // Transition density
  scalar_type partial_tred(0.0f);
  vec_type<scalar_type, 3> tredxyz = vec_type<scalar_type, 3>(0.0f,0.0f,0.0f);

  // Difference density
  scalar_type partial_diff(0.0f);
  vec_type<scalar_type, 3> diffxyz = vec_type<scalar_type, 3>(0.0f,0.0f,0.0f);


  if (valid_thread) {
    scalar_type Fi = function_values[(m)*point + i];
    vec_type<scalar_type, 3> Fgi;
    Fgi = vec_type<scalar_type, 3>(gradient_values[(m)*point + i]);
    // GS density
    partial_density = Fi * w;
    dxyz = Fgi * w + w3 * Fi;
    // Transition density
    partial_tred = Fi * z;
    tredxyz = Fgi * z + z3 * Fi;
    // Difference density
    partial_diff = Fi * x;
    diffxyz = Fgi * x + x3 * Fi;

    if (valid_thread2) {
      scalar_type Fi2 = function_values[(m)*point + i2];
      vec_type<scalar_type, 3> Fgi2;
      Fgi2 = vec_type<scalar_type, 3>(gradient_values[(m)*point + i2]);
      // GS density
      partial_density += Fi2 * w2;
      dxyz += Fgi2 * w2 + w32 * Fi2;
      // Transition density
      partial_tred += Fi2 * z2;
      tredxyz += Fgi2 * z2 + z32 * Fi2;
      // Difference density
      partial_diff += Fi2 * x2;
      diffxyz += Fgi2 * x2 + x32 * Fi2;
    } // valid thread 2
  } // valid thread

  __syncthreads();
// GS density
  fj_sh[position] = partial_density;
  fgj_sh[position] = dxyz;
// Transition density
  fj_sh_tred[position] = partial_tred;
  fgj_sh_tred[position] = tredxyz;
// Difference density
  fj_sh_diff[position] = partial_diff;
  fgj_sh_diff[position] = diffxyz;
  __syncthreads();

  for (int j = 2; j <= DENSITY_BLOCK_SIZE; j = j * 2) {
    int index = position + DENSITY_BLOCK_SIZE / j;
    if (position < DENSITY_BLOCK_SIZE / j) {
      // GS density
      fj_sh[position] += fj_sh[index];
      fgj_sh[position] += fgj_sh[index];

      // Transition density
      fj_sh_tred[position] += fj_sh_tred[index];
      fgj_sh_tred[position] += fgj_sh_tred[index];

      // Difference density
      fj_sh_diff[position] += fj_sh_diff[index];
      fgj_sh_diff[position] += fgj_sh_diff[index];
    }
  }

  if (threadIdx.x == 0) {
    const int myPoint = blockIdx.y * points + blockIdx.x;
    // GS density
    out_partial_density[myPoint] = fj_sh[position];
    out_dxyz[myPoint] = vec_type<scalar_type, 4>(fgj_sh[position]);
    // Transition density
    out_partial_tred[myPoint] = fj_sh_tred[position];
    out_tredxyz[myPoint] = vec_type<scalar_type, 4>(fgj_sh_tred[position]);
    // Difference density
    out_partial_diff[myPoint] = fj_sh_diff[position];
    out_diffxyz[myPoint] = vec_type<scalar_type, 4>(fgj_sh_diff[position]);
  }
}
