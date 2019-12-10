template <class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void ES_compute_partial(uint points,
                                const scalar_type* function_values, uint m,
                                const vec_type<scalar_type, 4>* gradient_values,
                                scalar_type* out_partial_tred, vec_type<scalar_type, 4>* out_tredxyz)
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

  int position = threadIdx.x;

  __shared__ scalar_type fj_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh[DENSITY_BLOCK_SIZE];

  __shared__ scalar_type fj_sh_tred[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh_tred[DENSITY_BLOCK_SIZE];

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
          rdm_this_thread = fetch(tred_gpu_tex, (float)(bj + j), (float)i);
          z  += rdm_this_thread * fjreg;
          z3 += fgjreg * rdm_this_thread;
        }

        if (valid_thread2 && ((bj + j) <= i2)) {
          scalar_type rdm_this_thread2;

          // Transition density
          rdm_this_thread2 = fetch(tred_gpu_tex, (float)(bj + j), (float)i2);
          z2  += rdm_this_thread2 * fjreg;
          z32 += fgjreg * rdm_this_thread2;
        }
      } // loop j
    } // valid tread
  } // loop bj

  // Transition density
  scalar_type partial_tred(0.0f);
  vec_type<scalar_type, 3> tredxyz = vec_type<scalar_type, 3>(0.0f,0.0f,0.0f);

  if (valid_thread) {
    scalar_type Fi = function_values[(m)*point + i];
    vec_type<scalar_type, 3> Fgi;
    Fgi = vec_type<scalar_type, 3>(gradient_values[(m)*point + i]);
    // Transition density
    partial_tred = Fi * z;
    tredxyz = Fgi * z + z3 * Fi;

    if (valid_thread2) {
      scalar_type Fi2 = function_values[(m)*point + i2];
      vec_type<scalar_type, 3> Fgi2, Fhi12, Fhi22;
      Fgi2 = vec_type<scalar_type, 3>(gradient_values[(m)*point + i2]);
      // Transition density
      partial_tred += Fi2 * z2;
      tredxyz += Fgi2 * z2 + z32 * Fi2;
    } // valid thread 2
  } // valid thread

  __syncthreads();
// Transition density
  fj_sh_tred[position] = partial_tred;
  fgj_sh_tred[position] = tredxyz;
  __syncthreads();

  for (int j = 2; j <= DENSITY_BLOCK_SIZE; j = j * 2) {
    int index = position + DENSITY_BLOCK_SIZE / j;
    if (position < DENSITY_BLOCK_SIZE / j) {

      // Transition density
      fj_sh_tred[position] += fj_sh_tred[index];
      fgj_sh_tred[position] += fgj_sh_tred[index];
    }
  }

  if (threadIdx.x == 0) {
    const int myPoint = blockIdx.y * points + blockIdx.x;
    // Transition density
    out_partial_tred[myPoint] = fj_sh_tred[position];
    out_tredxyz[myPoint] = vec_type<scalar_type, 4>(fgj_sh_tred[position]);
  }
}
