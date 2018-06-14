#if FULL_DOUBLE
/*
static __inline__ __device__ double fetch_double(texture<int2, 2> t, float x,
float y)
{
   int2 v = tex2D(t,x,y);
   return __hiloint2double(v.y, v.x);
}*/
#define fetch(t, x, y) fetch_double(t, x, y)
#else
#define fetch(t, x, y) tex2D(t, x, y)
#endif

template <class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_compute_density_opened(
    const scalar_type* point_weights, uint points,
    const scalar_type* function_values,
    const vec_type<scalar_type, 4>* gradient_values,
    const vec_type<scalar_type, 4>* hessian_values, uint m,
    scalar_type* out_partial_density_a, vec_type<scalar_type, 4>* out_dxyz_a,
    vec_type<scalar_type, 4>* out_dd1_a, vec_type<scalar_type, 4>* out_dd2_a,
    scalar_type* out_partial_density_b, vec_type<scalar_type, 4>* out_dxyz_b,
    vec_type<scalar_type, 4>* out_dd1_b, vec_type<scalar_type, 4>* out_dd2_b) {
  uint point = blockIdx.x;
  uint i = threadIdx.x + blockIdx.y * 2 * DENSITY_BLOCK_SIZE;
  uint i2 = i + DENSITY_BLOCK_SIZE;
  uint min_i = blockIdx.y * 2 * DENSITY_BLOCK_SIZE + DENSITY_BLOCK_SIZE;

  scalar_type partial_density_a(0.0f);
  scalar_type partial_density_b(0.0f);

  vec_type<scalar_type, 3> dxyz_a, dd1_a, dd2_a;
  vec_type<scalar_type, 3> dxyz_b, dd1_b, dd2_b;
  dxyz_a = dd1_a = dd2_a = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
  dxyz_b = dd1_b = dd2_b = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

  bool valid_thread = (i < m);
  bool valid_thread2 = (i2 < m);

  scalar_type w_a = 0.0f;
  scalar_type w_b = 0.0f;
  scalar_type w2_a = 0.0f;
  scalar_type w2_b = 0.0f;
  vec_type<scalar_type, 3> w3_a, ww1_a, ww2_a;
  vec_type<scalar_type, 3> w3_b, ww1_b, ww2_b;
  vec_type<scalar_type, 3> w32_a, ww12_a, ww22_a;
  vec_type<scalar_type, 3> w32_b, ww12_b, ww22_b;

  if (!lda) {
    w3_a = ww1_a = ww2_a = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    w3_b = ww1_b = ww2_b = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    w32_a = ww12_a = ww22_a = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    w32_b = ww12_b = ww22_b = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
  }

  scalar_type Fi, Fi2;
  vec_type<scalar_type, 3> Fgi, Fhi1, Fhi2, Fgi2, Fhi12, Fhi22;

  int position = threadIdx.x;

  __shared__ scalar_type fj_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fgj_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fh1j_sh[DENSITY_BLOCK_SIZE];
  __shared__ vec_type<scalar_type, 3> fh2j_sh[DENSITY_BLOCK_SIZE];

  uint min_i2 = min_i;

  // Si nos vamos a pasar del bloque con el segundo puntero, hacemos que haga la
  // misma cuenta
  if (min_i > m) {
    min_i2 = min_i - DENSITY_BLOCK_SIZE;
  }

  if (valid_thread2) {
    Fi2 = function_values[(m)*point + i2];
    if (!lda) {
      Fgi2 = vec_type<scalar_type, 3>(gradient_values[(m)*point + i2]);
      Fhi12 = vec_type<scalar_type, 3>(
          hessian_values[(m)*2 * point + (2 * i2 + 0)]);
      Fhi22 = vec_type<scalar_type, 3>(
          hessian_values[(m)*2 * point + (2 * i2 + 1)]);
    }
  }
  if (valid_thread) {
    Fi = function_values[(m)*point + i];
    if (!lda) {
      Fgi = vec_type<scalar_type, 3>(gradient_values[(m)*point + i]);
      Fhi1 =
          vec_type<scalar_type, 3>(hessian_values[(m)*2 * point + (2 * i + 0)]);
      Fhi2 =
          vec_type<scalar_type, 3>(hessian_values[(m)*2 * point + (2 * i + 1)]);
    }
  }

  for (int bj = 0; bj <= min_i2; bj += DENSITY_BLOCK_SIZE) {
    // Density deberia ser GET_DENSITY_BLOCK_SIZE

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
    scalar_type fjreg = 0.0f;
    vec_type<scalar_type, 3> fgjreg =
        vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    vec_type<scalar_type, 3> fh1jreg =
        vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    vec_type<scalar_type, 3> fh2jreg =
        vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);

    if (valid_thread) {
      for (int j = 0; j < DENSITY_BLOCK_SIZE; j++) {
        fjreg = fj_sh[j];
        if (!lda) {
          fgjreg = fgj_sh[j];
          fh1jreg = fh1j_sh[j];
          fh2jreg = fh2j_sh[j];
        }

        if ((bj + j) <= i) {
          // Fetch is a  macro for tex2D
          scalar_type rdm_this_thread_a =
              fetch(rmm_input_gpu_tex, (float)(bj + j), (float)i);
          scalar_type rdm_this_thread_b =
              fetch(rmm_input_gpu_tex2, (float)(bj + j), (float)i);

          w_a += rdm_this_thread_a * fjreg;
          w_b += rdm_this_thread_b * fjreg;

          if (!lda) {
            w3_a += fgjreg * rdm_this_thread_a;
            ww1_a += fh1jreg * rdm_this_thread_a;
            ww2_a += fh2jreg * rdm_this_thread_a;

            w3_b += fgjreg * rdm_this_thread_b;
            ww1_b += fh1jreg * rdm_this_thread_b;
            ww2_b += fh2jreg * rdm_this_thread_b;
          }
        }

        if (valid_thread2 && ((bj + j) <= i2)) {
          scalar_type rdm_this_thread2_a =
              fetch(rmm_input_gpu_tex, (float)(bj + j), (float)i2);
          scalar_type rdm_this_thread2_b =
              fetch(rmm_input_gpu_tex2, (float)(bj + j), (float)i2);
          w2_a += rdm_this_thread2_a * fjreg;
          w2_b += rdm_this_thread2_b * fjreg;
          if (!lda) {
            w32_a += fgjreg * rdm_this_thread2_a;
            ww12_a += fh1jreg * rdm_this_thread2_a;
            ww22_a += fh2jreg * rdm_this_thread2_a;

            w32_b += fgjreg * rdm_this_thread2_b;
            ww12_b += fh1jreg * rdm_this_thread2_b;
            ww22_b += fh2jreg * rdm_this_thread2_b;
          }
        }
      }
    }
  }
  if (valid_thread) {
    partial_density_a = Fi * w_a;
    partial_density_b = Fi * w_b;
    if (!lda) {
      dxyz_a = Fgi * w_a + w3_a * Fi;
      dd1_a = Fgi * w3_a * 2.0f + Fhi1 * w_a + ww1_a * Fi;

      dxyz_b = Fgi * w_b + w3_b * Fi;
      dd1_b = Fgi * w3_b * 2.0f + Fhi1 * w_b + ww1_b * Fi;

      vec_type<scalar_type, 3> FgXXY(Fgi.x, Fgi.x, Fgi.y);
      vec_type<scalar_type, 3> w3YZZ_a(w3_a.y, w3_a.z, w3_a.z);
      vec_type<scalar_type, 3> w3YZZ_b(w3_b.y, w3_b.z, w3_b.z);
      vec_type<scalar_type, 3> FgiYZZ(Fgi.y, Fgi.z, Fgi.z);
      vec_type<scalar_type, 3> w3XXY_a(w3_a.x, w3_a.x, w3_a.y);
      vec_type<scalar_type, 3> w3XXY_b(w3_b.x, w3_b.x, w3_b.y);

      dd2_a = FgXXY * w3YZZ_a + FgiYZZ * w3XXY_a + Fhi2 * w_a + ww2_a * Fi;
      dd2_b = FgXXY * w3YZZ_b + FgiYZZ * w3XXY_b + Fhi2 * w_b + ww2_b * Fi;
    }
  }
  if (valid_thread2) {
    partial_density_a += Fi2 * w2_a;
    partial_density_b += Fi2 * w2_b;
    if (!lda) {
      dxyz_a += Fgi2 * w2_a + w32_a * Fi2;
      dd1_a += Fgi2 * w32_a * 2.0f + Fhi12 * w2_a + ww12_a * Fi2;
      dxyz_b += Fgi2 * w2_b + w32_b * Fi2;
      dd1_b += Fgi2 * w32_b * 2.0f + Fhi12 * w2_b + ww12_b * Fi2;

      vec_type<scalar_type, 3> FgXXY(Fgi2.x, Fgi2.x, Fgi2.y);
      vec_type<scalar_type, 3> w3YZZ_a(w32_a.y, w32_a.z, w32_a.z);
      vec_type<scalar_type, 3> w3YZZ_b(w32_b.y, w32_b.z, w32_b.z);
      vec_type<scalar_type, 3> FgiYZZ(Fgi2.y, Fgi2.z, Fgi2.z);
      vec_type<scalar_type, 3> w3XXY_a(w32_a.x, w32_a.x, w32_a.y);
      vec_type<scalar_type, 3> w3XXY_b(w32_b.x, w32_b.x, w32_b.y);

      dd2_a += FgXXY * w3YZZ_a + FgiYZZ * w3XXY_a + Fhi22 * w2_a + ww22_a * Fi2;
      dd2_b += FgXXY * w3YZZ_b + FgiYZZ * w3XXY_b + Fhi22 * w2_b + ww22_b * Fi2;
    }
  }

  __syncthreads();
  // Reusing per-block shared memory in order to perform per-block accumulation.
  // Alpha density
  if (valid_thread) {
    fj_sh[position] = partial_density_a;
    fgj_sh[position] = dxyz_a;
    fh1j_sh[position] = dd1_a;
    fh2j_sh[position] = dd2_a;
  } else {
    fj_sh[position] = scalar_type(0.0f);
    fgj_sh[position] = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    fh1j_sh[position] = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    fh2j_sh[position] = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
  }
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
    out_partial_density_a[myPoint] = fj_sh[position];
    out_dxyz_a[myPoint] = vec_type<scalar_type, 4>(fgj_sh[position]);
    out_dd1_a[myPoint] = vec_type<scalar_type, 4>(fh1j_sh[position]);
    out_dd2_a[myPoint] = vec_type<scalar_type, 4>(fh2j_sh[position]);
  }

  __syncthreads();
  // Beta density.
  if (valid_thread) {
    fj_sh[position] = partial_density_b;
    fgj_sh[position] = dxyz_b;
    fh1j_sh[position] = dd1_b;
    fh2j_sh[position] = dd2_b;
  } else {
    fj_sh[position] = scalar_type(0.0f);
    fgj_sh[position] = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    fh1j_sh[position] = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
    fh2j_sh[position] = vec_type<scalar_type, 3>(0.0f, 0.0f, 0.0f);
  }
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
    out_partial_density_b[myPoint] = fj_sh[position];
    out_dxyz_b[myPoint] = vec_type<scalar_type, 4>(fgj_sh[position]);
    out_dd1_b[myPoint] = vec_type<scalar_type, 4>(fh1j_sh[position]);
    out_dd2_b[myPoint] = vec_type<scalar_type, 4>(fh2j_sh[position]);
  }
}
