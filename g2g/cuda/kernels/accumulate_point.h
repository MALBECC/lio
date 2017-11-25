#define WIDTH 4

// OPEN SHELL CASE
template <class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_accumulate_point_open(
    scalar_type* const energy, scalar_type* const energy_i,
    scalar_type* const energy_c, scalar_type* const energy_c1,
    scalar_type* const energy_c2, scalar_type* const factor_a,
    scalar_type* const factor_b, const scalar_type* const point_weights,
    uint points, int block_height, scalar_type* partial_density_a,
    vec_type<scalar_type, WIDTH>* dxyz_a, vec_type<scalar_type, WIDTH>* dd1_a,
    vec_type<scalar_type, WIDTH>* dd2_a, scalar_type* partial_density_b,
    vec_type<scalar_type, WIDTH>* dxyz_b, vec_type<scalar_type, WIDTH>* dd1_b,
    vec_type<scalar_type, WIDTH>* dd2_b) {
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  scalar_type point_weight = 0.0f;
  scalar_type exc_corr, v_a, v_b, exc, corr, corr1, corr2;
  exc_corr = v_a = v_b = exc = corr = corr1 = corr2 = 0.0f;

  scalar_type _partial_density_a(0.0f), _partial_density_b(0.0f);
  vec_type<scalar_type, WIDTH> _dxyz_a, _dd1_a, _dd2_a, _dxyz_b, _dd1_b, _dd2_b;
  _dxyz_a = _dd1_a = _dd2_a = _dxyz_b = _dd1_b = _dd2_b =
      vec_type<scalar_type, WIDTH>(0.0f, 0.0f, 0.0f, 0.0f);

  bool valid_thread = (point < points);
  if (valid_thread) point_weight = point_weights[point];

  if (valid_thread) {
    for (int j = 0; j < block_height; j++) {
      const int this_row = j * points + point;
      _partial_density_a += partial_density_a[this_row];
      _partial_density_b += partial_density_b[this_row];
      _dxyz_a += dxyz_a[this_row];
      _dxyz_b += dxyz_b[this_row];
      _dd1_a += dd1_a[this_row];
      _dd1_b += dd1_b[this_row];
      _dd2_a += dd2_a[this_row];
      _dd2_b += dd2_b[this_row];
    }
  }
  calc_ggaOS<scalar_type, WIDTH>(
      _partial_density_a, _partial_density_b, _dxyz_a, _dxyz_b, _dd1_a, _dd1_b,
      _dd2_a, _dd2_b, exc_corr, exc, corr, corr1, corr2, v_a, v_b, 9);

  if (compute_energy && valid_thread) {
    energy[point] =
        ((_partial_density_a + _partial_density_b) * point_weight) * exc_corr;
    energy_i[point] =
        ((_partial_density_a + _partial_density_b) * point_weight) * exc;
    energy_c[point] =
        ((_partial_density_a + _partial_density_b) * point_weight) * corr;
    energy_c1[point] =
        ((_partial_density_a + _partial_density_b) * point_weight) * corr1;
    energy_c2[point] =
        ((_partial_density_a + _partial_density_b) * point_weight) * corr2;
  }

  if (compute_factor && valid_thread) {
    factor_a[point] = point_weight * v_a;
    factor_b[point] = point_weight * v_b;
  }
}

template <class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_accumulate_point(
    scalar_type* const energy, scalar_type* const factor,
    const scalar_type* const point_weights, uint points, int block_height,
    scalar_type* partial_density, vec_type<scalar_type, WIDTH>* dxyz,
    vec_type<scalar_type, WIDTH>* dd1, vec_type<scalar_type, WIDTH>* dd2) {
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  scalar_type point_weight = 0.0f;
  scalar_type y2a, exc_corr, exc_c, exc_x;

  scalar_type _partial_density(0.0f);
  vec_type<scalar_type, WIDTH> _dxyz, _dd1, _dd2;

  _dxyz = _dd1 = _dd2 = vec_type<scalar_type, WIDTH>(0.0f, 0.0f, 0.0f, 0.0f);

  bool valid_thread = (point < points);
  if (valid_thread) point_weight = point_weights[point];

  if (valid_thread) {
    for (int j = 0; j < block_height; j++) {
      const int this_row = j * points + point;

      _partial_density += partial_density[this_row];
      _dxyz += dxyz[this_row];
      _dd1 += dd1[this_row];
      _dd2 += dd2[this_row];
    }
  }

  calc_ggaCS_in<scalar_type, 4>(_partial_density, _dxyz, _dd1, _dd2, exc_x,
                                exc_c, y2a, 9);
  exc_corr = exc_x + exc_c;

  if (compute_energy && valid_thread)
    energy[point] = (_partial_density * point_weight) * exc_corr;

  if (compute_factor && valid_thread) factor[point] = point_weight * y2a;
}
