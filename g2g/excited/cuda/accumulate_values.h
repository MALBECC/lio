template <class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void accumulate_values(uint points, int block_height,
                                  T* partial_density_in,
                                  G2G::vec_type<T, WIDTH>* dxyz_in,
                                  T* accumulated_density,
                                  G2G::vec_type<T, WIDTH>* dxyz_accum) {
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;
  // uint point = blockIdx.x * blockDim.x + threadIdx.x;

  T _partial_density(0.0f);
  G2G::vec_type<T, WIDTH> _dxyz;
  _dxyz = G2G::vec_type<T, WIDTH>(0.0f, 0.0f, 0.0f, 0.0f);

  bool valid_thread = (point < points);

  if (valid_thread) {
    for (int j = 0; j < block_height; j++) {
      const int this_row = j * points + point;
      _partial_density += partial_density_in[this_row];
      _dxyz += dxyz_in[this_row];
    }  // END J LOOP

    accumulated_density[point] = _partial_density;
    dxyz_accum[point] = _dxyz;
  }  // END VALID THREAD
}

template <class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void accumulate_values(uint points, int block_height,
                                  T* partial_tred_in, T* partial_diff_in,
                                  G2G::vec_type<T, WIDTH>* tredxyz_in,
                                  G2G::vec_type<T, WIDTH>* diffxyz_in,
                                  T* accumulated_tred, T* accumulated_diff,
                                  G2G::vec_type<T, WIDTH>* tredxyz_accum,
                                  G2G::vec_type<T, WIDTH>* diffxyz_accum) {
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;
  // uint point = blockIdx.x * blockDim.x + threadIdx.x;

  T _partial_tred(0.0f);
  T _partial_diff(0.0f);
  G2G::vec_type<T, WIDTH> _tredxyz;
  G2G::vec_type<T, WIDTH> _diffxyz;
  _tredxyz = _diffxyz = G2G::vec_type<T, WIDTH>(0.0f, 0.0f, 0.0f, 0.0f);

  bool valid_thread = (point < points);

  if (valid_thread) {
    for (int j = 0; j < block_height; j++) {
      const int this_row = j * points + point;
      // TRANSITION density
      _partial_tred += partial_tred_in[this_row];
      _tredxyz += tredxyz_in[this_row];
      // Difference density
      _partial_diff += partial_diff_in[this_row];
      _diffxyz += diffxyz_in[this_row];
    }  // END J LOOP

    // Accumulate TRANSITION density
    accumulated_tred[point] = _partial_tred;
    tredxyz_accum[point] = _tredxyz;
    // Accumulate DIFFERENCE density
    accumulated_diff[point] = _partial_diff;
    diffxyz_accum[point] = _diffxyz;
  }  // END VALID THREAD
}
