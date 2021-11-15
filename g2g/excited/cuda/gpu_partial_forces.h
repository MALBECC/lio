template <class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_partial_forces(
    const T* fv, const vec_type<T, 4>* gv, const vec_type<T, 4>* hv,
    const T* mat_dens, const T* mat_tred, const T* mat_diff, const T* dens,
    const T* tred, const T* diff, const G2G::vec_type<T, 4>* dxyz,
    const G2G::vec_type<T, 4>* tredxyz, const G2G::vec_type<T, 4>* diffxyz,
    G2G::vec_type<T, 4>* force, int basis_tot, uint npoints,
    const T* const weights, int dim_block) {
  int ii = blockIdx.x;
  int block_ii = threadIdx.x;
  int start_point = block_ii * POINTS_BLOCK_SIZE;
  int end_point = start_point + POINTS_BLOCK_SIZE;
  int dim = basis_tot;
  int element = ii * dim_block + block_ii;

  double fvj;
  double factor, temp[4];
  double DJII, PJII, VJII;
  double wp;
  double grdx, grdy, grdz;
  vec_type<T, 4> h1, h2;  // h1 = hxx, hyy, hzz // h2 = hxy, hxz, hyz //
  vec_type<T, 4> gfii, gfj;

  for (int point = start_point; point < end_point; point++) {
    if (point < npoints) {
      h1 = vec_type<T, 4>(hv[basis_tot * 2 * point + (2 * ii + 0)]);
      h2 = vec_type<T, 4>(hv[basis_tot * 2 * point + (2 * ii + 1)]);
      gfii = vec_type<T, 4>(gv[basis_tot * point + ii]);
      wp = weights[point];
      grdx = grdy = grdz = 0.0f;

      for (int j = 0; j < basis_tot; j++) {
        fvj = fv[basis_tot * point + j];
        gfj = vec_type<T, 4>(gv[basis_tot * point + j]);
        factor = (ii == j ? 2.0f : 1.0f);
        DJII = mat_dens[j + ii * dim] * factor * 0.5f;
        PJII = mat_diff[j + ii * dim] * factor;
        VJII = mat_tred[j + ii * dim] * factor;
        temp[0] = 2.0f * (dens[point] * DJII + diff[point] * PJII +
                          tred[point] * VJII);
        temp[1] = 2.0f * (dxyz[point].x * DJII + diffxyz[point].x * PJII +
                          tredxyz[point].x * VJII);
        temp[2] = 2.0f * (dxyz[point].y * DJII + diffxyz[point].y * PJII +
                          tredxyz[point].y * VJII);
        temp[3] = 2.0f * (dxyz[point].z * DJII + diffxyz[point].z * PJII +
                          tredxyz[point].z * VJII);

        grdx += temp[0] * gfii.x * fvj;
        grdx += temp[1] * (h1.x * fvj + gfii.x * gfj.x);
        grdx += temp[2] * (h2.x * fvj + gfii.x * gfj.y);
        grdx += temp[3] * (h2.y * fvj + gfii.x * gfj.z);

        grdy += temp[0] * gfii.y * fvj;
        grdy += temp[1] * (h2.x * fvj + gfii.y * gfj.x);
        grdy += temp[2] * (h1.y * fvj + gfii.y * gfj.y);
        grdy += temp[3] * (h2.z * fvj + gfii.y * gfj.z);

        grdz += temp[0] * gfii.z * fvj;
        grdz += temp[1] * (h2.y * fvj + gfii.z * gfj.x);
        grdz += temp[2] * (h2.z * fvj + gfii.z * gfj.y);
        grdz += temp[3] * (h1.z * fvj + gfii.z * gfj.z);
      }  // end j
      force[element] +=
          G2G::vec_type<T, WIDTH>(grdx * wp, grdy * wp, grdz * wp, 0.0f);
    }  // end if point
  }    // end point
}

template <class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_accum_forces(const G2G::vec_type<T, 4>* force_in,
                                 G2G::vec_type<T, 4>* force_out, int dim_block,
                                 int M) {
  int ii = blockIdx.x;
  for (int jj = 0; jj < dim_block; jj++)
    force_out[ii] += force_in[ii * dim_block + jj];
}
