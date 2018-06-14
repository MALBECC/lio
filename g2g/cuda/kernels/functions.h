
// -*- mode: c -*-

/* TODO: coalescear contractions y demas */
template <class scalar_type, bool do_forces, bool do_gga>
static __device__ __host__ void compute_function(
    uint m, uint idx, vec_type<scalar_type, 3> point_position,
    uint contractions, scalar_type* factor_a_sh, scalar_type* factor_c_sh,
    uint nuc, scalar_type& t, scalar_type& tg, scalar_type& th,
    vec_type<scalar_type, 3>& v) {
  vec_type<scalar_type, 3> atom_nuc_position(
      gpu_atom_positions[nuc]);  // TODO: ver si al usar memoria compartida para
                                 // esto, pago menos precio por todos los misses
  v = point_position - atom_nuc_position;
  scalar_type dist = length2(v);

  t = 0.0f;
  if (do_forces || do_gga) tg = 0.0f;
  if (do_gga) th = 0.0f;

  for (uint contraction = 0; contraction < contractions; contraction++) {
    scalar_type t0 =
        expf(-(factor_a_sh[contraction] * dist)) * factor_c_sh[contraction];
    t += t0;
    if (do_forces || do_gga) tg += t0 * factor_a_sh[contraction];
    if (do_gga)
      th += t0 * (factor_a_sh[contraction] * factor_a_sh[contraction]);
  }
}

/**
 * gpu_compute_functions
 */

template <class scalar_type, bool do_forces, bool do_gga>
__global__ void gpu_compute_functions(vec_type<scalar_type, 4>* point_positions,
                                      uint points, uint* contractions,
                                      vec_type<scalar_type, 2>* factor_ac,
                                      uint* nuc, scalar_type* function_values,
                                      vec_type<scalar_type, 4>* gradient_values,
                                      vec_type<scalar_type, 4>* hessian_values,
                                      uint4 functions) {
  dim3 pos = index(blockDim, blockIdx, threadIdx);
  uint point = pos.x;

  /**** Load Point Information ****/
  bool valid_thread = (point < points);
  vec_type<scalar_type, 3> point_position;
  if (valid_thread) {
    vec_type<scalar_type, 4> point_position4 = point_positions[point];
    point_position = vec_type<scalar_type, 3>(point_position4);
  }

  /** Compute functions ***/
  scalar_type t, tg, th;
  vec_type<scalar_type, 3> v;

  __shared__ uint nuc_sh[FUNCTIONS_BLOCK_SIZE];
  __shared__ uint contractions_sh[FUNCTIONS_BLOCK_SIZE];
  __shared__ scalar_type factor_a_sh[FUNCTIONS_BLOCK_SIZE][MAX_CONTRACTIONS];
  __shared__ scalar_type factor_c_sh[FUNCTIONS_BLOCK_SIZE][MAX_CONTRACTIONS];

  for (uint i = 0; i < functions.w; i += FUNCTIONS_BLOCK_SIZE) {
    if (i + threadIdx.x < functions.w) {
      nuc_sh[threadIdx.x] = nuc[i + threadIdx.x];
      contractions_sh[threadIdx.x] = contractions[i + threadIdx.x];
      for (uint contraction = 0; contraction < contractions_sh[threadIdx.x];
           contraction++) {
        vec_type<scalar_type, 2> factor_ac_local =
            factor_ac[COALESCED_DIMENSION(functions.w) * contraction +
                      (i + threadIdx.x)];
        factor_a_sh[threadIdx.x][contraction] = factor_ac_local.x;
        factor_c_sh[threadIdx.x][contraction] = factor_ac_local.y;
      }
    }

    __syncthreads();

    // TODO: se podrian evitar los modulos
    if (valid_thread) {
      for (uint ii = 0; ii < FUNCTIONS_BLOCK_SIZE && (i + ii < functions.w);
           ii++) {
        compute_function<scalar_type, do_forces, do_gga>(
            functions.w, ii, point_position, contractions_sh[ii],
            factor_a_sh[ii], factor_c_sh[ii], nuc_sh[ii], t, tg, th, v);
        uint idx = COALESCED_DIMENSION(points) * (i + ii) + point;
        uint hidx1, hidx2;
        if (do_gga) {
          hidx1 = COALESCED_DIMENSION(points) * (2 * (i + ii) + 0) + point;
          hidx2 = COALESCED_DIMENSION(points) * (2 * (i + ii) + 1) + point;
        }

        vec_type<scalar_type, 3> vxxy, vyzz;
        if (do_gga) {
          vxxy = vec_type<scalar_type, 3>(v.x, v.x, v.y);
          vyzz = vec_type<scalar_type, 3>(v.y, v.z, v.z);
        }

        if (i + ii < functions.x) {
          function_values[idx] = t;
          if (do_forces || do_gga)
            gradient_values[idx] = vec_type<scalar_type, 4>(v * (-2.0f * tg));
          if (do_gga) {
            hessian_values[hidx1] = vec_type<scalar_type, 4>(
                (v * v) * 4.0f * th - 2.0f * tg);  // Fxx, Fxy, Fxz
            hessian_values[hidx2] = vec_type<scalar_type, 4>(
                vxxy * vyzz * 4.0f * th);  // Fxy, Fxz, Fyz
          }
        } else if (i + ii < (functions.x + functions.y * 3)) {
          uint p_idx = ((i + ii) - functions.x) % 3;
          switch (p_idx) {
            case 0:
              function_values[idx] = v.x * t;
              if (do_forces || do_gga)
                gradient_values[idx] = vec_type<scalar_type, 4>(
                    vec_type<scalar_type, 3>(t, 0.0f, 0.0f) -
                    v * 2.0f * tg * v.x);
              if (do_gga) {
                hessian_values[hidx1] = vec_type<scalar_type, 4>(
                    (v * v) * 4.0f * th * v.x -
                    vec_type<scalar_type, 3>(6.0f, 2.0f, 2.0f) * tg * v.x);
                hessian_values[hidx2] = vec_type<scalar_type, 4>(
                    (vxxy * vyzz) * 4.0f * th * v.x -
                    vec_type<scalar_type, 3>(v.y, v.z, 0.0f) * 2.0f * tg);
              }
              break;
            case 1:
              function_values[idx] = v.y * t;
              if (do_forces || do_gga)
                gradient_values[idx] = vec_type<scalar_type, 4>(
                    vec_type<scalar_type, 3>(0.0f, t, 0.0f) -
                    v * 2.0f * tg * v.y);
              if (do_gga) {
                hessian_values[hidx1] = vec_type<scalar_type, 4>(
                    (v * v) * 4.0f * th * v.y -
                    vec_type<scalar_type, 3>(2.0f, 6.0f, 2.0f) * tg * v.y);
                hessian_values[hidx2] = vec_type<scalar_type, 4>(
                    (vxxy * vyzz) * 4.0f * th * v.y -
                    vec_type<scalar_type, 3>(v.x, 0.0f, v.z) * 2.0f * tg);
              }
              break;
            case 2:
              function_values[idx] = v.z * t;
              if (do_forces || do_gga)
                gradient_values[idx] = vec_type<scalar_type, 4>(
                    vec_type<scalar_type, 3>(0.0f, 0.0f, t) -
                    v * 2.0f * tg * v.z);
              if (do_gga) {
                hessian_values[hidx1] = vec_type<scalar_type, 4>(
                    (v * v) * 4.0f * th * v.z -
                    vec_type<scalar_type, 3>(2.0f, 2.0f, 6.0f) * tg * v.z);
                hessian_values[hidx2] = vec_type<scalar_type, 4>(
                    (vxxy * vyzz) * 4.0f * th * v.z -
                    vec_type<scalar_type, 3>(0.0f, v.x, v.y) * 2.0f * tg);
              }
              break;
          }
        } else {
          uint d_idx = ((i + ii) - functions.x - functions.y * 3) % 6;
          switch (d_idx) {
            case 0:
              function_values[idx] = t * v.x * v.x * gpu_normalization_factor;
              if (do_forces || do_gga)
                gradient_values[idx] = vec_type<scalar_type, 4>(
                    (vec_type<scalar_type, 3>(2.0f * v.x, 0.0f, 0.0f) * t -
                     v * 2.0f * tg * v.x * v.x) *
                    gpu_normalization_factor);
              if (do_gga) {
                hessian_values[hidx1] = vec_type<scalar_type, 4>(
                    ((v * v) * 4.0f * th * (v.x * v.x) -
                     vec_type<scalar_type, 3>(10.0f, 2.0f, 2.0f) * tg *
                         (v.x * v.x) +
                     vec_type<scalar_type, 3>(2.0f * t, 0.0f, 0.0f)) *
                    gpu_normalization_factor);
                hessian_values[hidx2] = vec_type<scalar_type, 4>(
                    ((vxxy * vyzz) * 4.0f * th * (v.x * v.x) -
                     vec_type<scalar_type, 3>(4.0f, 4.0f, 0.0f) *
                         (vxxy * vyzz) * tg) *
                    gpu_normalization_factor);
              }
              break;
            case 1:
              function_values[idx] = t * v.y * v.x;
              if (do_forces || do_gga)
                gradient_values[idx] = vec_type<scalar_type, 4>(
                    vec_type<scalar_type, 3>(v.y, v.x, 0.0f) * t -
                    v * 2.0f * tg * v.y * v.x);
              if (do_gga) {
                hessian_values[hidx1] = vec_type<scalar_type, 4>(
                    ((v * v) * 4.0f * th * (v.x * v.y) -
                     vec_type<scalar_type, 3>(6.0f, 6.0f, 2.0f) * tg *
                         (v.x * v.y)));
                hessian_values[hidx2] = vec_type<scalar_type, 4>(
                    ((vxxy * vyzz) * 4.0f * th * (v.x * v.y) -
                     vec_type<scalar_type, 3>(2.0f * (v.x * v.x + v.y * v.y),
                                              2.0f * v.y * v.z,
                                              2.0f * v.x * v.z) *
                         tg +
                     vec_type<scalar_type, 3>(t, 0.0f, 0.0f)));
              }
              break;
            case 2:
              function_values[idx] = t * v.y * v.y * gpu_normalization_factor;
              if (do_forces || do_gga)
                gradient_values[idx] = vec_type<scalar_type, 4>(
                    (vec_type<scalar_type, 3>(0.0f, 2.0f * v.y, 0.0f) * t -
                     v * 2.0f * tg * v.y * v.y) *
                    gpu_normalization_factor);
              if (do_gga) {
                hessian_values[hidx1] = vec_type<scalar_type, 4>(
                    ((v * v) * 4.0f * th * (v.y * v.y) -
                     vec_type<scalar_type, 3>(2.0f, 10.0f, 2.0f) * tg *
                         (v.y * v.y) +
                     vec_type<scalar_type, 3>(0.0f, 2.0f * t, 0.0f)) *
                    gpu_normalization_factor);
                hessian_values[hidx2] = vec_type<scalar_type, 4>(
                    ((vxxy * vyzz) * 4.0f * th * (v.y * v.y) -
                     vec_type<scalar_type, 3>(4.0f, 0.0f, 4.0f) *
                         (vxxy * vyzz) * tg) *
                    gpu_normalization_factor);
              }
              break;
            case 3:
              function_values[idx] = t * v.z * v.x;
              if (do_forces || do_gga)
                gradient_values[idx] = vec_type<scalar_type, 4>(
                    vec_type<scalar_type, 3>(v.z, 0.0f, v.x) * t -
                    v * 2.0f * tg * v.z * v.x);
              if (do_gga) {
                hessian_values[hidx1] = vec_type<scalar_type, 4>(
                    ((v * v) * 4.0f * th * (v.x * v.z) -
                     vec_type<scalar_type, 3>(6.0f, 2.0f, 6.0f) * tg *
                         (v.x * v.z)));
                hessian_values[hidx2] = vec_type<scalar_type, 4>(
                    ((vxxy * vyzz) * 4.0f * th * (v.x * v.z) -
                     vec_type<scalar_type, 3>(2.0f * v.y * v.z,
                                              2.0f * (v.x * v.x + v.z * v.z),
                                              2.0f * v.x * v.y) *
                         tg +
                     vec_type<scalar_type, 3>(0.0f, t, 0.0f)));
              }
              break;
            case 4:
              function_values[idx] = t * v.z * v.y;
              if (do_forces || do_gga)
                gradient_values[idx] = vec_type<scalar_type, 4>(
                    vec_type<scalar_type, 3>(0.0f, v.z, v.y) * t -
                    v * 2.0f * tg * v.z * v.y);
              if (do_gga) {
                hessian_values[hidx1] = vec_type<scalar_type, 4>(
                    ((v * v) * 4.0f * th * (v.y * v.z) -
                     vec_type<scalar_type, 3>(2.0f, 6.0f, 6.0f) * tg *
                         (v.y * v.z)));
                hessian_values[hidx2] = vec_type<scalar_type, 4>((
                    (vxxy * vyzz) * 4.0f * th * (v.y * v.z) -
                    vec_type<scalar_type, 3>(2.0f * v.x * v.z, 2.0f * v.x * v.y,
                                             2.0f * (v.y * v.y + v.z * v.z)) *
                        tg +
                    vec_type<scalar_type, 3>(0.0f, 0.0f, t)));
              }
              break;
            case 5:
              function_values[idx] = t * v.z * v.z * gpu_normalization_factor;
              if (do_forces || do_gga)
                gradient_values[idx] = vec_type<scalar_type, 4>(
                    (vec_type<scalar_type, 3>(0.0f, 0.0f, 2.0f * v.z) * t -
                     v * 2.0f * tg * v.z * v.z) *
                    gpu_normalization_factor);
              if (do_gga) {
                hessian_values[hidx1] = vec_type<scalar_type, 4>(
                    ((v * v) * 4.0f * th * (v.z * v.z) -
                     vec_type<scalar_type, 3>(2.0f, 2.0f, 10.0f) * tg *
                         (v.z * v.z) +
                     vec_type<scalar_type, 3>(0.0f, 0.0f, 2.0f * t)) *
                    gpu_normalization_factor);
                hessian_values[hidx2] = vec_type<scalar_type, 4>(
                    ((vxxy * vyzz) * 4.0f * th * (v.z * v.z) -
                     (vxxy * vyzz) *
                         vec_type<scalar_type, 3>(0.0f, 4.0f, 4.0f) * tg) *
                    gpu_normalization_factor);
              }
              break;
          }
        }
      }
    }

    __syncthreads();
  }
}
