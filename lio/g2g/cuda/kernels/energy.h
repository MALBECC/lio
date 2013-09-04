template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_compute_density(scalar_type* const energy, scalar_type* const factor, const scalar_type* const point_weights,
                                    uint points, const scalar_type* rdm, const scalar_type* function_values, const vec_type<scalar_type,4>* gradient_values,
                                    const vec_type<scalar_type,4>* hessian_values, uint m)
{
  uint point = index_x(blockDim, blockIdx, threadIdx);

  scalar_type partial_density = 0.0f;
  vec_type<scalar_type,4> dxyz, dd1, dd2;
  if (!lda) { dxyz = dd1 = dd2 = vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f); }

  bool valid_thread = (point < points);
  scalar_type point_weight;
  if (valid_thread) point_weight = point_weights[point];

  __shared__ scalar_type rdm_sh[DENSITY_BATCH_SIZE];

  /***** compute density ******/
  for (uint i = 0; i < m; i++) {
    scalar_type w = 0.0f;
    vec_type<scalar_type,4> w3, ww1, ww2;
    if (!lda) { w3 = ww1 = ww2 = vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f); }

    scalar_type Fi;
    vec_type<scalar_type,4> Fgi, Fhi1, Fhi2;

    if (valid_thread) {
      Fi = function_values[COALESCED_DIMENSION(points) * i + point];
      if (!lda) {
        Fgi = gradient_values[COALESCED_DIMENSION(points) * i + point];
        Fhi1 = hessian_values[COALESCED_DIMENSION(points) * (2 * i + 0) + point];
        Fhi2 = hessian_values[COALESCED_DIMENSION(points) * (2 * i + 1) + point];
      }
    }

    for (uint bj = 0; bj <= i; bj += DENSITY_BATCH_SIZE) {
      __syncthreads();
      if (threadIdx.x < DENSITY_BATCH_SIZE) {
        if (bj + threadIdx.x <= i) rdm_sh[threadIdx.x] = rdm[COALESCED_DIMENSION(m) * i + (bj + threadIdx.x)]; // TODO: uncoalesced. invertir triangulo?
        else rdm_sh[threadIdx.x] = 0.0f;
      }
      __syncthreads();

      if (valid_thread) {
        for (uint j = 0; j < DENSITY_BATCH_SIZE && (bj + j) <= i; j++) {
          float Fj = function_values[COALESCED_DIMENSION(points) * (bj + j) + point];
          w += rdm_sh[j] * Fj;

          if (!lda) {
            vec_type<scalar_type,4> Fgj = gradient_values[COALESCED_DIMENSION(points) * (bj + j) + point];
            w3 += Fgj * rdm_sh[j];

            vec_type<scalar_type,4> Fhj1 = hessian_values[COALESCED_DIMENSION(points) * (2 * (bj + j) + 0) + point];
            vec_type<scalar_type,4> Fhj2 = hessian_values[COALESCED_DIMENSION(points) * (2 * (bj + j) + 1) + point];
            ww1 += Fhj1 * rdm_sh[j];
            ww2 += Fhj2 * rdm_sh[j];
          }
        }
      }
    }

    partial_density += Fi * w;
    //TODO: Insertar aca funcion que convierte <,4> a <,3>
    if (!lda) {
      dxyz += Fgi * w + w3 * Fi;
      dd1 += Fgi * w3 * 2.0f + Fhi1 * w + ww1 * Fi;

      vec_type<scalar_type,4> FgXXY(Fgi.x, Fgi.x, Fgi.y, 0.0f);
      vec_type<scalar_type,4> w3YZZ(w3.y, w3.z, w3.z, 0.0f);
      vec_type<scalar_type,4> FgiYZZ(Fgi.y, Fgi.z, Fgi.z, 0.0f);
      vec_type<scalar_type,4> w3XXY(w3.x, w3.x, w3.y, 0.0f);

      dd2 += FgXXY * w3YZZ + FgiYZZ * w3XXY + Fhi2 * w + ww2 * Fi;
    }
  }

  gpu_accumulate_point<scalar_type, compute_energy, compute_factor, lda>(energy, factor, point_weights,points, 
                                            partial_density, dxyz, dd1, dd2);

}
