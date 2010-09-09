template<bool compute_energy, bool compute_derivs, bool lda>
__global__ void gpu_compute_density(float* const energy, float* const factor, const float* const point_weights,
  uint points, const float* rdm, const float* function_values, const float4* gradient_values, const float4* hessian_values, uint m)
{
  uint point = index_x(blockDim, blockIdx, threadIdx);

  float partial_density = 0.0f;
  float4 dxyz, dd1, dd2;
  if (!lda) { dxyz = dd1 = dd2 = make_float4(0.0f,0.0f,0.0f,0.0f); }

  bool valid_thread = (point < points);
  float point_weight;
  if (valid_thread) point_weight = point_weights[point];

  __shared__ float rdm_sh[DENSITY_BATCH_SIZE];

  /***** compute density ******/
  for (uint i = 0; i < m; i++) {
    float w = 0.0f;
    float4 w3, ww1, ww2;
    if (!lda) { w3 = ww1 = ww2 = make_float4(0.0f,0.0f,0.0f,0.0f); }

    float Fi;
    float4 Fgi, Fhi1, Fhi2;

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
            float4 Fgj = gradient_values[COALESCED_DIMENSION(points) * (bj + j) + point];
            w3 += Fgj * rdm_sh[j];

            float4 Fhj1 = hessian_values[COALESCED_DIMENSION(points) * (2 * (bj + j) + 0) + point];
            float4 Fhj2 = hessian_values[COALESCED_DIMENSION(points) * (2 * (bj + j) + 1) + point];
            ww1 += Fhj1 * rdm_sh[j];
            ww2 += Fhj2 * rdm_sh[j];
          }
        }
      }
    }

    partial_density += Fi * w;

    if (!lda) {
      dxyz += Fgi * w + w3 * Fi;
      dd1 += Fgi * w3 * 2.0f + Fhi1 * w + ww1 * Fi;

      float4 FgXXY = make_float4(Fgi.x, Fgi.x, Fgi.y, 0.0f);
      float4 w3YZZ = make_float4(w3.y, w3.z, w3.z, 0.0f);
      float4 FgiYZZ = make_float4(Fgi.y, Fgi.z, Fgi.z, 0.0f);
      float4 w3XXY = make_float4(w3.x, w3.x, w3.y, 0.0f);

      dd2 += FgXXY * w3YZZ + FgiYZZ * w3XXY + Fhi2 * w + ww2 * Fi;
    }
  }

  /***** compute energy / factor *****/
  float y2a, exc_corr;
  if (compute_energy) {
    if (compute_derivs) {
			gpu_pot<true, true, lda>(partial_density, dxyz, dd1, dd2, exc_corr, y2a);
      if (valid_thread) factor[point] = point_weight * y2a;
		}
		else gpu_pot<true, false, lda>(partial_density, dxyz, dd1, dd2, exc_corr, y2a);

    if (valid_thread) energy[point] = (partial_density * point_weight) * exc_corr;
	}
	else {
		gpu_pot<false, true, lda>(partial_density, dxyz, dd1, dd2, exc_corr, y2a);
    if (valid_thread) factor[point] = point_weight * y2a;
  }
}
