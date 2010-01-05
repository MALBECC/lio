__global__ void gpu_compute_density_derivs(uint points, float* rdmt, float4* gradient_values, float4* density_deriv, uint* nuc,
                                           uint nucleii_count, uint m, float* w)
{
  uint point = index_x(blockDim, blockIdx, threadIdx);
  bool valid_thread = (point < points);

  __shared__ uint nuc_sh[DENSITY_DERIV_BLOCK_SIZE];
  __shared__ float rdm_sh[DENSITY_DERIV_BLOCK_SIZE];

  if (valid_thread) { for (uint i = 0; i < nucleii_count; i++) density_deriv[COALESCED_DIMENSION(points) * i + point] = make_float4(0.0f,0.0f,0.0f,0.0f); }

  for (uint j = 0; j < m; j += DENSITY_DERIV_BLOCK_SIZE) {
    if (j + threadIdx.x < m) nuc_sh[threadIdx.x] = nuc[j + threadIdx.x];

    for (uint jj = 0; jj < DENSITY_DERIV_BLOCK_SIZE && (j + jj < m); jj++) {
      float wrdm = 0.0f;

      for (uint i = 0; i < gpu_nco; i += DENSITY_DERIV_BLOCK_SIZE) {
        if (i + threadIdx.x < gpu_nco) rdm_sh[threadIdx.x] = rdmt[COALESCED_DIMENSION(gpu_nco) * (j + jj) + i + threadIdx.x];

        __syncthreads();
        if (valid_thread) {
          for (uint ii = 0; ii < DENSITY_DERIV_BLOCK_SIZE  && (i + ii < gpu_nco); ii++) {
            wrdm += rdm_sh[ii] * w[COALESCED_DIMENSION(points) * (i + ii) + point];
          }
        }
        __syncthreads();
      }

      if (valid_thread) {
        uint this_nuc = nuc_sh[jj]; // Parece ser necesario para que el compilador coalescee la escritura de abajo
        density_deriv[COALESCED_DIMENSION(points) * this_nuc + point] += gradient_values[COALESCED_DIMENSION(points) * (j + jj) + point] * wrdm;
      }
    }

    __syncthreads();
  }
}

