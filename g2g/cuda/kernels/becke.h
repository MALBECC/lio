#define WIDTH 4
template <class scalar_type>
__global__ void gpu_compute_becke_os(scalar_type* becke_dens, const scalar_type* partial_density_a,
                                     const scalar_type* partial_density_b ,
                                     const scalar_type* point_weights, const scalar_type* becke_w,
                                     uint points, uint atoms, int block_height) {
  
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  scalar_type _partial_density(0.0f);

  // Checks if we are in a valid thread.
  if (!(point < points)) return;
    
  for (int j = 0; j < block_height; j++) {
    const int this_row = j * points + point;
    _partial_density += partial_density_a[this_row];
    _partial_density += partial_density_b[this_row];
  }
  
  for (int j = 0; j < atoms; j++) {
    becke_dens[point*atoms +j] += _partial_density * point_weights[point] * becke_w[point*atoms +j];
  }
}

template <class scalar_type>
__global__ void gpu_compute_becke_cs(scalar_type* becke_dens, const scalar_type* partial_density,
                                     const scalar_type* point_weights, const scalar_type* becke_w,
                                     uint points, uint atoms, int block_height) {
  
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  scalar_type _partial_density(0.0f);

  // Checks if we are in a valid thread.
  if (!(point < points)) return;
    
  for (int j = 0; j < block_height; j++) {
    const int this_row = j * points + point;
    _partial_density += partial_density[this_row];
  }
  
  for (int j = 0; j < atoms; j++) {
    becke_dens[point*atoms +j] += _partial_density * point_weights[point] * becke_w[point*atoms +j];
  }
}
