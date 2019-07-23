#define WIDTH 4
template <class scalar_type>
__global__ void gpu_compute_becke_os(scalar_type* becke_dens, scalar_type* becke_spin,
                                     const scalar_type* partial_density_a,
                                     const scalar_type* partial_density_b,
                                     const scalar_type* point_weights,
                                     const scalar_type* becke_w, uint points,
                                     uint atoms, int block_height) {
  
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  scalar_type _partial_density_a(0.0f);
  scalar_type _partial_density_b(0.0f);

  // Checks if we are in a valid thread.
  if (!(point < points)) return;
    
  for (int j = 0; j < block_height; j++) {
    const int this_row = j * points + point;
    _partial_density_a += partial_density_a[this_row];
    _partial_density_b += partial_density_b[this_row];
  }

  for (int j = 0; j < atoms; j++) {
    becke_dens[point*atoms +j] += (_partial_density_b + _partial_density_a)
                                  * point_weights[point] * becke_w[point*atoms +j];
    becke_spin[point*atoms +j] += (_partial_density_b - _partial_density_a) 
                                  * point_weights[point] * becke_w[point*atoms +j];
  }
}

template <class scalar_type>
__global__ void gpu_compute_becke_cs(scalar_type* becke_dens,
                                     const scalar_type* partial_density,
                                     const scalar_type* point_weights, 
                                     const scalar_type* becke_w, uint points,
                                     uint atoms, int block_height) {
  
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

template <class scalar_type>
__global__ void gpu_cdft_factors(scalar_type* factors, const uint* reg_natom,
                                 const uint* reg_atoms,
                                 const scalar_type* point_weights, 
                                 const scalar_type* becke_w, uint points,
                                 uint atoms, uint regions, uint max_nat) {
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  scalar_type _accum_factor;
  // Checks if thread is valid.
  if (!(point < points)) return;
  for (int j = 0; j < regions     ; j++) {
    _accum_factor = (scalar_type) 0.0f;
    for (int i = 0; i < reg_natom[j]; i++) {
      _accum_factor += becke_w[point*atoms + reg_atoms[j*max_nat + i]];
    }
    factors[point*regions +j] = point_weights[point] * _accum_factor;
  }
}

template <class scalar_type>
__global__ void gpu_cdft_factors_accum(scalar_type* factors, uint points,
                                       uint regions, scalar_type* pot,
                                       scalar_type* newfacs) {
  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  // Checks if thread is valid.
  if (!(point < points)) return;
  for (int j = 0; j < regions; j++) {
    newfacs[point] += factors[point*regions + j] * pot[j];
  }
}
