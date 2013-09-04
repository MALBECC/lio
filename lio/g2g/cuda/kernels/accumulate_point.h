#define WIDTH 4

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__device__ void gpu_accumulate_point(scalar_type* const energy, scalar_type* const factor, const scalar_type* const point_weights,
                                    uint points, scalar_type partial_density, vec_type<scalar_type,WIDTH>  dxyz, 
                                    vec_type<scalar_type,WIDTH> dd1, vec_type<scalar_type,WIDTH> dd2){

  //gpu_accumulate_point<scalar_type, compute_energy, true, lda>(energy, factor, point_weights,points, 
  //                                        partial_density, dxyz, dd1, dd2);

  uint point = index_x(blockDim, blockIdx, threadIdx);
  scalar_type point_weight = 0.0f;
  scalar_type y2a, exc_corr;

  bool valid_thread = (point < points);
  if (valid_thread) point_weight = point_weights[point];

  gpu_pot<scalar_type, compute_energy, true, lda>(partial_density, dxyz, dd1, dd2, exc_corr, y2a); // TODO: segundo parametro, que tenga en cuenta RMM+energy
  
  if (compute_energy && valid_thread)
    energy[point] = (partial_density * point_weight) * exc_corr;
  
  if (compute_factor && valid_thread)
    factor[point] = point_weight * y2a;

}
