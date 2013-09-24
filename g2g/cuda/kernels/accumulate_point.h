#define WIDTH 4

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_accumulate_point(scalar_type* const energy, scalar_type* const factor, const scalar_type* const point_weights,
                                    uint points, int block_height, scalar_type* partial_density, vec_type<scalar_type,WIDTH>*  dxyz, 
                                    vec_type<scalar_type,WIDTH>* dd1, vec_type<scalar_type,WIDTH>* dd2){


  uint point = blockIdx.x * DENSITY_ACCUM_BLOCK_SIZE + threadIdx.x;

  
  scalar_type point_weight = 0.0f;
  scalar_type y2a, exc_corr;

  scalar_type _partial_density(0.0f);
  vec_type<scalar_type,WIDTH> _dxyz, _dd1, _dd2;
  _dxyz=_dd1=_dd2=vec_type<scalar_type,WIDTH>(0.0f,0.0f,0.0f,0.0f);

  bool valid_thread = (point < points);
  if (valid_thread) point_weight = point_weights[point];

  if (valid_thread)
  {
     for(int j =0 ; j<block_height; j++)
     {
         const int this_row=j*points+point;
         _partial_density += partial_density[this_row];
         _dxyz += dxyz[this_row];
         _dd1 += dd1[this_row];
         _dd2 += dd2[this_row];
     }
  }

  gpu_pot<scalar_type, compute_energy, true, lda>(_partial_density, _dxyz, _dd1, _dd2, exc_corr, y2a); // TODO: segundo parametro, que tenga en cuenta RMM+energy

  if (compute_energy && valid_thread)
    energy[point] = (_partial_density * point_weight) * exc_corr;
  if (compute_factor && valid_thread)
    factor[point] = point_weight * y2a;

}
