#define WIDTH 4

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_all_density(scalar_type* const energy, scalar_type* const factor, const scalar_type* const point_weights,
                                    uint points, const scalar_type* rdm, const scalar_type* function_values, const vec_type<scalar_type,4>* gradient_values,
                                    const vec_type<scalar_type,4>* hessian_values, uint m)
{
    vec_type<scalar_type,WIDTH>  dxyz, dd1, dd2;
    scalar_type partial_density;
    gpu_compute_density<scalar_type, compute_energy, compute_factor, lda>(energy, factor, point_weights,points, rdm, function_values, gradient_values, hessian_values, m, partial_density, dxyz, dd1, dd2);
    gpu_accumulate_point<scalar_type, compute_energy, compute_factor, lda>(energy, factor, point_weights,points, partial_density, dxyz, dd1, dd2);
}
