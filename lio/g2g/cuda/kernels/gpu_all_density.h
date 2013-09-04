#define WIDTH 4

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_all_density(scalar_type* const energy, scalar_type* const factor, const scalar_type* const point_weights,
                                    uint points, const scalar_type* rdm, const scalar_type* function_values, const vec_type<scalar_type,4>* gradient_values,
                                    const vec_type<scalar_type,4>* hessian_values, uint m)
{

    //Paralelizacion:
    //threads x : 0 a m = numero de funciones
    //Blocks x: 0 a cantidad de puntos = points

    vec_type<scalar_type,WIDTH>  dxyz, dd1, dd2;
    scalar_type partial_density;

    gpu_compute_density<scalar_type, compute_energy, compute_factor, lda>(energy, factor, point_weights,points, rdm, function_values, gradient_values, hessian_values, m, partial_density, dxyz, dd1, dd2);
/*
    acum_dxyz[threadIdx.x]=dxyz;
    acum_dd1[threadIdx.x]=dd1;
    acum_dd2[threadIdx.x]=dd2;
    acum_pd[threadIdx.x]=partial_density;

    dxyz=reduce(acum_dxyz,sum);
    dd1=reduce(acum_dd1,sum);
    dd2=reduce(acum_dd2,sum);
    partial_density=reduce(acum_pd,sum);

    __syncThreads();
    if(threadIdx.x ==1)
    */
        gpu_accumulate_point<scalar_type, compute_energy, compute_factor, lda>(energy, factor, point_weights,points, partial_density, dxyz, dd1, dd2);
}
