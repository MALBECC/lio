#define WIDTH 4

    extern  __shared__ void* acum[];
template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_all_density(scalar_type* const energy, scalar_type* const factor, const scalar_type* const point_weights,
                                    uint points, const scalar_type* rdm, const scalar_type* function_values, const vec_type<scalar_type,4>* gradient_values,
                                    const vec_type<scalar_type,4>* hessian_values, uint m)
{

    //Paralelizacion:
    //threads x : 0 a m = numero de funciones
    //Blocks x: 0 a cantidad de puntos = points

    vec_type<scalar_type, WIDTH>* acum_density = (vec_type<scalar_type, WIDTH>*) acum;

    vec_type<scalar_type,WIDTH>  dxyz, dd1, dd2;
    scalar_type partial_density;

    gpu_compute_density<scalar_type, compute_energy, compute_factor, lda>(energy, factor, point_weights,points, rdm, function_values, gradient_values, hessian_values, m, partial_density, dxyz, dd1, dd2);

    acum_density[threadIdx.x*4+0]=dxyz;
    acum_density[threadIdx.x*4+1]=dd1;
    acum_density[threadIdx.x*4+2]=dd2;
    acum_density[threadIdx.x*4+3]=vec_type<scalar_type,WIDTH>(partial_density,0.0f,0.0f,0.0f);
    __syncthreads();
    
    
    if(threadIdx.x ==0)  /* TODO: Esta acumulacion deberia ser realizada por varios threads, no por solo el 1 */
    {       
        dxyz=dd1=dd2=vec_type<scalar_type,WIDTH>(0.0f,0.0f,0.0f,0.0f);
        partial_density=scalar_type(0.0f);
        for(int i =0; i<blockDim.x*4; i+=4)
        {
            dxyz+=acum_density[i+0];
            dd1+=acum_density[i+1];
            dd2+=acum_density[i+2];
            partial_density+=acum_density[i+3].x;
            
        }
        gpu_accumulate_point<scalar_type, compute_energy, compute_factor, lda>(energy, factor, point_weights,points, partial_density, dxyz, dd1, dd2);

    }
    
}
