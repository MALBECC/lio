
#if FULL_DOUBLE
static __inline__ __device__ double fetch_double(texture<int2, 2> t, float x, float y)
{
    int2 v = tex2D(t,x,y);
    return __hiloint2double(v.y, v.x);
}   
#define fetch(t,x,y) fetch_double(t,x,y)
#else
#define fetch(t,x,y) tex2D(t,x,y)
#endif

template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_compute_density(scalar_type* const energy, scalar_type* const factor, const scalar_type* const point_weights,
                                    uint points, const scalar_type* function_values, const vec_type<scalar_type,4>* gradient_values,
                                    const vec_type<scalar_type,4>* hessian_values, uint m, scalar_type* out_partial_density, vec_type<scalar_type,4>* out_dxyz, vec_type<scalar_type,4>* out_dd1, vec_type<scalar_type,4>*  out_dd2
                                    ,int* out_ops)
{

    uint point = blockIdx.x;
    uint i     = threadIdx.x + blockIdx.y * DENSITY_BLOCK_SIZE;

    uint count_sums=0;
    uint count_muls=0;
    uint count_fma =1;

    scalar_type partial_density (0.0f);
    vec_type<scalar_type,WIDTH> dxyz, dd1, dd2;
    dxyz=dd1=dd2 =vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f);

    if (!lda)
    {
        dxyz = dd1 = dd2 = vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f);
    }

    bool valid_thread = (point < points) && ( i < m );


    scalar_type w = 0.0f;
    vec_type<scalar_type,4> w3, ww1, ww2;
    if (!lda)
    {
        w3 = ww1 = ww2 = vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f);
    }

    scalar_type Fi;
    vec_type<scalar_type,4> Fgi, Fhi1, Fhi2;

    int position = threadIdx.x;

    __shared__ scalar_type fj_sh[DENSITY_BLOCK_SIZE];
    __shared__ vec_type<scalar_type, WIDTH> fgj_sh [DENSITY_BLOCK_SIZE];
    __shared__ vec_type<scalar_type, WIDTH> fh1j_sh [DENSITY_BLOCK_SIZE];
    __shared__ vec_type<scalar_type, WIDTH> fh2j_sh [DENSITY_BLOCK_SIZE];

    __shared__ uint ops [DENSITY_BLOCK_SIZE];

    for (int bj = 0; bj <= i; bj += DENSITY_BLOCK_SIZE)
    {
        //Density deberia ser GET_DENSITY_BLOCK_SIZE

        if( (bj+position<m) && (point<points) )
        {
            __syncthreads();

             fj_sh[position] = //  function_values[(m) * point + (bj+position)];
                               fetch(function_values_gpu_tex, (float)(point),(float)(bj+position));

            if(!lda)
            {
                fgj_sh[position] = gradient_values[(m) * point + (bj+position)];

                fh1j_sh[position] = hessian_values[(m)*2 * point +(2 * (bj + position) + 0)];
                fh2j_sh[position] = hessian_values[(m)*2 * point +(2 * (bj + position) + 1)];
            }
        }
        __syncthreads();

        if(valid_thread)
        {
            for(int j=0; j<DENSITY_BLOCK_SIZE && bj+j <= i; j++)
            {
                //fetch es una macro para tex2D
                scalar_type rdm_this_thread = fetch(rmm_input_gpu_tex, (float)(bj+j), (float)i);
                w += rdm_this_thread * fj_sh[j];
                count_sums++;
                count_fma++;
                Fi=fj_sh[j];
                if(!lda)
                {
                    w3 += fgj_sh[j]* rdm_this_thread ;
                    ww1 += fh1j_sh[j] * rdm_this_thread;
                    ww2 += fh2j_sh[j] * rdm_this_thread;
                    Fgi = fgj_sh[j];
                    Fhi1 = fh1j_sh[j] ;
                    Fhi2 = fh2j_sh[j] ;


                    count_fma+=(3*4);

                }
            }
        }
    }
    if(valid_thread)
    {
        partial_density = Fi * w;
        count_muls++;
        //TODO: Insertar aca funcion que convierte <,4> a <,3>
        if (!lda)
        {
            dxyz += Fgi * w + w3 * Fi;
            dd1 += Fgi * w3 * 2.0f + Fhi1 * w + ww1 * Fi;

            vec_type<scalar_type,4> FgXXY(Fgi.x, Fgi.x, Fgi.y, 0.0f);
            vec_type<scalar_type,4> w3YZZ(w3.y, w3.z, w3.z, 0.0f);
            vec_type<scalar_type,4> FgiYZZ(Fgi.y, Fgi.z, Fgi.z, 0.0f);
            vec_type<scalar_type,4> w3XXY(w3.x, w3.x, w3.y, 0.0f);

            dd2 += FgXXY * w3YZZ + FgiYZZ * w3XXY + Fhi2 * w + ww2 * Fi;

            count_sums+=(1+2)*4;
            count_muls+=(1+3)*4;
            count_fma+=(1+3+4)*4;
        }
    }


    __syncthreads();
    //Estamos reutilizando la memoria shared por block para hacer el acumulado por block.
    if(valid_thread)
    {
        fj_sh[position]=partial_density;
        fgj_sh[position]=dxyz;
        fh1j_sh[position]=dd1;
        fh2j_sh[position]=dd2;
        ops[position]=count_sums+count_muls+2*count_fma;
    }
    else
    {
        fj_sh[position]=scalar_type(0.0f);
        fgj_sh[position]=vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f);
        fh1j_sh[position]=vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f);
        fh2j_sh[position]=vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f);
        ops[position]=0;
    }
    __syncthreads();

    for(int j=2;  j <= DENSITY_BLOCK_SIZE ; j=j*2) // 
    {
        int index=position + DENSITY_BLOCK_SIZE/j;
        if( position < DENSITY_BLOCK_SIZE/j)
        {
            fj_sh[position]      += fj_sh[index];
            fgj_sh[position]     += fgj_sh[index];
            fh1j_sh[position]    += fh1j_sh[index];
            fh2j_sh[position]    += fh2j_sh[index];

            ops[position]        += ops[index];
        }
    }
    if(threadIdx.x==0)
    {
        const int myPoint = blockIdx.y*points + blockIdx.x;
        out_partial_density[myPoint] = fj_sh[position];
        out_dxyz[myPoint]            = fgj_sh[position];
        out_dd1[myPoint]             = fh1j_sh[position];
        out_dd2[myPoint]             = fh2j_sh[position];
        
        out_ops[myPoint]            = ops[position];
    }
}

