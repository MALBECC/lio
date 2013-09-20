/*
__device__ doubleatomicAdd(vec_type<scalar_type,WIDTH>* address, doubleval)
{
    unsigned long long int* address_as_ull =(unsigned long long int*)address;
    unsigned long long intold = *address_as_ull, assumed;
    do{
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val +__longlong_as_double(assumed)));
    }while(assumed != old);
    return __longlong_as_double(old);
}
__device__ doubleatomicAdd(scalar_type* address, scalar_type doubleval)
{
    unsigned long long int* address_as_ull =(unsigned long long int*)address;
    unsigned long long intold = *address_as_ull, assumed;
    do{
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val +__longlong_as_double(assumed)));
    }while(assumed != old);
    return __longlong_as_double(old);
}
*/


template<class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_compute_density(scalar_type* const energy, scalar_type* const factor, const scalar_type* const point_weights,
                                    uint points, const scalar_type* rdm, const scalar_type* function_values, const vec_type<scalar_type,4>* gradient_values,
                                    const vec_type<scalar_type,4>* hessian_values, uint m, scalar_type* out_partial_density, vec_type<scalar_type,4>* out_dxyz, vec_type<scalar_type,4>* out_dd1, vec_type<scalar_type,4>*  out_dd2)
{

  uint point = blockIdx.x;
  /* Old */
  //uint i     = threadIdx.x;
  /* New */
  uint i     = threadIdx.x + blockIdx.y * DENSITY_BLOCK_SIZE;

  scalar_type partial_density (0.0f);
  vec_type<scalar_type,WIDTH> dxyz, dd1, dd2;
  dxyz=dd1=dd2 =vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f); 

  if (!lda) { dxyz = dd1 = dd2 = vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f); }

    bool valid_thread = (point < points) && ( i < m );
 

    scalar_type w = 0.0f;
    vec_type<scalar_type,4> w3, ww1, ww2;
    if (!lda) { w3 = ww1 = ww2 = vec_type<scalar_type,4>(0.0f,0.0f,0.0f,0.0f); }

    scalar_type Fi;
    vec_type<scalar_type,4> Fgi, Fhi1, Fhi2;

    //TODO: Cada thread del güarp trae su Fi.
    if (valid_thread) {
      Fi = function_values[COALESCED_DIMENSION(points) * i + point]; //Con la paralelizacion a nivel de thread, esta coalescencia desaparece. Hay que darlo vuelta / transpose.
      if (!lda) {
        Fgi = gradient_values[COALESCED_DIMENSION(points) * i + point];  //Deberia ser: Coalesced_dimension(i) * point + i 
        Fhi1 = hessian_values[COALESCED_DIMENSION(points) * (2 * i + 0) + point];   //Hay que cambiarlo de functions.h
        Fhi2 = hessian_values[COALESCED_DIMENSION(points) * (2 * i + 1) + point];
      }
    }
    int position = threadIdx.x;

    __shared__ scalar_type fj_sh[DENSITY_BLOCK_SIZE];
    __shared__ vec_type<scalar_type, WIDTH> fgj_sh [DENSITY_BLOCK_SIZE];
    __shared__ vec_type<scalar_type, WIDTH> fh1j_sh [DENSITY_BLOCK_SIZE];
    __shared__ vec_type<scalar_type, WIDTH> fh2j_sh [DENSITY_BLOCK_SIZE];
    
    if(valid_thread)
    {
        for (int bj = 0; bj <= i; bj += DENSITY_BLOCK_SIZE) 
        { //Density deberia ser GET_DENSITY_BLOCK_SIZE
     
            scalar_type valid_position = (1.0f * (bj+position<=i)) ;
            fj_sh[position] = function_values[COALESCED_DIMENSION(points) * (bj + position) + point] * valid_position;
            if(!lda)
            {
                fgj_sh[position] = vec_type<scalar_type, WIDTH> (gradient_values[COALESCED_DIMENSION(points) * (bj + position) + point] *valid_position);

                fh1j_sh[position] = vec_type<scalar_type,WIDTH>(hessian_values[COALESCED_DIMENSION(points) * (2 * (bj + position) + 0) + point]*valid_position);
                fh2j_sh[position] = vec_type<scalar_type,WIDTH>(hessian_values[COALESCED_DIMENSION(points) * (2 * (bj + position) + 1) + point]*valid_position);
            }
            __syncthreads();

            for(int j=0; j<DENSITY_BLOCK_SIZE && bj+j <= i; j++){            
                scalar_type rdm_this_thread = rdm[COALESCED_DIMENSION(m) * i + (bj+j)];
                w += rdm_this_thread * fj_sh[j];

                if(!lda)
                {
                    w3 += fgj_sh[j]* rdm_this_thread ;
                    ww1 += fh1j_sh[j] * rdm_this_thread;
                    ww2 += fh2j_sh[j] * rdm_this_thread;
                }
            }

        }

        partial_density = Fi * w;
        //TODO: Insertar aca funcion que convierte <,4> a <,3>
        if (!lda) {
          dxyz += Fgi * w + w3 * Fi;
          dd1 += Fgi * w3 * 2.0f + Fhi1 * w + ww1 * Fi;

          vec_type<scalar_type,4> FgXXY(Fgi.x, Fgi.x, Fgi.y, 0.0f);
          vec_type<scalar_type,4> w3YZZ(w3.y, w3.z, w3.z, 0.0f);
          vec_type<scalar_type,4> FgiYZZ(Fgi.y, Fgi.z, Fgi.z, 0.0f);
          vec_type<scalar_type,4> w3XXY(w3.x, w3.x, w3.y, 0.0f);

          dd2 += FgXXY * w3YZZ + FgiYZZ * w3XXY + Fhi2 * w + ww2 * Fi;
        }
    }

/*
    __shared__ vec_type<scalar_type,WIDTH> _dxyz, _dd1, _dd2;
    __shared__ scalar_type _pd (0.0f);
    dxyz=dd1=dd2=vec_type<scalar_type,WIDTH>(0.0f,0.0f,0.0f,0.0f);
    
    doubleatomicAdd(_pd, partial_density);
    vec_typeatomicAdd(_dxyz, dxyz);
    vec_typeatomicAdd(_dd1, dd1);
    vec_typeatomicAdd(_dd2, dd2);
*/
    __syncthreads();
    //Estamos reutilizando la memoria shared por block para hacer el acumulado por block.
    fj_sh[position]=partial_density;
    fgj_sh[position]=dxyz;
    fh1j_sh[position]=dd1;
    fh2j_sh[position]=dd2;

    __syncthreads();

    if(threadIdx.x==0)
    {
        dxyz=dd1=dd2=vec_type<scalar_type,WIDTH>(0.0f,0.0f,0.0f,0.0f);
        partial_density=scalar_type(0.0f);

        for(int j=0; (blockIdx.y < m/DENSITY_BLOCK_SIZE && j<DENSITY_BLOCK_SIZE) 
                         || 
                    (blockIdx.y == m/DENSITY_BLOCK_SIZE && j<(m % DENSITY_BLOCK_SIZE)); j++)
        {
            
            partial_density     += fj_sh[j];
            dxyz                += fgj_sh[j];
            dd1                 += fh1j_sh[j];
            dd2                 += fh2j_sh[j];
        }        
        const int myPoint = blockIdx.y*points + blockIdx.x;
        //printf("BLKX: %d  BLKY: %d   --- PD: %e\n", blockIdx.x, blockIdx.y, partial_density);
        out_partial_density[myPoint] = partial_density;
        out_dxyz[myPoint] = dxyz;
        out_dd1[myPoint] = dd1;
        out_dd2[myPoint] = dd2;
    }
}

