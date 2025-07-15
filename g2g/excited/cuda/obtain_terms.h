template <class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_obtain_term(uint npoints, uint M, const scalar_type* const weights,
                               const scalar_type* function_values,
                               const vec_type<scalar_type, 4>* gradient_values,
                               const scalar_type* factorF,const vec_type<scalar_type, 4>* factorD,
                               const vec_type<scalar_type, 4>* factorT,
                               scalar_type* terms)
{
   uint point = blockIdx.x;
   uint func  = threadIdx.x;
   bool valid_thread = (point < npoints) && (func < M);
 
   __shared__ scalar_type wp;
   scalar_type term1, term2, term3, term4;

   if (func == 0) wp = weights[point];
   __syncthreads();

   if ( valid_thread ) {
      term1  = factorF[point]   * function_values[point*M+func] * 0.5f;
      term1 += factorD[point].x * gradient_values[point*M+func].x;

      term2  = factorD[point].y * gradient_values[point*M+func].y;
      term2 += factorD[point].z * gradient_values[point*M+func].z;

      term3  = factorT[point].x * gradient_values[point*M+func].x;
      term3 += factorT[point].y * gradient_values[point*M+func].y;

      term4  = factorT[point].z * gradient_values[point*M+func].z;

      terms[point*M+func]   = (term1+term2+term3+term4) * wp;
   }
}

