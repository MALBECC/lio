template <class scalar_type, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_obtain_fock(uint npoints, uint M, const scalar_type* const weights,
                               const scalar_type* function_values, 
                               const scalar_type* factor,
                               scalar_type* Fock)
{

   int temp = int(sqrtf(1 + 8 * blockIdx.x));
   temp -= (1 - temp % 2);
   int block_row = ( temp - 1 ) / 2;
   int block_col = blockIdx.x - ( (block_row+1) * block_row / 2 );
   int first_row = block_row * blockDim.x;
   int first_col = block_col * blockDim.y;
   int row = first_row + threadIdx.x;
   int col = first_col + threadIdx.y;
   bool valid_thread = ( row < M && col < M && row >= col );

   scalar_type term1, term2, term3, term4;
   scalar_type result=0.0f;
   int tot_p = COALESCED_DIMENSION(npoints);

   if ( valid_thread ) {
     for(uint point=0; point<npoints;point++) {
        term1 = function_values[point*M+col]*factor[point*M+row];
        term2 = function_values[point*M+row]*factor[point*M+col];
        result += (term1+term2);
     }
     Fock[row*M+col] = result;
   } // end valid thread
}

