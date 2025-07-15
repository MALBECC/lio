template<class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_calc_FXC(uint points, T* dens, T* trad, G2G::vec_type<T,WIDTH>* dxyz,
                G2G::vec_type<T,WIDTH>* tradxyz, T* vsigma,
                T* v2rho2, T* v2rhosigma, T* v2sigma2, T* v3rho3,
                T* v3rho2sigma, T* v3rhosigma2, T* v3sigma3, 
                // OUTPUTS //
                T* dfac, T* tfac, G2G::vec_type<T,WIDTH>* dfacxyz, G2G::vec_type<T,WIDTH>* tfacxyz)
{
   uint idx = blockDim.x * blockIdx.x + threadIdx.x;
   bool valid_thread = ( idx < points);

   if ( valid_thread ) {
      double cruz, contr;
      cruz  = tradxyz[idx].x*dxyz[idx].x + tradxyz[idx].y*dxyz[idx].y + tradxyz[idx].z*dxyz[idx].z;
      contr = tradxyz[idx].x*tradxyz[idx].x + tradxyz[idx].y*tradxyz[idx].y + tradxyz[idx].z*tradxyz[idx].z;

      double term0;
      term0 = (2.0f*v2rho2[idx]*trad[idx]+ 4.0f*v2rhosigma[idx]*cruz) * 2.0f;

      double derm0, derm1;
      derm0  = trad[idx]*trad[idx]*v3rho3[idx]+4.0f*cruz*trad[idx]*v3rho2sigma[idx];
      derm0 += 4.0f*cruz*cruz*v3rhosigma2[idx] + 2.0f*contr*v2rhosigma[idx];
      derm0 *= 4.0f;

      derm1  = trad[idx]*trad[idx]*v3rho2sigma[idx] + 4.0f*v3rhosigma2[idx]*cruz*trad[idx];
      derm1 += 4.0f* v3sigma3[idx]*cruz*cruz + 2.0f*v2sigma2[idx]*contr;
      derm1 *= 8.0f;

      double temp0, temp1;
      temp0  = (4.0f*v2rhosigma[idx]*trad[idx] + 8.0f*v2sigma2[idx]*cruz) * 4.0f;
      temp1  = (4.0f*vsigma[idx]) * 4.0f;

      // DENSITY FACTOR
      dfac[idx] += derm0;
      dfacxyz[idx].x += derm1*dxyz[idx].x + temp0*tradxyz[idx].x;
      dfacxyz[idx].y += derm1*dxyz[idx].y + temp0*tradxyz[idx].y;
      dfacxyz[idx].z += derm1*dxyz[idx].z + temp0*tradxyz[idx].z;

      // TRANSITION DENSITY FACTOR
      tfac[idx] = term0;
      tfacxyz[idx].x = 0.5f*(temp0*dxyz[idx].x+temp1*tradxyz[idx].x);
      tfacxyz[idx].y = 0.5f*(temp0*dxyz[idx].y+temp1*tradxyz[idx].y);
      tfacxyz[idx].z = 0.5f*(temp0*dxyz[idx].z+temp1*tradxyz[idx].z);
   } // end valid thread
}
