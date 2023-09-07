template<class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_calc_VXC(uint points, T* diff, G2G::vec_type<T,WIDTH>* dxyz, 
              G2G::vec_type<T,WIDTH>* diffxyz, T* vrho, T* vsigma, T* v2rho2, T* v2rhosigma, T* v2sigma2,
              // OUTPUTS //
              T* dfac, T* pfac, G2G::vec_type<T,WIDTH>* dfacxyz, G2G::vec_type<T,WIDTH>* pfacxyz,
              int gamma)
{
   uint idx = blockDim.x * blockIdx.x + threadIdx.x;
   bool valid_thread = ( idx < points);

   if ( valid_thread ) {
      double gdens, gdensG;
      gdens  = vrho[idx];
      gdensG = 2.0f*vsigma[idx];
      if (gamma == 1) {
        gdens  = 0.0f;
        gdensG = 0.0f;
      }
      double cruz;
      cruz = dxyz[idx].x*diffxyz[idx].x + dxyz[idx].y*diffxyz[idx].y + dxyz[idx].z*diffxyz[idx].z;

      double perm0, derm0;
      perm0 = vrho[idx];
      derm0 = (v2rho2[idx]*diff[idx] + 2.0f*v2rhosigma[idx]*cruz) * 2.0f;

      double temp;
      temp = 2.0f*vsigma[idx];

      double perm1;
      perm1 = (v2rhosigma[idx]*diff[idx] + 2.0f*v2sigma2[idx]*cruz) * 4.0f;

      // Ground + Diff Excited Densities
      dfac[idx] = gdens + derm0;
      dfacxyz[idx].x = gdensG * dxyz[idx].x + 2.0f*temp*diffxyz[idx].x + perm1 * dxyz[idx].x;
      dfacxyz[idx].y = gdensG * dxyz[idx].y + 2.0f*temp*diffxyz[idx].y + perm1 * dxyz[idx].y;
      dfacxyz[idx].z = gdensG * dxyz[idx].z + 2.0f*temp*diffxyz[idx].z + perm1 * dxyz[idx].z;

      // Diff Excited Density
      pfac[idx] = perm0;
      pfacxyz[idx].x = temp * dxyz[idx].x;
      pfacxyz[idx].y = temp * dxyz[idx].y;
      pfacxyz[idx].z = temp * dxyz[idx].z;
   } // end valid thread
}
