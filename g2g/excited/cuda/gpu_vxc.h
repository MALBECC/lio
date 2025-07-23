template<class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_calc_VXC(uint points, T* diff, G2G::vec_type<T,WIDTH>* dxyz, 
              G2G::vec_type<T,WIDTH>* diffxyz, G2G::vec_type<T,2>* vrho, 
              G2G::vec_type<T,2>* vsigma, G2G::vec_type<T,2>* v2rho2,
              G2G::vec_type<T,2>* v2rhosigma, G2G::vec_type<T,2>* v2sigma2,
              // OUTPUTS //
              T* dfac, T* pfac, G2G::vec_type<T,WIDTH>* dfacxyz, G2G::vec_type<T,WIDTH>* pfacxyz,
              int gamma)
{
   uint idx = blockDim.x * blockIdx.x + threadIdx.x;
   if ( idx < points ) {

      // DIFFERENCE RELAXED DENSITY
      double coef[6];

      double dd = diff[idx];
      double diff_derX, diff_derY, diff_derZ;
      double dens_derX, dens_derY, dens_derZ;

      diff_derX = diffxyz[idx].x;
      diff_derY = diffxyz[idx].y;
      diff_derZ = diffxyz[idx].z;
      dens_derX = dxyz[idx].x * 0.5f;
      dens_derY = dxyz[idx].y * 0.5f;
      dens_derZ = dxyz[idx].z * 0.5f;

      double grad;
      grad=diff_derX*dens_derX+diff_derY*dens_derY+diff_derZ*dens_derZ;

      // GROUND STATE
      double gdens,gdensG1,gdensG2;
      gdens=vrho[idx].x+vrho[idx].y;
      gdensG1=2.0f*(vsigma[idx].x*2.0f+vsigma[idx].y);
      gdensG2=vsigma[idx].y*2.0f;

      // CONTRACTION
      double gdensX,gdensY,gdensZ;
      gdensX=gdensG1*dens_derX+gdensG2*dens_derX;
      gdensY=gdensG1*dens_derY+gdensG2*dens_derY;
      gdensZ=gdensG1*dens_derZ+gdensG2*dens_derZ;

      // V NON CORE CONTRIBUTION
      double dncV1;
      double grad1,grad2,grad3;
      dncV1=vrho[idx].x+vrho[idx].y;
      grad1=2.0f*(vsigma[idx].x*2.0f+vsigma[idx].y);
      grad2=vsigma[idx].y*2.0f;
      grad3=vsigma[idx].y*2.0f;

      double Vncdens,VncdensX,VncdensY,VncdensZ;
      double VncdiffX,VncdiffY,VncdiffZ;
      Vncdens=dncV1;
      VncdensX=grad1*dens_derX+grad2*dens_derX;
      VncdensY=grad1*dens_derY+grad2*dens_derY;
      VncdensZ=grad1*dens_derZ+grad2*dens_derZ;
     
      VncdiffX=grad1*diff_derX+grad3*diff_derX;
      VncdiffY=grad1*diff_derY+grad3*diff_derY;
      VncdiffZ=grad1*diff_derZ+grad3*diff_derZ;

      // V CORE CONTRIBUTION
      coef[0]=dd;
      coef[1]=2.0f*grad;
      coef[2]=grad;
      coef[3]=dd;
      coef[4]=2.0f*grad;
      coef[5]=grad;

      double term1,term2,term3;
      term1=coef[0]*v2rho2[idx].x*2.0f;
      term2=2.0f*coef[0]*(v2rhosigma[idx].x*4.0f+v2rhosigma[idx].y);
      term3=coef[0]*v2rhosigma[idx].y*2.0f;
      term1=term1+coef[1]*(v2rhosigma[idx].x*4.0f+v2rhosigma[idx].y);
      term2=term2+2.0f*coef[1]*(v2sigma2[idx].x*8.0f+v2sigma2[idx].y);
      term3=term3+coef[1]*v2sigma2[idx].y*2.0f;
      term1=term1+coef[2]*v2rhosigma[idx].y*2.0f;
      term2=term2+2.0f*coef[2]*v2sigma2[idx].y*2.0f;
      term3=term3+coef[2]*v2sigma2[idx].y*4.0f;
      term1=term1+coef[3]*v2rho2[idx].y*2.0f;
      term2=term2+2.0f*coef[3]*v2rhosigma[idx].y;
      term3=term3+coef[3]*v2rhosigma[idx].y*2.0f;
      term1=term1+coef[4]*v2rhosigma[idx].y;
      term2=term2+2.0f*coef[4]*v2sigma2[idx].y;
      term3=term3+coef[4]*v2sigma2[idx].y*2.0f;
      term1=term1+coef[5]*v2rhosigma[idx].y*2.0f;
      term2=term2+2.0f*coef[5]*v2sigma2[idx].y*2.0f;
      term3=term3+coef[5]*v2sigma2[idx].y*4.0f;

      // CONTRACTION OF VC
      double Vcdiff,VcdiffX,VcdiffY,VcdiffZ;
      Vcdiff=term1;
      VcdiffX=term2*dens_derX+term3*dens_derX;
      VcdiffY=term2*dens_derY+term3*dens_derY;
      VcdiffZ=term2*dens_derZ+term3*dens_derZ;
      // END V CORE

      if (gamma == 1) {
        gdens=0.0f;
        gdensX=0.0f;
        gdensY=0.0f;
        gdensZ=0.0f;
      }

      dfac[idx]=gdens+Vcdiff;
      dfacxyz[idx].x=gdensX+VncdiffX+VcdiffX;
      dfacxyz[idx].y=gdensY+VncdiffY+VcdiffY;
      dfacxyz[idx].z=gdensZ+VncdiffZ+VcdiffZ;

      // PA
      pfac[idx]=Vncdens;
      pfacxyz[idx].x=VncdensX;
      pfacxyz[idx].y=VncdensY;
      pfacxyz[idx].z=VncdensZ;
   }
}
