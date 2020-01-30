template<class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_calc_VXC(uint points, T* diff, G2G::vec_type<T,WIDTH>* dxyz, 
              G2G::vec_type<T,WIDTH>* diffxyz, G2G::vec_type<T,2>* vrho, 
              G2G::vec_type<T,2>* vsigma, G2G::vec_type<T,2>* v2rho2,
              G2G::vec_type<T,2>* v2rhosigma, G2G::vec_type<T,2>* v2sigma2,
              // OUTPUTS //
              T* DDUM, T* PDUM, G2G::vec_type<T,WIDTH>* DDUMxyz, G2G::vec_type<T,WIDTH>* PDUMxyz,
              int gamma)
{
   uint idx = blockDim.x * blockIdx.x + threadIdx.x;
   if ( idx < points ) {

      // DIFFERENCE RELAXED DENSITY
      double C[6];

      double dd = diff[idx];
      double ddx, ddy, ddz;
      double pdx, pdy, pdz;

      ddx = diffxyz[idx].x;
      ddy = diffxyz[idx].y;
      ddz = diffxyz[idx].z;
      pdx = dxyz[idx].x * 0.5f;
      pdy = dxyz[idx].y * 0.5f;
      pdz = dxyz[idx].z * 0.5f;

      double DUMGRV[4], DUMXX[4];
      DUMGRV[0]=ddx*pdx+ddy*pdy+ddz*pdz;
      DUMGRV[1]=ddx*pdx+ddy*pdy+ddz*pdz;
      DUMGRV[2]=ddx*pdx+ddy*pdy+ddz*pdz;
      DUMGRV[3]=ddx*pdx+ddy*pdy+ddz*pdz;

      DUMXX[0]=ddx*ddx+ddy*ddy+ddz*ddz;
      DUMXX[1]=ddx*ddx+ddy*ddy+ddz*ddz;
      DUMXX[2]=ddx*ddx+ddy*ddy+ddz*ddz;
      DUMXX[3]=DUMXX[2];

      // GROUND STATE
      double GDUMA,GDUMAG1,GDUMAG2;
      GDUMA=vrho[idx].x+vrho[idx].y;
      GDUMAG1=2.0f*(vsigma[idx].x*2.0f+vsigma[idx].y);
      GDUMAG2=vsigma[idx].y*2.0f;

      // CONTRACTION
      double GDUMAX,GDUMAY,GDUMAZ;
      GDUMAX=GDUMAG1*pdx+GDUMAG2*pdx;
      GDUMAY=GDUMAG1*pdy+GDUMAG2*pdy;
      GDUMAZ=GDUMAG1*pdz+GDUMAG2*pdz;

      // V NON CORE CONTRIBUTION
      double DUMNV1;
      double DUMGRV1,DUMGRV3,DUMGRV4;
      DUMNV1=vrho[idx].x+vrho[idx].y;
      DUMGRV1=2.0f*(vsigma[idx].x*2.0f+vsigma[idx].y);
      DUMGRV3=vsigma[idx].y*2.0f;
      DUMGRV4=vsigma[idx].y*2.0f;

      double VNCDOMA,VNCDOMAX,VNCDOMAY,VNCDOMAZ;
      double VNCDUMAX,VNCDUMAY,VNCDUMAZ;
      VNCDOMA=DUMNV1;
      VNCDOMAX=DUMGRV1*pdx+DUMGRV3*pdx;
      VNCDOMAY=DUMGRV1*pdy+DUMGRV3*pdy;
      VNCDOMAZ=DUMGRV1*pdz+DUMGRV3*pdz;
     
      VNCDUMAX=DUMGRV1*ddx+DUMGRV4*ddx;
      VNCDUMAY=DUMGRV1*ddy+DUMGRV4*ddy;
      VNCDUMAZ=DUMGRV1*ddz+DUMGRV4*ddz;

      // V CORE CONTRIBUTION
      C[0]=dd;
      C[1]=2.0f*DUMGRV[0];
      C[2]=DUMGRV[2];
      C[3]=dd;
      C[4]=2.0f*DUMGRV[1];
      C[5]=DUMGRV[3];

      double DUMA,DUMAG1,DUMAG2;
      DUMA=C[0]*v2rho2[idx].x*2.0f;
      DUMAG1=2.0f*C[0]*(v2rhosigma[idx].x*4.0f+v2rhosigma[idx].y);
      DUMAG2=C[0]*v2rhosigma[idx].y*2.0f;
      DUMA=DUMA+C[1]*(v2rhosigma[idx].x*4.0f+v2rhosigma[idx].y);
      DUMAG1=DUMAG1+2.0f*C[1]*(v2sigma2[idx].x*8.0f+v2sigma2[idx].y);
      DUMAG2=DUMAG2+C[1]*v2sigma2[idx].y*2.0f;
      DUMA=DUMA+C[2]*v2rhosigma[idx].y*2.0f;
      DUMAG1=DUMAG1+2.0f*C[2]*v2sigma2[idx].y*2.0f;
      DUMAG2=DUMAG2+C[2]*v2sigma2[idx].y*4.0f;
      DUMA=DUMA+C[3]*v2rho2[idx].y*2.0f;
      DUMAG1=DUMAG1+2.0f*C[3]*v2rhosigma[idx].y;
      DUMAG2=DUMAG2+C[3]*v2rhosigma[idx].y*2.0f;
      DUMA=DUMA+C[4]*v2rhosigma[idx].y;
      DUMAG1=DUMAG1+2.0f*C[4]*v2sigma2[idx].y;
      DUMAG2=DUMAG2+C[4]*v2sigma2[idx].y*2.0f;
      DUMA=DUMA+C[5]*v2rhosigma[idx].y*2.0f;
      DUMAG1=DUMAG1+2.0f*C[5]*v2sigma2[idx].y*2.0f;
      DUMAG2=DUMAG2+C[5]*v2sigma2[idx].y*4.0f;

      // CONTRACTION OF VC
      double VCDUMA,VCDUMAX,VCDUMAY,VCDUMAZ;
      VCDUMA=DUMA;
      VCDUMAX=DUMAG1*pdx+DUMAG2*pdx;
      VCDUMAY=DUMAG1*pdy+DUMAG2*pdy;
      VCDUMAZ=DUMAG1*pdz+DUMAG2*pdz;
      // END V CORE

      if (gamma == 1) {
        GDUMA=0.0f;
        GDUMAX=0.0f;
        GDUMAY=0.0f;
        GDUMAZ=0.0f;
      }

      DDUM[idx]=GDUMA+VCDUMA;
      DDUMxyz[idx].x=GDUMAX+VNCDUMAX+VCDUMAX;
      DDUMxyz[idx].y=GDUMAY+VNCDUMAY+VCDUMAY;
      DDUMxyz[idx].z=GDUMAZ+VNCDUMAZ+VCDUMAZ;

      // PA
      PDUM[idx]=VNCDOMA;
      PDUMxyz[idx].x=VNCDOMAX;
      PDUMxyz[idx].y=VNCDOMAY;
      PDUMxyz[idx].z=VNCDOMAZ;
   }
}
