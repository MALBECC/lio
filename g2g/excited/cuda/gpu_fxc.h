template<class T, bool compute_energy, bool compute_factor, bool lda>
__global__ void gpu_calc_FXC(uint points, T* dens, T* trad, G2G::vec_type<T,WIDTH>* dxyz,
                G2G::vec_type<T,WIDTH>* tradxyz, G2G::vec_type<T,2>* vsigma,
                G2G::vec_type<T,2>* v2rho2, G2G::vec_type<T,2>* v2rhosigma,
                G2G::vec_type<T,2>* v2sigma2, G2G::vec_type<T,2>* v3rho3,
                G2G::vec_type<T,2>* v3rho2sigma, G2G::vec_type<T,2>* v3rhosigma2,
                G2G::vec_type<T,2>* v3sigma3, 
                // OUTPUTS //
                T* DDUM, T* VDUM, G2G::vec_type<T,WIDTH>* DDUMxyz, G2G::vec_type<T,WIDTH>* VDUMxyz)
{
   uint idx = blockDim.x * blockIdx.x + threadIdx.x;

   if ( idx < points ) {
      // TRANSITION DENSITY
      double DUMNV[2], DXV[2], DYV[2], DZV[2];
      DUMNV[0]=trad[idx];
      DXV[0]=tradxyz[idx].x;
      DYV[0]=tradxyz[idx].y;
      DZV[0]=tradxyz[idx].z;
      DUMNV[1]=trad[idx];
      DXV[1]=tradxyz[idx].x;
      DYV[1]=tradxyz[idx].y;
      DZV[1]=tradxyz[idx].z;
      double ddx, ddy, ddz;
      ddx = dxyz[idx].x*0.5f;
      ddy = dxyz[idx].y*0.5f;
      ddz = dxyz[idx].z*0.5f;

      double DUMGRV[4], DUMXX[4];
      DUMGRV[0]=DXV[0]*ddx+DYV[0]*ddy+DZV[0]*ddz;
      DUMGRV[1]=DXV[1]*ddx+DYV[1]*ddy+DZV[1]*ddz;
      DUMGRV[2]=DXV[0]*ddx+DYV[0]*ddy+DZV[0]*ddz;
      DUMGRV[3]=DXV[1]*ddx+DYV[1]*ddy+DZV[1]*ddz;
      DUMXX[0]=DXV[0]*DXV[0]+DYV[0]*DYV[0]+DZV[0]*DZV[0];
      DUMXX[1]=DXV[1]*DXV[1]+DYV[1]*DYV[1]+DZV[1]*DZV[1];
      DUMXX[2]=DXV[0]*DXV[1]+DYV[0]*DYV[1]+DZV[0]*DZV[1];
      DUMXX[3]=DUMXX[2];
   
      // F NON CORE
      double C[20];
      C[0]=2.0f*(vsigma[idx].x*2.0f+vsigma[idx].y);
      C[1]=v2rho2[idx].x*2.0f+v2rho2[idx].y;
      C[2]=2.0f*(v2rhosigma[idx].x*4.0f+v2rhosigma[idx].y);
      C[3]=2.0f*(v2rhosigma[idx].x*4.0f+v2rhosigma[idx].y);
      C[4]=v2rhosigma[idx].y*2.0f;
      C[5]=v2rhosigma[idx].y*2.0f;
      C[6]=4.0f*(v2sigma2[idx].x*8.0f+v2sigma2[idx].y);
      C[7]=2.0f*v2sigma2[idx].y*2.0f;
      C[8]=2.0f*v2sigma2[idx].y*2.0f;
      C[9]=v2sigma2[idx].y*4.0f;

      double DUMC, DUMXX1;
      DUMC=C[0];
      DUMXX1=DUMC*2.0f;
      double DUMNV1;
      DUMC=C[1];
      DUMNV1=DUMC*DUMNV[0]*2.0f;

      double DUMGRV1;
      DUMC=C[2]+C[3];
      DUMNV1=DUMNV1+DUMC*DUMGRV[0];
      DUMGRV1=DUMC*DUMNV[0];

      double DUMGRV3;
      DUMC=C[4]+C[5];
      DUMNV1=DUMNV1+DUMC*DUMGRV[2];
      DUMGRV3=DUMC*DUMNV[0];

      DUMC=C[6];
      DUMGRV1=DUMGRV1+C[6]*DUMGRV[0]*2.0f;
      DUMC=C[7]+C[8];
      DUMGRV1=DUMGRV1+DUMC*DUMGRV[2];
      DUMGRV3=DUMGRV3+DUMC*DUMGRV[0];
      DUMC=C[9];
      DUMGRV3=DUMGRV3+DUMC*DUMGRV[2]*2.0f;
  
      // --AB AND BA COMPONENTS--
      C[10]=vsigma[idx].y*2.0f;
      C[11]=v2rho2[idx].y;
      C[12]=2.0f*v2rhosigma[idx].y;
      C[13]=2.0f*v2rhosigma[idx].y;
      C[14]=v2rhosigma[idx].y*2.0f;
      C[15]=v2rhosigma[idx].y*2.0f;
      C[16]=4.0f*v2sigma2[idx].y;
      C[17]=2.0f*v2sigma2[idx].y*2.0f;
      C[18]=2.0f*v2sigma2[idx].y*2.0f;
      C[19]=v2sigma2[idx].y*4.0f;

      double DUMXX3, DUMNV2;
      DUMC=C[10]*2.0f;
      DUMXX3=DUMC;
      DUMC=C[11]*2.0f;
      DUMNV1=DUMNV1+DUMC*DUMNV[1];
      DUMNV2=DUMC*DUMNV[0];

      double DUMGRV2;
      DUMC=C[12]*2.0f;
      DUMNV1=DUMNV1+DUMC*DUMGRV[1];
      DUMGRV2=DUMC*DUMNV[0];
      DUMC=C[13]*2.0f;
      DUMNV2=DUMNV2+DUMC*DUMGRV[0];
      DUMGRV1=DUMGRV1+DUMC*DUMNV[1];

      double DUMGRV4;
      DUMC=C[14]*2.0f;
      DUMNV1=DUMNV1+DUMC*DUMGRV[3];
      DUMGRV4=DUMC*DUMNV[0];
      DUMC=C[15]*2.0f;
      DUMNV2=DUMNV2+DUMC*DUMGRV[2];
      DUMGRV3=DUMGRV3+DUMC*DUMNV[1];

      DUMC=C[16]*2.0f;
      DUMGRV1=DUMGRV1+DUMC*DUMGRV[1];
      DUMGRV2=DUMGRV2+DUMC*DUMGRV[0];
      DUMC=C[17]*2.0f;
      DUMGRV1=DUMGRV1+DUMC*DUMGRV[3];
      DUMGRV4=DUMGRV4+DUMC*DUMGRV[0];
      DUMC=C[18]*2.0f;
      DUMGRV2=DUMGRV2+DUMC*DUMGRV[2];
      DUMGRV3=DUMGRV3+DUMC*DUMGRV[1];
      DUMC=C[19]*2.0f;
      DUMGRV3=DUMGRV3+DUMC*DUMGRV[3];
      DUMGRV4=DUMGRV4+DUMC*DUMGRV[2];

      // --BB COMPONENTS--
      C[0]=2.0f*(vsigma[idx].x*2.0f+vsigma[idx].y);
      C[1]=v2rho2[idx].x*2.0f+v2rho2[idx].y;
      C[2]=2.0f*(v2rhosigma[idx].x*4.0f+v2rhosigma[idx].y);
      C[3]=2.0f*(v2rhosigma[idx].x*4.0f+v2rhosigma[idx].y);
      C[4]=v2rhosigma[idx].y*2.0f;
      C[5]=v2rhosigma[idx].y*2.0f;
      C[6]=4.0f*(v2sigma2[idx].x*8.0f+v2sigma2[idx].y);
      C[7]=2.0f*v2sigma2[idx].y*2.0f;
      C[8]=2.0f*v2sigma2[idx].y*2.0f;
      C[9]=v2sigma2[idx].y*4.0f;

      DUMC=C[0];
      DUMC=C[1];
      DUMNV2=DUMNV2+DUMC*DUMNV[1]*2.0f;
      DUMC=C[2]+C[3];
      DUMNV2=DUMNV2+DUMC*DUMGRV[1];
      DUMGRV2=DUMGRV2+DUMC*DUMNV[1];
      DUMC=C[4]+C[5];
      DUMNV2=DUMNV2+DUMC*DUMGRV[3];
      DUMGRV4=DUMGRV4+DUMC*DUMNV[1];
      DUMC=C[6];
      DUMGRV2=DUMGRV2+DUMC*DUMGRV[1]*2.0f;
      DUMC=C[7]+C[8];
      DUMGRV2=DUMGRV2+DUMC*DUMGRV[3];
      DUMGRV4=DUMGRV4+DUMC*DUMGRV[1];
      DUMC=C[9];
      DUMGRV4=DUMGRV4+DUMC*DUMGRV[3]*2.0f;

      // CONTRACTION OF FNC
      double FNCDOMA,FNCDOMAX,FNCDOMAY,FNCDOMAZ;
      FNCDOMA=DUMNV1;
      FNCDOMAX=DUMGRV1*ddx+DUMGRV3*ddx+
               DUMXX1*DXV[0]+DUMXX3*DXV[1];
      FNCDOMAY=DUMGRV1*ddy+DUMGRV3*ddy+
               DUMXX1*DYV[0]+DUMXX3*DYV[1];
      FNCDOMAZ=DUMGRV1*ddz+DUMGRV3*ddz+
               DUMXX1*DZV[0]+DUMXX3*DZV[1];

      double FNCDUMAX,FNCDUMAY,FNCDUMAZ;
      FNCDUMAX=DUMGRV1*DXV[0]+DUMGRV4*DXV[1];
      FNCDUMAY=DUMGRV1*DYV[0]+DUMGRV4*DYV[1];
      FNCDUMAZ=DUMGRV1*DZV[0]+DUMGRV4*DZV[1];

      // F CORE
      C[0]=2.0f*DUMXX[0];
      C[1]=DUMNV[0]*DUMNV[0];
      C[2]=2.0f*DUMNV[0]*DUMGRV[0];
      C[3]=2.0f*DUMGRV[0]*DUMNV[0];
      C[4]=DUMGRV[0]*DUMNV[0];
      C[5]=DUMNV[0]*DUMGRV[0];
      C[6]=4.0f*DUMGRV[0]*DUMGRV[0];
      C[7]=2.0f*DUMGRV[0]*DUMGRV[0];
      C[8]=2.0f*DUMGRV[0]*DUMGRV[0];
      C[9]=DUMGRV[0]*DUMGRV[0];

      // EXCHANGE
      double XDUMA,XDUMAG;
      XDUMA=C[0]*v2rhosigma[idx].x*4.0f;
      XDUMAG=C[0]*2.0f*v2sigma2[idx].x*8.0f;
      XDUMA=XDUMA+C[1]*v3rho3[idx].x*4.0f;
      XDUMAG=XDUMAG+C[1]*2.0f*v3rho2sigma[idx].x*8.0f;
      XDUMA=XDUMA+C[2]*v3rho2sigma[idx].x*8.0f;
      XDUMAG=XDUMAG+C[2]*2.0f*v3rhosigma2[idx].x*16.0f;
      XDUMA=XDUMA+C[3]*v3rho2sigma[idx].x*8.0f;
      XDUMAG=XDUMAG+C[3]*2.0f*v3rhosigma2[idx].x*16.0f;
      XDUMA=XDUMA+C[6]*v3rhosigma2[idx].x*16.0f;
      XDUMAG=XDUMAG+C[6]*2.0f*v3sigma3[idx].x*32.0f;

      // CORRELATION
      double CDUMA,CDUMAG1,CDUMAG2;
      CDUMA=C[0]*v2rhosigma[idx].y;
      CDUMAG1=C[0]*2.0f*v2sigma2[idx].y;
      CDUMAG2=C[0]*v2sigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[1]*v3rho3[idx].y;
      CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma[idx].y;
      CDUMAG2=CDUMAG2+C[1]*v3rho2sigma[idx].y*2.0f;
      CDUMA=CDUMA+C[2]*v3rho2sigma[idx].y;
      CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2[idx].y;
      CDUMAG2=CDUMAG2+C[2]*v3rhosigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[3]*v3rho2sigma[idx].y;
      CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2[idx].y;
      CDUMAG2=CDUMAG2+C[3]*v3rhosigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[4]*v3rho2sigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[4]*v3rhosigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[5]*v3rho2sigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[5]*v3rhosigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[6]*v3rhosigma2[idx].y;
      CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3[idx].y;
      CDUMAG2=CDUMAG2+C[6]*v3sigma3[idx].y*2.0f;
      CDUMA=CDUMA+C[7]*v3rhosigma2[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[7]*v3sigma3[idx].y*4.0f;
      CDUMA=CDUMA+C[8]*v3rhosigma2[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[8]*v3sigma3[idx].y*4.0f;
      CDUMA=CDUMA+C[9]*v3rhosigma2[idx].y*4.0f;
      CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3[idx].y*4.0f;
      CDUMAG2=CDUMAG2+C[9]*v3sigma3[idx].y*8.0f;

      C[0]=2.0f*DUMXX[1];
      C[1]=DUMNV[1]*DUMNV[1];
      C[2]=2.0f*DUMNV[1]*DUMGRV[1];
      C[3]=2.0f*DUMGRV[1]*DUMNV[1];
      C[4]=DUMGRV[1]*DUMNV[1];
      C[5]=DUMNV[1]*DUMGRV[1];
      C[6]=4.0f*DUMGRV[1]*DUMGRV[1];
      C[7]=2.0f*DUMGRV[1]*DUMGRV[1];
      C[8]=2.0f*DUMGRV[1]*DUMGRV[1];
      C[9]=DUMGRV[1]*DUMGRV[1];

      CDUMA=CDUMA+C[0]*v2rhosigma[idx].y;
      CDUMAG1=CDUMAG1+C[0]*2.0f*v2sigma2[idx].y;
      CDUMAG2=CDUMAG2+C[0]*v2sigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[1]*v3rho3[idx].y;
      CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma[idx].y;
      CDUMAG2=CDUMAG2+C[1]*v3rho2sigma[idx].y*2.0f;
      CDUMA=CDUMA+C[2]*v3rho2sigma[idx].y;
      CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2[idx].y;
      CDUMAG2=CDUMAG2+C[2]*v3rhosigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[3]*v3rho2sigma[idx].y;
      CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2[idx].y;
      CDUMAG2=CDUMAG2+C[3]*v3rhosigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[4]*v3rho2sigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[4]*v3rhosigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[5]*v3rho2sigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[5]*v3rhosigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[6]*v3rhosigma2[idx].y;
      CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3[idx].y;
      CDUMAG2=CDUMAG2+C[6]*v3sigma3[idx].y*2.0f;
      CDUMA=CDUMA+C[7]*v3rhosigma2[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[7]*v3sigma3[idx].y*4.0f;
      CDUMA=CDUMA+C[8]*v3rhosigma2[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[8]*v3sigma3[idx].y*4.0f;
      CDUMA=CDUMA+C[9]*v3rhosigma2[idx].y*4.0f;
      CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3[idx].y*4.0f;
      CDUMAG2=CDUMAG2+C[9]*v3sigma3[idx].y*8.0f;

      C[10]=DUMXX[2];
      C[11]=DUMNV[0]*DUMNV[1];
      C[12]=2.0f*DUMNV[0]*DUMGRV[1];
      C[13]=2.0f*DUMGRV[0]*DUMNV[1];
      C[14]=DUMNV[0]*DUMGRV[3];
      C[15]=DUMGRV[2]*DUMNV[1];
      C[16]=4.0f*DUMGRV[0]*DUMGRV[1];
      C[17]=2.0f*DUMGRV[0]*DUMGRV[3];
      C[18]=2.0f*DUMGRV[2]*DUMGRV[1];
      C[19]=DUMGRV[2]*DUMGRV[3];

      CDUMA=CDUMA+C[10]*v2rhosigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[10]*v2sigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[11]*v3rho3[idx].y;
      CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma[idx].y;
      CDUMAG2=CDUMAG2+C[11]*v3rho2sigma[idx].y*2.0f;
      CDUMA=CDUMA+C[12]*v3rho2sigma[idx].y;
      CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2[idx].y;
      CDUMAG2=CDUMAG2+C[12]*v3rhosigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[13]*v3rho2sigma[idx].y;
      CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2[idx].y;
      CDUMAG2=CDUMAG2+C[13]*v3rhosigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[14]*v3rho2sigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[14]*v3rhosigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[15]*v3rho2sigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[15]*v3rhosigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[16]*v3rhosigma2[idx].y;
      CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3[idx].y;
      CDUMAG2=CDUMAG2+C[16]*v3sigma3[idx].y*2.0f;
      CDUMA=CDUMA+C[17]*v3rhosigma2[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[17]*v3sigma3[idx].y*4.0f;
      CDUMA=CDUMA+C[18]*v3rhosigma2[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[18]*v3sigma3[idx].y*4.0f;
      CDUMA=CDUMA+C[19]*v3rhosigma2[idx].y*4.0f;
      CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3[idx].y*4.0f;
      CDUMAG2=CDUMAG2+C[19]*v3sigma3[idx].y*8.0f;

      // **BA_TERM=AB_TERM**
      C[10]=DUMXX[2];
      C[11]=DUMNV[0]*DUMNV[1];
      C[12]=2.0f*DUMNV[0]*DUMGRV[1];
      C[13]=2.0f*DUMGRV[0]*DUMNV[1];
      C[14]=DUMNV[0]*DUMGRV[3];
      C[15]=DUMGRV[2]*DUMNV[1];
      C[16]=4.0f*DUMGRV[0]*DUMGRV[1];
      C[17]=2.0f*DUMGRV[0]*DUMGRV[3];
      C[18]=2.0f*DUMGRV[2]*DUMGRV[1];
      C[19]=DUMGRV[2]*DUMGRV[3];

      CDUMA=CDUMA+C[10]*v2rhosigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[10]*v2sigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[11]*v3rho3[idx].y;
      CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma[idx].y;
      CDUMAG2=CDUMAG2+C[11]*v3rho2sigma[idx].y*2.0f;
      CDUMA=CDUMA+C[12]*v3rho2sigma[idx].y;
      CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2[idx].y;
      CDUMAG2=CDUMAG2+C[12]*v3rhosigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[13]*v3rho2sigma[idx].y;
      CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2[idx].y;
      CDUMAG2=CDUMAG2+C[13]*v3rhosigma2[idx].y*2.0f;
      CDUMA=CDUMA+C[14]*v3rho2sigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[14]*v3rhosigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[15]*v3rho2sigma[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[15]*v3rhosigma2[idx].y*4.0f;
      CDUMA=CDUMA+C[16]*v3rhosigma2[idx].y;
      CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3[idx].y;
      CDUMAG2=CDUMAG2+C[16]*v3sigma3[idx].y*2.0f;
      CDUMA=CDUMA+C[17]*v3rhosigma2[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[17]*v3sigma3[idx].y*4.0f;
      CDUMA=CDUMA+C[18]*v3rhosigma2[idx].y*2.0f;
      CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3[idx].y*2.0f;
      CDUMAG2=CDUMAG2+C[18]*v3sigma3[idx].y*4.0f;
      CDUMA=CDUMA+C[19]*v3rhosigma2[idx].y*4.0f;
      CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3[idx].y*4.0f;
      CDUMAG2=CDUMAG2+C[19]*v3sigma3[idx].y*8.0f;

      // CONTRACTION OF FC
      double FCDUMA,FCDUMAG1,FCDUMAG2,FCDUMAX,FCDUMAY,FCDUMAZ;
      FCDUMA=XDUMA+CDUMA;
      FCDUMAG1=XDUMAG+CDUMAG1;
      FCDUMAG2=CDUMAG2;
      FCDUMAX=FCDUMAG1*ddx+FCDUMAG2*ddx;
      FCDUMAY=FCDUMAG1*ddy+FCDUMAG2*ddy;
      FCDUMAZ=FCDUMAG1*ddz+FCDUMAG2*ddz;

      // DENSITY FACTOR
      DDUM[idx]+=FCDUMA;
      DDUMxyz[idx].x+=FNCDUMAX+FCDUMAX;
      DDUMxyz[idx].y+=FNCDUMAY+FCDUMAY;
      DDUMxyz[idx].z+=FNCDUMAZ+FCDUMAZ;

      // TRANSITION DENSITY FACTOR
      VDUM[idx]=FNCDOMA;
      VDUMxyz[idx].x=FNCDOMAX;
      VDUMxyz[idx].y=FNCDOMAY;
      VDUMxyz[idx].z=FNCDOMAZ;
   }
}
