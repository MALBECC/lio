#ifndef _ERI_H
#define _ERI_H
  #if FULL_DOUBLE
  void eri(double*,int,int,int,uint,uint*,double*,double*,double*,uint*,
           int,int,int);
  #else
  void eri(float*,int,int,int,uint,uint*,float*,float*,float*,uint*,
           int,int,int);
  #endif
#endif
