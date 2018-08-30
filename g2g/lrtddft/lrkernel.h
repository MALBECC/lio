#ifndef _LRKERNEL_H
#define _LRKERNEL_H
 #if FULL_DOUBLE
 void lrkernel(double*,double*,double*,double*,double*,double*,int&);
 #else
 void (lrkernel(float*,float*,float*,float*,float*,float*,int&);
 #endif
#endif
