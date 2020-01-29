#ifndef LIBXCPROXY_H
#define LIBXCPROXY_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <xc.h>
#include <vector>
#include "../scalar_vector_types.h"
#include "print_utils.h"
#include "../fix_compile.h"
//#include <cuda_runtime.h>
#include "../timer.h"

extern "C" void g2g_timer_sum_start_(const char* timer_name, unsigned int length_arg);
extern "C" void g2g_timer_sum_stop_(const char* timer_name, unsigned int length_arg);
extern "C" void g2g_timer_sum_pause_(const char* timer_name, unsigned int length_arg);

template <class T, int width>
class LibxcProxy
{
private:

    // The libxc components
    xc_func_type funcForExchange;
    xc_func_type funcForCorrelation;

    // Functional ids
    int funcIdForExchange;
    int funcIdForCorrelation;
    int nspin;

    // Is inited
    bool inited;
    double fact_exchange;
    void printFunctionalInformation (xc_func_type* func);

    // 3rd order derivative
    void Zv_exchange(double td, double* dxyz, double* txyz,
           double* Coef, double v2rho2, double v2rhosigma, double v2sigma2,
           double v3rho3,double v3rho2sigma,double v3rhosigma2,double v3sigma3,double fact_ex);

    void Zv_coulomb(double td,double* dxyz,double* txyz, double* Coef, 
                double v2rhosigma, double v2sigma2,double v3rho3,
                double v3rho2sigma,double v3rhosigma2,double v3sigma3);


public:
    LibxcProxy ();
    LibxcProxy (int exchangeFunctionId, int correlationFuncionalId, int nspin, double fexc);
    ~LibxcProxy ();

    void doGGA (T dens,
                const G2G::vec_type<T,width>& grad,
                const G2G::vec_type<T,width>& hess1,
                const G2G::vec_type<T,width>& hess2,
                T& ex,
                T& ec,
                T& y2a);

    void doGGA (T* dens,
                const int number_of_points,
		const G2G::vec_type<T,width>* grad,
                const G2G::vec_type<T,width>* hess1,
                const G2G::vec_type<T,width>* hess2,
                T* ex,
                T* ec,
                T* y2a);

    void doGGA (T dens,
                T sigma,
                T* v2rho2,
                T v2rhosigma,
                T v2sigma2);

    void doGGA (T* dens,
		T* sigma,
                const int number_of_points,
		T* v2rho2,
		T* v2rhosigma,
		T* v2sigma2);

    void doGGA (T* dens,
                const int number_of_points,
		const T* contracted_grad,
		const G2G::vec_type<T,width>* grad,
                const G2G::vec_type<T,width>* hess1,
                const G2G::vec_type<T,width>* hess2,
                T* ex,
                T* ec,
                T* y2a);

    void coefLR(double* rho, double* sigma, double red,
                double cruz, double* lrCoef);

    void coefZv(double pd, double sigma, double pdx, double pdy,
                double pdz, double red, double redx, double redy,
                double redz, double* zcoef);

    void terms_derivs(double* dens, double sigma,
                double* vrho, double* vsigma, double* v2rho2,
                double* v2rhosigma, double* v2sigma2, double* v3rho3,
                double* v3rho2sigma, double* v3rhosigma2, double* v3sigma3);

    void doLDA (T dens,
                const G2G::vec_type<T,width>& grad,
                const G2G::vec_type<T,width>& hess1,
                const G2G::vec_type<T,width>& hess2,
                T& ex,
                T& ec,
                T& y2a);

    void init (int exId, int xcId, int nspin);
    void closeProxy ();
    void printFunctionalsInformation (int exchangeFunctionalId, int correlationFunctionalId);
};

template <class T, int width>
LibxcProxy <T, width>::LibxcProxy()
{
    funcIdForExchange = 0;
    funcIdForCorrelation = 0;
    nspin = 0;
    inited = false;
}

template <class T, int width>
LibxcProxy <T, width>::LibxcProxy (int exchangeFunctionalId, int correlationFuncionalId, int nSpin, double fexc)
{
//    printf("LibxcProxy::LibxcProxy (%u, %u, %u) \n", exchangeFunctionalId, correlationFuncionalId, nSpin);
/*
    funcIdForExchange = exchangeFunctionalId;
    funcIdForCorrelation = correlationFuncionalId;
    nspin = nSpin;

    if (xc_func_init (&funcForExchange, funcIdForExchange, nspin) != 0) {
        fprintf (stderr, "Functional '%d' not found\n", funcIdForExchange);
	exit(-1);
    }

    if (xc_func_init (&funcForCorrelation, funcIdForCorrelation, nspin) != 0){
	fprintf (stderr, "Functional '%d' not found\n", funcIdForCorrelation);
	exit(-1);
    }
*/
    fact_exchange = fexc;
    init (exchangeFunctionalId, correlationFuncionalId, nSpin);
}

template <class T, int width>
void LibxcProxy<T, width>::init (int exchangeFunctionalId, int correlationFunctionalId, int nSpin)
{
//    printf("LibxcProxy::init(%u, %u, %u)\n", exchangeFunctionalId, correlationFunctionalId, nSpin);
    funcIdForExchange = exchangeFunctionalId;
    funcIdForCorrelation = correlationFunctionalId;
    nspin = nSpin;

    if (xc_func_init (&funcForExchange, funcIdForExchange, nspin) != 0) {
        fprintf (stderr, "Functional '%d' not found\n", funcIdForExchange);
	exit(-1);
    }

    if (xc_func_init (&funcForCorrelation, funcIdForCorrelation, nspin) != 0){
	fprintf (stderr, "Functional '%d' not found\n", funcIdForCorrelation);
	exit(-1);
    }

    inited = true;
}

template <class T, int width>
LibxcProxy <T, width>::~LibxcProxy ()
{
    //xc_func_end (&funcForExchange);
    //xc_func_end (&funcForCorrelation);
//    printf("LibxcProxy::~LibxcProxy()\n");
    closeProxy ();
}

template <class T, int width>
void LibxcProxy <T, width>::closeProxy ()
{
//    printf("LibxcProxy::closeProxy()\n");
    if (inited) {
	xc_func_end (&funcForExchange);
        xc_func_end (&funcForCorrelation);
	inited = false;
    }
}

template <class T, int width>
void LibxcProxy<T, width>::printFunctionalsInformation (int exchangeFunctionalId, int correlationFunctionalId) 
{
    if (!inited) 
    {
	init (exchangeFunctionalId, correlationFunctionalId, 1);
    }

    printFunctionalInformation (&funcForExchange);
    printFunctionalInformation (&funcForCorrelation);

    if (inited) 
    {
	closeProxy ();
    }
}

template <class T, int width>
void LibxcProxy<T, width>::printFunctionalInformation (xc_func_type* func)
{
    printf("The functional '%s' is ", func->info->name);
    switch (func->info->kind) {
	case (XC_EXCHANGE):
	    printf("an exchange functional");
	break;
	case (XC_CORRELATION):
	    printf("a correlation functional");
	break;
	case (XC_EXCHANGE_CORRELATION):
	    printf("an exchange-correlation functional");
	break;
	case (XC_KINETIC):
	    printf("a kinetic energy functional");
	break;
	default:
	    printf("of unknown kind");
	break;
    }

    printf(", it belongs to the '", func->info->name);
    switch (func->info->family) {
	case (XC_FAMILY_LDA):
	    printf("LDA");
        break;
	case (XC_FAMILY_GGA):
	    printf("GGA");
        break;
	case (XC_FAMILY_HYB_GGA):
	    printf("Hybrid GGA");
	break;
	case (XC_FAMILY_MGGA):
	    printf("MGGA");
        break;
	case (XC_FAMILY_HYB_MGGA):
	    printf("Hybrid MGGA");
        break;
	default:
	    printf("unknown");
        break;
    }
    printf("' family and is defined in the reference(s):\n");

    for (int ii = 0; func->info->refs[ii] != NULL; ii++) {
	printf ("[%d] %s\n", ii+1, func->info->refs[ii]->ref);
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doGGA - CPU Version 1
// Calls the XC_GGA function from libxc for multiple points.
// dens: pointer for the density array
// grad:.
// hess1:
// hess2:
// ex: here goes the results after calling xc_gga from libxc for the exchange functional
// ec: here goes the results after calling xc_gga from libxc for the correlation functional
// y2a:
//
template <class T, int width>
void LibxcProxy <T, width>::doGGA(T dens,
    const G2G::vec_type<T, width> &grad,
    const G2G::vec_type<T, width> &hess1,
    const G2G::vec_type<T, width> &hess2,
    T &ex, T &ec, T &y2a)
{
    //printf("LibxcProxy::doGGA cpu simple(...) \n");

    const double rho[1] = {dens};
    // Libxc needs the 'contracted gradient'
    double sigma[1] = {(grad.x * grad.x) + (grad.y * grad.y) + (grad.z * grad.z)};
    double exchange[1];
    double correlation[1];

    // The outputs for exchange
    double vrho [1];
    double vsigma [1];
    double v2rho [1];
    double v2rhosigma[1];
    double v2sigma [1];

    // The outputs for correlation
    double vrhoC [1];
    double vsigmaC [1];
    double v2rhoC [1];
    double v2rhosigmaC [1];
    double v2sigmaC [1];

    // The exchange values
    xc_gga (&funcForExchange, 1,
                rho,
                sigma,
                exchange,
                vrho,
                vsigma,
                v2rho,
                v2rhosigma,
                v2sigma,
                NULL, NULL, NULL, NULL);

    // Now the correlation value.
    xc_gga (&funcForCorrelation, 1,
                rho,
                sigma,
                correlation,
                vrhoC,
                vsigmaC,
                v2rhoC,
                v2rhosigmaC,
                v2sigmaC,
                NULL, NULL, NULL, NULL);

    ex = fact_exchange * exchange[0];
    ec = correlation[0];

    // Merge the results for the derivatives.
    vrho[0]       = vrho[0] * fact_exchange + vrhoC[0];
    vsigma[0]     = vsigma[0] * fact_exchange + vsigmaC[0];
    v2rho[0]      = v2rho[0] * fact_exchange + v2rhoC[0];
    v2rhosigma[0] = v2rhosigma[0] * fact_exchange + v2rhosigmaC[0];
    v2sigma[0]    = v2sigma[0] * fact_exchange + v2sigmaC[0];

    // Now, compute y2a value.
    y2a = vrho[0] - (2 * sigma[0] * v2rhosigma[0]
            + 2 * (hess1.x + hess1.y + hess1.z) * vsigma[0]
            + 4 * v2sigma[0] * (grad.x * grad.x * hess1.x + grad.y * grad.y * hess1.y + grad.z * grad.z * hess1.z + 2 * grad.x * grad.y * hess2.x + 2 * grad.x * grad.z * hess2.y + 2 * grad.y * grad.z * hess2.z));

    return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doGGA - CPU Version 2
// Calls the XC_GGA function from libxc for multiple points.
// dens: pointer for the density array
// number_of_points: the size of all the input arrays
// grad: gradient value in each point
// hess1:
// hess2:
// ex: here goes the results after calling xc_gga from libxc for the exchange functional
// ec: here goes the results after calling xc_gga from libxc for the correlation functional
// y2a:
//
template <class T, int width>
void LibxcProxy <T, width>::doGGA(T* dens,
    const int number_of_points,
    const G2G::vec_type<T, width>* grad,
    const G2G::vec_type<T, width>* hess1,
    const G2G::vec_type<T, width>* hess2,
    T* ex,
    T* ec,
    T* y2a)
{
    //printf("LibxcProxy::doGGA cpu multiple (...) \n");

    int array_size = sizeof(double)*number_of_points;
    double* rho = (double*)malloc(array_size);
    for (int i=0; i<number_of_points; i++) {
	rho[i] = (double)dens[i];
    }

    // Libxc needs the 'contracted gradient'
    double* sigma = (double*)malloc(array_size);
    for (int i=0; i< number_of_points; i++) {
	sigma[i] = (double)((grad[i].x * grad[i].x) + (grad[i].y * grad[i].y) + (grad[i].z * grad[i].z));
    }
    double* exchange = (double*)malloc(array_size);
    double* correlation = (double*)malloc(array_size);

    // The outputs for exchange
    double* vrho = (double*)malloc(array_size);
    double* vsigma = (double*)malloc(array_size);
    double* v2rho = (double*)malloc(array_size);
    double* v2rhosigma = (double*)malloc(array_size);
    double* v2sigma = (double*)malloc(array_size);

    // The outputs for correlation
    double* vrhoC = (double*)malloc(array_size);
    double* vsigmaC = (double*)malloc(array_size);
    double* v2rhoC = (double*)malloc(array_size);
    double* v2rhosigmaC = (double*)malloc(array_size);
    double* v2sigmaC = (double*)malloc(array_size);

    // Exchange values
    xc_gga (&funcForExchange, number_of_points,
                rho,
                sigma,
                exchange,
                vrho,
                vsigma,
                v2rho,
                v2rhosigma,
                v2sigma,
                NULL, NULL, NULL, NULL);

    // Now the correlation value.
    xc_gga (&funcForCorrelation, number_of_points,
                rho,
                sigma,
                correlation,
                vrhoC,
                vsigmaC,
                v2rhoC,
                v2rhosigmaC,
                v2sigmaC,
                NULL, NULL, NULL, NULL);

    for (int i=0; i<number_of_points; i++) {
	ex[i] = exchange[i];
        ec[i] = correlation[i];
        // Merge the results for the derivatives.
        vrho[i] += vrhoC[i];
        vsigma[i] += vsigmaC[i];
        v2rho[i] += v2rhoC[i];
        v2rhosigma[i] += v2rhosigmaC[i];
        v2sigma[i] += v2sigmaC[i];

	// Now, compute y2a value.
        y2a[i] = vrho[i] - (2 * sigma[i] * v2rhosigma[i]
	        + 2 * (hess1[i].x + hess1[i].y + hess1[i].z) * vsigma[i]
    		+ 4 * v2sigma[i] * (grad[i].x * grad[i].x * hess1[i].x +
				    grad[i].y * grad[i].y * hess1[i].y +
				    grad[i].z * grad[i].z * hess1[i].z +
				    2 * grad[i].x * grad[i].y * hess2[i].x +
				    2 * grad[i].x * grad[i].z * hess2[i].y +
				    2 * grad[i].y * grad[i].z * hess2[i].z));
    }

    // Free memory.
    free(rho);
    free(sigma);
    free(exchange);
    free(correlation);

    // The outputs for exchange
    free(vrho);
    free(vsigma);
    free(v2rho);
    free(v2rhosigma);
    free(v2sigma);

    // The outputs for correlation
    free(vrhoC);
    free(vsigmaC);
    free(v2rhoC);
    free(v2rhosigmaC);
    free(v2sigmaC);

    return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doGGA - CPU Version 3
// Calls the XC_GGA_FXC function from libxc for multiple points.
// rho: pointer for the density array.
// sigma: contracted gradient array.
// number_of_points: the size of all the input arrays.
// v2rho2: second partial derivative of the energy per unit volume in terms of the density.
// v2rhosigma: second partial derivative of the energy per unit volume in terms of the density and sigma.
// v2sigma2: second partial derivative of the energy per unit volume in terms of sigma.
//
template <class T, int width>
void LibxcProxy <T, width>::doGGA (T* rho,
    T* sigma,
    const int number_of_points,
    T* v2rho2,
    T* v2rhosigma,
    T* v2sigma2)
{
    int array_size = sizeof(double)*number_of_points;

    // The outputs for exchange
    double* v2rho2X = (double*)malloc(array_size);
    double* v2rhosigmaX = (double*)malloc(array_size);
    double* v2sigma2X = (double*)malloc(array_size);

    // The outputs for correlation
    double* v2rho2C = (double*)malloc(array_size);
    double* v2rhosigmaC = (double*)malloc(array_size);
    double* v2sigma2C = (double*)malloc(array_size);

    // Exchange values
    xc_gga_fxc (&funcForExchange, number_of_points,
                rho,
                sigma,
                v2rho2X,
                v2rhosigmaX,
                v2sigma2X);

    // Now the correlation value.
    xc_gga_fxc (&funcForCorrelation, number_of_points,
                rho,
                sigma,
                v2rho2C,
                v2rhosigmaC,
                v2sigma2C);

    for (int i=0; i<number_of_points; i++) {
        // Merge the results for the derivatives.
        v2rho2[i] = v2rho2X[i] + v2rho2C[i];
        v2rhosigma[i] = v2rhosigmaX[i] + v2rhosigmaC[i];
        v2sigma2[i] += v2sigma2X[i] + v2sigma2C[i];
    }

    // Free memory
    // The outputs for exchange
    free(v2rho2X);
    free(v2rhosigmaX);
    free(v2sigma2X);

    // The outputs for correlation
    free(v2rho2C);
    free(v2rhosigmaC);
    free(v2sigma2C);

    return;
}

template <class T, int width>
void LibxcProxy <T, width>::doGGA (T rho,
               T sigma,
               T* v2rho2,
               T v2rhosigma,
               T v2sigma2)
{
    const double dens = rho;
    const double sig = sigma;

    // The outputs for exchange
    double v2rho2X = 0.0;
    double v2rhosigmaX = 0.0;
    double v2sigma2X = 0.0;

    // The outputs for correlation
    double v2rho2C = 0.0;
    double v2rhosigmaC = 0.0;
    double v2sigma2C = 0.0;

    // Exchange values
    xc_gga_fxc (&funcForExchange,1,
                &dens,
                &sig,
                &v2rho2X,
                &v2rhosigmaX,
                &v2sigma2X);

    // Correlation values
    xc_gga_fxc (&funcForCorrelation,1,
                &dens,
                &sig,
                &v2rho2C,
                &v2rhosigmaC,
                &v2sigma2C);


    *v2rho2 = v2rho2C + v2rho2X;
    return;
}

#ifdef __CUDACC__

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::joinResults
// Kernel for join the partial results that we obtain from
// calling Libxc.
//
template <class T, int width>
__global__ void joinResults(
		    double* ex, double* exchange,
		    double* ec, double* correlation,
		    double* vrho, double* vrhoC,
		    double* vsigma, double* vsigmaC,
		    double* v2rho, double* v2rhoC,
		    double* v2rhosigma, double* v2rhosigmaC,
		    double* v2sigma, double* v2sigmaC,
		    double* y2a,
		    const double* sigma,
		    const G2G::vec_type<double, width>* grad,
		    const G2G::vec_type<double, width>* hess1,
		    const G2G::vec_type<double, width>* hess2,
		    int numElements,double fEE)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
	ex[i] = fEE * exchange[i];
	ec[i] = correlation[i];

	// Merge the results for the derivatives.
	vrho[i]       = fEE * vrho[i] + vrhoC[i];
        vsigma[i]     = fEE * v2sigma[i] + vsigmaC[i];
        v2rho[i]      = fEE * v2rho[i] + v2rhoC[i];
        v2rhosigma[i] = fEE * v2rhosigma[i] + v2rhosigmaC[i];
        v2sigma[i]    = fEE * v2sigma[i] + v2sigmaC[i];
        // Now, compute y2a value.
	y2a[i] = vrho[i] - (2 * sigma[i] * v2rhosigma[i]
            + 2 * (hess1[i].x + hess1[i].y + hess1[i].z) * vsigma[i]
            + 4 * v2sigma[i] * (grad[i].x * grad[i].x * hess1[i].x +
					    grad[i].y * grad[i].y * hess1[i].y +
					    grad[i].z * grad[i].z * hess1[i].z +
					    2 * grad[i].x * grad[i].y * hess2[i].x +
					    2 * grad[i].x * grad[i].z * hess2[i].y +
					    2 * grad[i].y * grad[i].z * hess2[i].z));
    }
}

/////////////////////////////////////
// Conversion KERNELS
//
// Utils for data type conversion from lio to libxc
__global__ void convertFloatToDouble(const float* input, double* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	output[i] = (double)input[i];
    }
}

__global__ void convertDoubleToFloat(const double* input, float* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	output[i] = (float)input[i];
    }
}

__global__ void convertFloatToDouble(const G2G::vec_type<float,4>* input, G2G::vec_type<double,4>* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	//float x, y, z, _w;
	output[i].x = (double)(input[i].x);
	output[i].y = (double)(input[i].y);
	output[i].z = (double)(input[i].z);
	output[i].w = (double)(input[i].w);
    }
}

// Esto es para engañar al compilador pq el FLAG FULL_DOUBLE
// a veces permite que T sea double y se rompe todo.
__global__ void convertDoubleToFloat(const double* input, double* output, int numElements)
{
    return;
}

// Esto es para engañar al compilador pq el FLAG FULL_DOUBLE
// a veces permite que T sea double y se rompe todo.
__global__ void convertFloatToDouble(const double* input, double* output, int numElements)
{
    return;
}

// Esto es para engañar al compilador pq el FLAG FULL_DOUBLE
// a veces permite que T sea double y se rompe todo.
__global__ void convertFloatToDouble(const G2G::vec_type<double,4>* input, G2G::vec_type<double,4>* output, int numElements)
{
    return;
}

__global__ void convertDoubleToFloat(const G2G::vec_type<double,4>* input, G2G::vec_type<float,4>* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	//float x, y, z, _w;
	output[i].x = (float)(input[i].x);
	output[i].y = (float)(input[i].y);
	output[i].z = (float)(input[i].z);
	output[i].w = (float)(input[i].w);
    }
}

// end Conversion KERNELS
////////////////////////////////

#endif


template <class T, int width>
void LibxcProxy <T, width>::coefLR (double *rho,
                double* sgm,
                double red,
                double cruz,
                double* lrCoef)
{
   // The otputs for exchange
   double vrhoX, vsigmaX, v2rho2X, v2rhosigmaX, v2sigma2X;

   // The ouputs for correlation
   double vrhoC, vsigmaC, v2rho2C, v2rhosigmaC, v2sigma2C;

   // NOT Refence
   double exc = 0.0f;

   // convert to alfa and beta;
   double dens, sigma; dens = *rho; sigma = *sgm;

   // Exchange values
   xc_gga(&funcForExchange,1,&dens,&sigma,&exc,&vrhoX,&vsigmaX,
          &v2rho2X,&v2rhosigmaX,&v2sigma2X,NULL,NULL,NULL,NULL);

   // Correlation values
   xc_gga(&funcForCorrelation,1,&dens,&sigma,&exc,&vrhoC,&vsigmaC,
          &v2rho2C,&v2rhosigmaC,&v2sigma2C,NULL,NULL,NULL,NULL);

   // Results
   lrCoef[0] = (2.0f * red * v2rho2X + 8.0f * cruz * v2rhosigmaX) * fact_exchange;
   lrCoef[1] = (8.0f * red * v2rhosigmaX + 32.0f * cruz * v2sigma2X) * fact_exchange;
   lrCoef[2] = 4.0f * vsigmaX * fact_exchange;

   lrCoef[0] += 2.0f * red * v2rho2C + 8.0f * cruz * v2rhosigmaC;
   lrCoef[1] += 8.0f * red * v2rhosigmaC + 32.0f * cruz * v2sigma2C;
   lrCoef[2] += 4.0f * vsigmaC;

   return;
}

template <class T, int width>
void LibxcProxy<T, width>::Zv_exchange(double td, double* dxyz, double* txyz,
           double* Coef, double v2rho2, double v2rhosigma, double v2sigma2,
           double v3rho3,double v3rho2sigma,double v3rhosigma2,double v3sigma3,double fact_ex)
{

   double fex = fact_ex;
   double DUMNV[2],DXV[2],DYV[2],DZV[2],DUMGRV[4],DUMXX[4];
   double C[10];

   DUMNV[0] = DUMNV[1] = td;

   DUMGRV[0]=DUMGRV[1]=txyz[0]*dxyz[0]*0.5f + txyz[1]*dxyz[1]*0.5f + txyz[2]*dxyz[2]*0.5f;
   DUMGRV[2]=DUMGRV[3]=DUMGRV[0];

   DUMXX[0]=DUMXX[1]=txyz[0]*txyz[0] + txyz[1]*txyz[1] + txyz[2]*txyz[2];
   DUMXX[2]=DUMXX[3]=DUMXX[0];

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

   double XDUMA=0.0f;
   double XDUMAG=0.0f;
   XDUMA  = C[0] * 4.00f * v2rhosigma;
   XDUMAG = C[0] * 16.0f * v2sigma2;
   XDUMA += C[1] * 4.00f * v3rho3;
   XDUMAG+= C[1] * 16.0f * v3rho2sigma;
   XDUMA += C[2] * 8.00f * v3rho2sigma;
   XDUMAG+= C[2] * 32.0f * v3rhosigma2;
   XDUMA += C[3] * 8.00f * v3rho2sigma;
   XDUMAG+= C[3] * 32.0f * v3rhosigma2;
   XDUMA += C[6] * 16.0f * v3rhosigma2;
   XDUMAG+= C[6] * 64.0f * v3sigma3;

   double XDUMAGEA=0.0f;
   XDUMAGEA  = 4.0f * DUMNV[0]  * 4.0f * v2rhosigma;
   XDUMAGEA += 8.0f * DUMGRV[1] * 8.0f * v2sigma2;

   Coef[0] = XDUMA * fex;
   Coef[1] = XDUMAG * fex;
   Coef[2] = XDUMAGEA * fex;

}

template <class T, int width>
void LibxcProxy<T, width>::Zv_coulomb(double td,double* dxyz,double* txyz, double* Coef, 
                double v2rhosigma, double v2sigma2,double v3rho3,
                double v3rho2sigma,double v3rhosigma2,double v3sigma3)
{
   double DUMNV[2],DXV[2],DYV[2],DZV[2],DUMGRV[4],DUMXX[4];
   double C[20];

   DUMNV[0] = DUMNV[1] = td;

   DUMGRV[0]=DUMGRV[1]=txyz[0]*dxyz[0]*0.5f + txyz[1]*dxyz[1]*0.5f + txyz[2]*dxyz[2]*0.5f;
   DUMGRV[2]=DUMGRV[3]=DUMGRV[0];

   DUMXX[0]=DUMXX[1]=txyz[0]*txyz[0] + txyz[1]*txyz[1] + txyz[2]*txyz[2];
   DUMXX[2]=DUMXX[3]=DUMXX[0];

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

   double CDUMA=0.0f;
   double CDUMAG1=0.0f;
   double CDUMAG2=0.0f;
   CDUMA=C[0]*v2rhosigma;
   CDUMAG1=C[0]*2.0f*v2sigma2;
   CDUMAG2=C[0]*v2sigma2*2.0f;
   CDUMA=CDUMA+C[1]*v3rho3;
   CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma;
   CDUMAG2=CDUMAG2+C[1]*v3rho2sigma*2.0f;
   CDUMA=CDUMA+C[2]*v3rho2sigma;
   CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2;
   CDUMAG2=CDUMAG2+C[2]*v3rhosigma2*2.0f;
   CDUMA=CDUMA+C[3]*v3rho2sigma;
   CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2;
   CDUMAG2=CDUMAG2+C[3]*v3rhosigma2*2.0f;
   CDUMA=CDUMA+C[4]*v3rho2sigma*2.0f;
   CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2*2.0f;
   CDUMAG2=CDUMAG2+C[4]*v3rhosigma2*4.0f;
   CDUMA=CDUMA+C[5]*v3rho2sigma*2.0f;
   CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2*2.0f;
   CDUMAG2=CDUMAG2+C[5]*v3rhosigma2*4.0f;
   CDUMA=CDUMA+C[6]*v3rhosigma2;
   CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3;
   CDUMAG2=CDUMAG2+C[6]*v3sigma3*2.0f;
   CDUMA=CDUMA+C[7]*v3rhosigma2*2.0f;
   CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3*2.0f;
   CDUMAG2=CDUMAG2+C[7]*v3sigma3*4.0f;
   CDUMA=CDUMA+C[8]*v3rhosigma2*2.0f;
   CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3*2.0f;
   CDUMAG2=CDUMAG2+C[8]*v3sigma3*4.0f;
   CDUMA=CDUMA+C[9]*v3rhosigma2*4.0f;
   CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3*4.0f;
   CDUMAG2=CDUMAG2+C[9]*v3sigma3*8.0f;

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

   CDUMA=CDUMA+C[0]*v2rhosigma;
   CDUMAG1=CDUMAG1+C[0]*2.0f*v2sigma2;
   CDUMAG2=CDUMAG2+C[0]*v2sigma2*2.0f;
   CDUMA=CDUMA+C[1]*v3rho3;
   CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma;
   CDUMAG2=CDUMAG2+C[1]*v3rho2sigma*2.0f;
   CDUMA=CDUMA+C[2]*v3rho2sigma;
   CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2;
   CDUMAG2=CDUMAG2+C[2]*v3rhosigma2*2.0f;
   CDUMA=CDUMA+C[3]*v3rho2sigma;
   CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2;
   CDUMAG2=CDUMAG2+C[3]*v3rhosigma2*2.0f;
   CDUMA=CDUMA+C[4]*v3rho2sigma*2.0f;
   CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2*2.0f;
   CDUMAG2=CDUMAG2+C[4]*v3rhosigma2*4.0f;
   CDUMA=CDUMA+C[5]*v3rho2sigma*2.0f;
   CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2*2.0f;
   CDUMAG2=CDUMAG2+C[5]*v3rhosigma2*4.0f;
   CDUMA=CDUMA+C[6]*v3rhosigma2;
   CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3;
   CDUMAG2=CDUMAG2+C[6]*v3sigma3*2.0f;
   CDUMA=CDUMA+C[7]*v3rhosigma2*2.0f;
   CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3*2.0f;
   CDUMAG2=CDUMAG2+C[7]*v3sigma3*4.0f;
   CDUMA=CDUMA+C[8]*v3rhosigma2*2.0f;
   CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3*2.0f;
   CDUMAG2=CDUMAG2+C[8]*v3sigma3*4.0f;
   CDUMA=CDUMA+C[9]*v3rhosigma2*4.0f;
   CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3*4.0f;
   CDUMAG2=CDUMAG2+C[9]*v3sigma3*8.0f;

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

   CDUMA=CDUMA+C[10]*v2rhosigma*2.0f;
   CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2*2.0f;
   CDUMAG2=CDUMAG2+C[10]*v2sigma2*4.0f;
   CDUMA=CDUMA+C[11]*v3rho3;
   CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma;
   CDUMAG2=CDUMAG2+C[11]*v3rho2sigma*2.0f;
   CDUMA=CDUMA+C[12]*v3rho2sigma;
   CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2;
   CDUMAG2=CDUMAG2+C[12]*v3rhosigma2*2.0f;
   CDUMA=CDUMA+C[13]*v3rho2sigma;
   CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2;
   CDUMAG2=CDUMAG2+C[13]*v3rhosigma2*2.0f;
   CDUMA=CDUMA+C[14]*v3rho2sigma*2.0f;
   CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2*2.0f;
   CDUMAG2=CDUMAG2+C[14]*v3rhosigma2*4.0f;
   CDUMA=CDUMA+C[15]*v3rho2sigma*2.0f;
   CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2*2.0f;
   CDUMAG2=CDUMAG2+C[15]*v3rhosigma2*4.0f;
   CDUMA=CDUMA+C[16]*v3rhosigma2;
   CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3;
   CDUMAG2=CDUMAG2+C[16]*v3sigma3*2.0f;
   CDUMA=CDUMA+C[17]*v3rhosigma2*2.0f;
   CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3*2.0f;
   CDUMAG2=CDUMAG2+C[17]*v3sigma3*4.0f;
   CDUMA=CDUMA+C[18]*v3rhosigma2*2.0f;
   CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3*2.0f;
   CDUMAG2=CDUMAG2+C[18]*v3sigma3*4.0f;
   CDUMA=CDUMA+C[19]*v3rhosigma2*4.0f;
   CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3*4.0f;
   CDUMAG2=CDUMAG2+C[19]*v3sigma3*8.0f;

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

   CDUMA=CDUMA+C[10]*v2rhosigma*2.0f;
   CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2*2.0f;
   CDUMAG2=CDUMAG2+C[10]*v2sigma2*4.0f;
   CDUMA=CDUMA+C[11]*v3rho3;
   CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma;
   CDUMAG2=CDUMAG2+C[11]*v3rho2sigma*2.0f;
   CDUMA=CDUMA+C[12]*v3rho2sigma;
   CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2;
   CDUMAG2=CDUMAG2+C[12]*v3rhosigma2*2.0f;
   CDUMA=CDUMA+C[13]*v3rho2sigma;
   CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2;
   CDUMAG2=CDUMAG2+C[13]*v3rhosigma2*2.0f;
   CDUMA=CDUMA+C[14]*v3rho2sigma*2.0f;
   CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2*2.0f;
   CDUMAG2=CDUMAG2+C[14]*v3rhosigma2*4.0f;
   CDUMA=CDUMA+C[15]*v3rho2sigma*2.0f;
   CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2*2.0f;
   CDUMAG2=CDUMAG2+C[15]*v3rhosigma2*4.0f;
   CDUMA=CDUMA+C[16]*v3rhosigma2;
   CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3;
   CDUMAG2=CDUMAG2+C[16]*v3sigma3*2.0f;
   CDUMA=CDUMA+C[17]*v3rhosigma2*2.0f;
   CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3*2.0f;
   CDUMAG2=CDUMAG2+C[17]*v3sigma3*4.0f;
   CDUMA=CDUMA+C[18]*v3rhosigma2*2.0f;
   CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3*2.0f;
   CDUMAG2=CDUMAG2+C[18]*v3sigma3*4.0f;
   CDUMA=CDUMA+C[19]*v3rhosigma2*4.0f;
   CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3*4.0f;
   CDUMAG2=CDUMAG2+C[19]*v3sigma3*8.0f;

   double CDUMAGEA=0.0f;
   CDUMAGEA=CDUMAGEA+4.0f*DUMNV[0]*v2rhosigma;
   CDUMAGEA=CDUMAGEA+8.0f*DUMGRV[0]*v2sigma2;
   CDUMAGEA=CDUMAGEA+4.0f*DUMGRV[2]*v2sigma2*2.0f;
   CDUMAGEA=CDUMAGEA+4.0f*DUMNV[1]*v2rhosigma;
   CDUMAGEA=CDUMAGEA+8.0f*DUMGRV[1]*v2sigma2;
   CDUMAGEA=CDUMAGEA+4.0f*DUMGRV[3]*v2sigma2*2.0f;

   double CDUMAGEC=0.0f;
   CDUMAGEC=CDUMAGEC+2.0f*DUMNV[0]*v2rhosigma*2.0f;
   CDUMAGEC=CDUMAGEC+4.0f*DUMGRV[0]*v2sigma2*2.0f;
   CDUMAGEC=CDUMAGEC+2.0f*DUMGRV[2]*v2sigma2*4.0f;
   CDUMAGEC=CDUMAGEC+2.0f*DUMNV[1]*v2rhosigma*2.0f;
   CDUMAGEC=CDUMAGEC+4.0f*DUMGRV[1]*v2sigma2*2.0f;
   CDUMAGEC=CDUMAGEC+2.0f*DUMGRV[3]*v2sigma2*4.0f;

   Coef[0] += CDUMA;
   Coef[1] += CDUMAG1+CDUMAG2;
   Coef[2] += CDUMAGEA+CDUMAGEC;

}

template <class T, int width>
void LibxcProxy<T,width>::terms_derivs(double* dens, double sigma,
               double* vrho, double* vsigma, double* v2rho2,
               double* v2rhosigma, double* v2sigma2, double* v3rho3,
               double* v3rho2sigma, double* v3rhosigma2, double* v3sigma3)
{
   double pd = dens[0];
   double fex = fact_exchange;
   // The otputs libxc
   double lvrho, lvsigma, lv2rho2, lv2rhosigma, lv2sigma2;
   double lv3rho3, lv3rho2sigma, lv3rhosigma2, lv3sigma3, exc;
   lvrho = lvsigma = 0.0f; lv2rho2 = lv2rhosigma = lv2sigma2 = 0.0f;
   lv3rho3 = lv3rho2sigma = lv3rhosigma2 = lv3sigma3 = exc = 0.0f;

   // Exchange Values
   xc_gga(&funcForExchange,1,&pd,&sigma,&exc,&lvrho,&lvsigma,
          &lv2rho2,&lv2rhosigma,&lv2sigma2,&lv3rho3,&lv3rho2sigma,
          &lv3rhosigma2,&lv3sigma3);

   // Copy Exchange Values
   vrho[0] = lvrho * fex; vsigma[0] = lvsigma * fex;
   v2rho2[0] = lv2rho2 * fex; v2rhosigma[0] = lv2rhosigma * fex;
   v2sigma2[0] = lv2sigma2 * fex; v3rho3[0] = lv3rho3 * fex; 
   v3rho2sigma[0] = lv3rho2sigma * fex; v3sigma3[0] = lv3sigma3 * fex;
   v3rhosigma2[0] = lv3rhosigma2 * fex;

   lvrho = lvsigma = 0.0f; lv2rho2 = lv2rhosigma = lv2sigma2 = 0.0f;
   lv3rho3 = lv3rho2sigma = lv3rhosigma2 = lv3sigma3 = exc = 0.0f;

   // Correlation Values
   xc_gga(&funcForCorrelation,1,&pd,&sigma,&exc,&lvrho,&lvsigma,
          &lv2rho2,&lv2rhosigma,&lv2sigma2,&lv3rho3,&lv3rho2sigma,
          &lv3rhosigma2,&lv3sigma3);

   // Copy Correlation Values
   vrho[1] = lvrho; vsigma[1] = lvsigma;
   v2rho2[1] = lv2rho2; v2rhosigma[1] = lv2rhosigma;
   v2sigma2[1] = lv2sigma2; v3rho3[1] = lv3rho3; 
   v3rho2sigma[1] = lv3rho2sigma; v3sigma3[1] = lv3sigma3;
   v3rhosigma2[1] = lv3rhosigma2;
}

template <class T, int width>
void LibxcProxy <T, width>::coefZv(double dens, double sigma, double pdx, double pdy,
                double pdz, double red, double redx, double redy, double redz,
                double* zcoef)
{
   double* dxyz = (double*)malloc(3*sizeof(double));
   double* txyz = (double*)malloc(3*sizeof(double));

   dxyz[0] = pdx; dxyz[1] = pdy; dxyz[2] = pdz;
   txyz[0] = redx; txyz[1] = redy; txyz[2] = redz;

   // The otputs libxc
   double vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma;
   double v3rhosigma2, v3sigma3, exc;
   vrho = vsigma = 0.0f; v2rho2 = v2rhosigma = v2sigma2 = 0.0f;
   v3rho3 = v3rho2sigma = v3rhosigma2 = v3sigma3 = exc = 0.0f;
   

   // Exchange Values
   xc_gga(&funcForExchange,1,&dens,&sigma,&exc,&vrho,&vsigma,
          &v2rho2,&v2rhosigma,&v2sigma2,&v3rho3,&v3rho2sigma,
          &v3rhosigma2,&v3sigma3);

   // Obtain Z coef of exchange
   Zv_exchange(red,dxyz,txyz,zcoef,v2rho2,v2rhosigma,v2sigma2,
               v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3,fact_exchange);

   vrho = vsigma = 0.0f; v2rho2 = v2rhosigma = v2sigma2 = 0.0f;
   v3rho3 = v3rho2sigma = v3rhosigma2 = v3sigma3 = exc = 0.0f;

   // Correlation Values
   xc_gga(&funcForCorrelation,1,&dens,&sigma,&exc,&vrho,&vsigma,
          &v2rho2,&v2rhosigma,&v2sigma2,&v3rho3,&v3rho2sigma,
          &v3rhosigma2,&v3sigma3);

   // Obtain Z coef of correlation
   Zv_coulomb(red,dxyz,txyz,zcoef,v2rhosigma,v2sigma2,
              v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3);

   return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doGGA - GPU version
// Calls the xc_gga function from libxc
// dens: pointer for the density array
// number_of_points: the size of all the input arrays
// contracted_grad: the contracted grad for libxc
// grad: 
// hess1:
// hess2:
// ex: here goes the results after calling xc_gga from libxc for the exchange functional
// ec: here goes the results after calling xc_gga from libxc for the correlation functional
//
// Note: all the pointer data are pointers in CUDA memory.
//
template <class T, int width>
void LibxcProxy <T, width>::doGGA(T* dens,
    const int number_of_points,
    const T* contracted_grad,
    const G2G::vec_type<T, width>* grad,
    const G2G::vec_type<T, width>* hess1,
    const G2G::vec_type<T, width>* hess2,
    T* ex,
    T* ec,
    T* y2a)
{
#ifdef __CUDACC__
    //printf("doGGA - GPU \n");
    //printf("Number of points: %u\n", number_of_points);

    // Este flag esta asi ya que a veces lio utiliza precision mixta
    // y solo en tiempo de ejecucion podemos saber que tipos
    // de datos esta utilizando.
    bool full_double = (sizeof(T) == 8);

    // Variables for the Kernels
    int threadsPerBlock = 256;
    int blocksPerGrid = (number_of_points + threadsPerBlock - 1) / threadsPerBlock;

    cudaError_t err = cudaSuccess;

    // All the arrays for libxc must be of double*
    int array_size = sizeof(double) * number_of_points;
    int vec_size = sizeof(G2G::vec_type<double,width>) * number_of_points;

    double* rho = NULL;
    double* sigma = NULL;

    double* ex_double = NULL;
    double* ec_double = NULL;
    double* y2a_double = NULL;
    G2G::vec_type<double, width>* grad_double = NULL;
    G2G::vec_type<double, width>* hess1_double = NULL;
    G2G::vec_type<double, width>* hess2_double = NULL;

    err = cudaMalloc((void**)&rho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device rho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void**)&sigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device sigma!\n");
	exit(EXIT_FAILURE);
    }

    // Si el tipo de datos es float, creamos los arrays para copiar
    // los inputs y convertirlos a floats.
    if (!full_double) {
	err = cudaMalloc((void**)&ex_double, array_size);
        if (err != cudaSuccess)
        {
	    fprintf(stderr, "Failed to allocate device ex_double!\n");
    	    exit(EXIT_FAILURE);
        }
	cudaMemset(ex_double,0,array_size);

        err = cudaMalloc((void**)&ec_double, array_size);
	if (err != cudaSuccess)
        {
	    fprintf(stderr, "Failed to allocate device ec_double!\n");
	    exit(EXIT_FAILURE);
        }
	cudaMemset(ec_double,0,array_size);

        err = cudaMalloc((void**)&y2a_double, array_size);
	if (err != cudaSuccess)
        {
    	    fprintf(stderr, "Failed to allocate device y2a_double!\n");
	    exit(EXIT_FAILURE);
        }
	cudaMemset(y2a_double,0,array_size);

        err = cudaMalloc((void**)&grad_double, vec_size);
	if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to allocate device grad_double!\n");
    	    exit(EXIT_FAILURE);
        }

	err = cudaMalloc((void**)&hess1_double, vec_size);
        if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to allocate device hess1_double!\n");
	    exit(EXIT_FAILURE);
        }

	err = cudaMalloc((void**)&hess2_double, vec_size);
        if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to allocate device hess2_double!\n");
	    exit(EXIT_FAILURE);
        }
    }

    // Preparamos los datos.
    if (full_double) {
	err = cudaMemcpy(rho, dens, array_size, cudaMemcpyDeviceToDevice);
        if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to copy data from dens->rho\n");
        }

	err = cudaMemcpy(sigma, contracted_grad, array_size, cudaMemcpyDeviceToDevice);
        if (err != cudaSuccess)
	{
	    fprintf(stderr, "Failed to copy data from contracted_grad->sigma\n");
	}

        // Usamos los datos como vienen ya que son todos doubles.
	ex_double = (double*)ex;
	ec_double = (double*)ec;
        y2a_double = (double*)y2a;
	grad_double = (G2G::vec_type<double,4>*)grad;
        hess1_double = (G2G::vec_type<double,4>*)hess1;
        hess2_double = (G2G::vec_type<double,4>*)hess2;

    } else {
	// Como los inputs son float, los convertimos para libxc
	convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(dens, rho, number_of_points);
	convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(contracted_grad, sigma, number_of_points);
        convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(grad, grad_double, number_of_points);
	convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(hess1, hess1_double, number_of_points);
        convertFloatToDouble<<<blocksPerGrid, threadsPerBlock>>>(hess2, hess2_double, number_of_points);
    }

    // Preparamos los arrays de salida.
    double* exchange = NULL;
    err = cudaMalloc((void **)&exchange, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device exchange!\n");
        exit(EXIT_FAILURE);
    }

    double* correlation = NULL;
    err = cudaMalloc((void **)&correlation, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device correlation!\n");
        exit(EXIT_FAILURE);
    }

    cudaMemset(exchange,0,array_size);
    cudaMemset(correlation,0,array_size);

    // The outputs for exchange
    double* vrho = NULL;
    double* vsigma = NULL;
    double* v2rho = NULL;
    double* v2rhosigma = NULL;
    double* v2sigma = NULL;

    err = cudaMalloc((void **)&vrho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vrho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&vsigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vsigma!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rho, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rho!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhosigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhosigma!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2sigma, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2sigma!\n");
        exit(EXIT_FAILURE);
    }

    // Clear arrays
    cudaMemset(vrho, 0, array_size);
    cudaMemset(vsigma, 0, array_size);
    cudaMemset(v2rho, 0, array_size);
    cudaMemset(v2rhosigma, 0, array_size);
    cudaMemset(v2sigma, 0, array_size);

    // The outputs for correlation
    double* vrhoC = NULL;
    double* vsigmaC = NULL;
    double* v2rhoC = NULL;
    double* v2rhosigmaC = NULL;
    double* v2sigmaC = NULL;

    err = cudaMalloc((void **)&vrhoC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vrhoC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&vsigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vsigmaC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhoC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhoC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2rhosigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2rhosigmaC!\n");
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **)&v2sigmaC, array_size);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device v2sigmaC!\n");
        exit(EXIT_FAILURE);
    }

    ///////////////////////////////////
    // Clear arrays for correlation
    cudaMemset(vrhoC, 0, array_size);
    cudaMemset(vsigmaC, 0, array_size);
    cudaMemset(v2rhoC, 0, array_size);
    cudaMemset(v2rhosigmaC, 0, array_size);
    cudaMemset(v2sigmaC, 0, array_size);

    /////////////////////////////
    // Call LIBXC for exchange
    try {
        xc_gga (&funcForExchange, number_of_points,
                rho,
                sigma,
                exchange,
                vrho,
                vsigma,
                v2rho,
                v2rhosigma,
                v2sigma,
                NULL, NULL, NULL, NULL);
    } catch (int exception) {
        fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
        return;
    }

    ////////////////////////////////
    // Call LIBXC for correlation
    try {
        // Now the correlation value.
        xc_gga (&funcForCorrelation, number_of_points,
                rho,
                sigma,
                correlation,
                vrhoC,
                vsigmaC,
                v2rhoC,
                v2rhosigmaC,
                v2sigmaC,
                NULL, NULL, NULL, NULL);
    } catch (int exception) {
        fprintf (stderr, "Exception ocurred calling xc_gga for Correlation '%d' \n", exception);
        return;
    }

    ////////////////////////
    // Gather the results
    joinResults<T, width><<<blocksPerGrid, threadsPerBlock>>>(
	ex_double, exchange,
	ec_double, correlation,
	vrho, vrhoC,
	vsigma, vsigmaC,
	v2rho, v2rhoC,
	v2rhosigma, v2rhosigmaC,
	v2sigma, v2sigmaC,
	y2a_double,
	sigma,
	grad_double,
	hess1_double,
	hess2_double,
	number_of_points,fact_exchange);

    //////////////////////////
    // Convert if necessary
    if (!full_double) {
	convertDoubleToFloat<<<blocksPerGrid, threadsPerBlock>>> (ex_double, ex, number_of_points);
	convertDoubleToFloat<<<blocksPerGrid, threadsPerBlock>>> (ec_double, ec, number_of_points);
        convertDoubleToFloat<<<blocksPerGrid, threadsPerBlock>>> (y2a_double, y2a, number_of_points);
    }

    /////////////////////////
    // Free device memory
    if (rho != NULL) {
	cudaFree(rho);
    }
    if (sigma != NULL) {
	cudaFree(sigma);
    }
    if (exchange != NULL) {
	cudaFree(exchange);
    }
    if (correlation != NULL) {
	cudaFree(correlation);
    }
    if (vrho != NULL) {
        cudaFree(vrho);
    }
    if (vsigma != NULL) {
	cudaFree(vsigma);
    }
    if (v2rho != NULL) {
	cudaFree(v2rho);
    }
    if (v2rhosigma != NULL) {
	cudaFree(v2rhosigma);
    }
    if (v2sigma != NULL) {
	cudaFree(v2sigma);
    }
    if (vrhoC != NULL) {
        cudaFree(vrhoC);
    }
    if (vsigmaC != NULL) {
	cudaFree(vsigmaC);
    }
    if (v2rhoC != NULL) {
	cudaFree(v2rhoC);
    }
    if (v2rhosigmaC != NULL) {
	cudaFree(v2rhosigmaC);
    }
    if (v2sigmaC != NULL) {
	cudaFree(v2sigmaC);
    }

    if (!full_double) {
        if (ex_double != NULL) {
            cudaFree(ex_double);
        }
        if (ec_double != NULL) {
            cudaFree(ec_double);
        }
        if (y2a_double != NULL) {
        cudaFree(y2a_double);
        }
        if (grad_double != NULL) {
            cudaFree((void*)grad_double);
    	}
    	if (hess1_double != NULL) {
        cudaFree((void*)hess1_double);
    	}
        if (hess2_double != NULL) {
    	    cudaFree((void*)hess2_double);
        }
    }
#endif
    return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy::doLDA - GPU version
// Calls the xc_lda function from libxc
// dens: pointer for the density array
// number_of_points: the size of all the input arrays
// contracted_grad: the contracted grad for libxc
// grad:.
// hess1:
// hess2:
// ex: here goes the results after calling xc_gga from libxc for the exchange functional
// ec: here goes the results after calling xc_gga from libxc for the correlation functional
//
// Note: all the pointer data are pointers in CUDA memory.
//

template <class T, int width>
void LibxcProxy <T, width>::doLDA(T dens, const G2G::vec_type<T, width> &grad, const G2G::vec_type<T, width> &hess1, const G2G::vec_type<T, width> &hess2, T &ex, T &ec, T &y2a)
{
    //TODO: not implemented yet!
    return;
}

#endif // LIBXCPROXY_H
