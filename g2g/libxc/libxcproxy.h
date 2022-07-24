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
    // Inner Variables
    xc_func_type* funcsId;
    int ntotal_funcs, nxcoef, nccoef, nspin, nsrid;
    double omega;
    double* funcsCoef;

    // Is inited
    bool inited;

//  TODO: borrar esto, es solo pa que compile por partes
    xc_func_type funcForExchange, funcForCorrelation;
    double fact_exchange;
///////////////////

public:

    // Constructors
    LibxcProxy ();
    LibxcProxy (int* func_id,double* func_coef, int nx_coef, int nc_coef,
                                int nsr_id, double screen, int nSpin);
    void init (int* func_id,double* func_coef, int nx_coef, int nc_coef,
                                int nsr_id, double screen, int nSpin);
    // Destructors
    ~LibxcProxy ();
    void closeProxy ();

    // Libxc Calculations
    void doSCF(T& dens, const G2G::vec_type<T, width>& grad,
               const G2G::vec_type<T, width>& hess1,
               const G2G::vec_type<T, width>& hess2,
               T& ex, T& ec, T& y2a);

    // Open doSCF
    void doSCF(T& dens_a, T& dens_b, 
               const G2G::vec_type<T, width>& grad_a, const G2G::vec_type<T, width>& grad_b,
               T& ex, T& ec, T* coef_a, T* coef_b);

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

    void terms_derivs(double* dens, double* sigma,
                double &vrho, double &vsigma, double &v2rho2,
                double &v2rhosigma, double &v2sigma2, double &v3rho3,
                double &v3rho2sigma, double &v3rhosigma2, double &v3sigma3);

    void doLDA (T dens,
                const G2G::vec_type<T,width>& grad,
                const G2G::vec_type<T,width>& hess1,
                const G2G::vec_type<T,width>& hess2,
                T& ex,
                T& ec,
                T& y2a);
};

template <class T, int width>
LibxcProxy <T, width>::LibxcProxy()
{
   ntotal_funcs = -1;
   nxcoef = -1; nccoef = -1;
   nspin = -1;
   nsrid = -1;
   omega = -1.0f;
   inited = false;
}

template <class T, int width>
LibxcProxy <T, width>::LibxcProxy (int* func_id,double* func_coef, int nx_coef, int nc_coef,
                                   int nsr_id, double screen, int nSpin)
{
    init(func_id,func_coef,nx_coef,nc_coef,nsr_id,screen,nSpin);
}

template <class T, int width>
void LibxcProxy<T, width>::init(int* func_id,double* func_coef, int nx_coef, int nc_coef,
                                int nsr_id, double screen, int nSpin)
{
    // Set inner libxc variables
    ntotal_funcs = nx_coef + nc_coef;
    nxcoef = nx_coef; nccoef = nc_coef;
    nspin  = nSpin; nsrid = nsr_id; omega = screen;
    funcsId   = (xc_func_type*) malloc(ntotal_funcs*sizeof(xc_func_type));
    funcsCoef = (double*) malloc(ntotal_funcs*sizeof(double));
    double threshold = 1e-20;

    // Init of differents functionals needed in the simulation
    for (int ii=0; ii<ntotal_funcs; ii++) {
       if (xc_func_init (&funcsId[ii], func_id[ii], nspin) != 0) {
           fprintf (stderr, "Functional '%d' not found\n", func_id[ii]);
           exit(-1);
       }
       if ( ii == nsrid ) {
          xc_func_set_ext_params(&funcsId[ii],&omega);
       } else {
          xc_func_set_dens_threshold(&funcsId[ii],threshold);
       }
       funcsCoef[ii] = func_coef[ii];
    }

    inited = true;
}

template <class T, int width>
LibxcProxy <T, width>::~LibxcProxy ()
{
    closeProxy ();
}

template <class T, int width>
void LibxcProxy <T, width>::closeProxy ()
{
    if (inited) {
        for (int ii=0; ii<ntotal_funcs; ii++)
           xc_func_end(&funcsId[ii]);

        free(funcsId);   funcsId   = NULL;
        free(funcsCoef); funcsCoef = NULL;
	inited = false;
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
void LibxcProxy <T, width>::doSCF(T& dens,
    const G2G::vec_type<T, width>& grad,
    const G2G::vec_type<T, width>& hess1,
    const G2G::vec_type<T, width>& hess2,
    T& ex, T& ec, T& y2a)
{
    const double rho[1] = {dens};
    double sigma[1] = {(grad.x * grad.x) + (grad.y * grad.y) + (grad.z * grad.z)};

    //All outputs
    ex = ec = 0.0f;
    double vrho=0.0f;
    double vsigma=0.0f;
    double v2rho2=0.0f;
    double v2rhosigma=0.0f;
    double v2sigma2=0.0f;

    // Functional outputs
    double lenergy[1];
    double lvrho[1];
    double lvsigma[1];
    double lv2rho2[1];
    double lv2rhosigma[1];
    double lv2sigma2[1];

    double cc = 0.0f;

    // Exchange Calculation
    for (int ii=0; ii<nxcoef; ii++) {
        // Set zero values
        lenergy[0]  = 0.0f;
        lvrho[0]    = 0.0f;
        lv2rho2[0]  = 0.0f;
        lvsigma[0]  = 0.0f;
        lv2rhosigma[0] = 0.0f;
        lv2sigma2[0]   = 0.0f;

        switch( (&funcsId[ii])->info->family ) {
           case (XC_FAMILY_LDA):
              xc_lda(&funcsId[ii],1,rho,lenergy,lvrho,lv2rho2,
                     NULL,
                     NULL); break;
           case (XC_FAMILY_GGA):
              xc_gga(&funcsId[ii],1,rho,sigma,lenergy,lvrho,lvsigma,
                     lv2rho2,lv2rhosigma,lv2sigma2,
                     NULL,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL); break;
           default:
             printf("Unidentified Family Functional\n");
             exit(-1); break;

        } // end switch

        cc = funcsCoef[ii];
        ex         += cc*lenergy[0];
        vrho       += cc*lvrho[0];
        v2rho2     += cc*lv2rho2[0];
        vsigma     += cc*lvsigma[0];
        v2rhosigma += cc*lv2rhosigma[0];
        v2sigma2   += cc*lv2sigma2[0];
    } // end exchange

    // Correlation Calculation
    for (int ii=nxcoef; ii<ntotal_funcs; ii++) {
        // Set zero values
        lenergy[0]  = 0.0f;
        lvrho[0]    = 0.0f;
        lv2rho2[0]  = 0.0f;
        lvsigma[0]  = 0.0f;
        lv2rhosigma[0] = 0.0f;
        lv2sigma2[0]   = 0.0f;

        switch( (&funcsId[ii])->info->family ) {
           case (XC_FAMILY_LDA):
              xc_lda(&funcsId[ii],1,rho,lenergy,lvrho,lv2rho2,
                     NULL,
                     NULL); break;
           case (XC_FAMILY_GGA):
              xc_gga(&funcsId[ii],1,rho,sigma,lenergy,lvrho,lvsigma,
                     lv2rho2,lv2rhosigma,lv2sigma2,
                     NULL,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL); break;
           default:
             printf("Unidentified Family Functional\n");
             exit(-1); break;

        } // end switch

        cc = funcsCoef[ii];
        ec         += cc*lenergy[0];
        vrho       += cc*lvrho[0];
        v2rho2     += cc*lv2rho2[0];
        vsigma     += cc*lvsigma[0];
        v2rhosigma += cc*lv2rhosigma[0];
        v2sigma2   += cc*lv2sigma2[0];
    } // end correlation

    // Now, compute y2a value.
    y2a = vrho - (2.0f * sigma[0] * v2rhosigma + 2.0f * (hess1.x + hess1.y + hess1.z) * vsigma
               + 4.0f * v2sigma2 * (grad.x * grad.x * hess1.x + 
                                   grad.y * grad.y * hess1.y +
                                   grad.z * grad.z * hess1.z + 
                            2.0f * grad.x * grad.y * hess2.x + 
                            2.0f * grad.x * grad.z * hess2.y +
                            2.0f * grad.y * grad.z * hess2.z));
    return;
}

// OPEN version of doSCF
template <class T, int width>
void LibxcProxy <T, width>::doSCF(T& dens_a, T& dens_b,
    const G2G::vec_type<T, width>& grad_a, const G2G::vec_type<T, width>& grad_b,
    // Outputs
    T& ex, T& ec, T* coef_a, T* coef_b)
{
    int size = sizeof(double);

    double rho[2], sigma[3];
    rho[0] = dens_a; rho[1] = dens_b;
    sigma[0] = grad_a.x*grad_a.x + grad_a.y*grad_a.y + grad_a.z*grad_a.z; // alpha
    sigma[1] = grad_a.x*grad_b.x + grad_a.y*grad_b.y + grad_a.z*grad_b.z; // alpha-beta
    sigma[2] = grad_b.x*grad_b.x + grad_b.y*grad_b.y + grad_b.z*grad_b.z; // beta

    // Final Ouputs Libxc
    double *vrho, *vsigma;
    vrho   = (double*)malloc(2*size); memset(vrho,0.0f,2*size);
    vsigma = (double*)malloc(3*size); memset(vsigma,0.0f,3*size);

    // Local outputs libxc
    double *lenergy, *lvrho, *lvsigma, *lv2rho2, *lv2sigma2, *lv2rhosigma;
    lenergy     = (double*)malloc(1*size);
    lvrho       = (double*)malloc(2*size);
    lvsigma     = (double*)malloc(3*size);
    lv2rho2     = (double*)malloc(3*size);
    lv2sigma2   = (double*)malloc(6*size);
    lv2rhosigma = (double*)malloc(6*size);
    double cc = 0.0f; ex = ec = 0.0f;

    // Exchange Calculation
    for (int ii=0; ii<nxcoef; ii++) {
        // Set zero values
        memset(lenergy,0.0f,1*size);
        memset(lvrho,0.0f,2*size);
        memset(lvsigma,0.0f,3*size);
        memset(lv2rho2,0.0f,3*size);
        memset(lv2sigma2,0.0f,6*size);
        memset(lv2rhosigma,0.0f,6*size);

        switch( (&funcsId[ii])->info->family ) {
           case (XC_FAMILY_LDA):
              xc_lda(&funcsId[ii],1,rho,lenergy,lvrho,NULL,
                     NULL,
                     NULL); break;
           case (XC_FAMILY_GGA):
              xc_gga(&funcsId[ii],1,rho,sigma,lenergy,lvrho,lvsigma,
                     lv2rho2,lv2rhosigma,lv2sigma2,
                     NULL,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL); break;
           default:
             printf("Unidentified Family Functional\n");
             exit(-1); break;
        } // end switch
        cc  = funcsCoef[ii];
        ex += cc*lenergy[0];

        // there is no cross alpha-beta terms in exchange functional
        vrho[0] += cc*lvrho[0]; vrho[1] += cc*lvrho[1];
        vsigma[0] += cc*lvsigma[0]; vsigma[2] += cc*lvsigma[2];
    } // end exchange

    // Correlation Calculation
    for (int ii=nxcoef; ii<ntotal_funcs; ii++) {
        // Set zero values
        memset(lenergy,0.0f,1*size);
        memset(lvrho,0.0f,2*size);
        memset(lvsigma,0.0f,3*size);
        memset(lv2rho2,0.0f,3*size);
        memset(lv2sigma2,0.0f,6*size);
        memset(lv2rhosigma,0.0f,6*size);

        switch( (&funcsId[ii])->info->family ) {
           case (XC_FAMILY_LDA):
              xc_lda(&funcsId[ii],1,rho,lenergy,lvrho,lv2rho2,
                     NULL,
                     NULL); break;
           case (XC_FAMILY_GGA):
              xc_gga(&funcsId[ii],1,rho,sigma,lenergy,lvrho,lvsigma,
                     lv2rho2,lv2rhosigma,lv2sigma2,
                     NULL,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL); break;
           default:
             printf("Unidentified Family Functional\n");
             exit(-1); break;
        } // end switch
        cc  = funcsCoef[ii];
        ec += cc*lenergy[0];

        // alpha and beta terms in correlation
        vrho[0] += cc*lvrho[0]; vrho[1] += cc*lvrho[1];
        vsigma[0] += cc*lvsigma[0]; vsigma[2] += cc*lvsigma[2];

        // cross alpha-beta terms in correlation
        vsigma[1] += cc*lvsigma[1];

    } // end correlation

    // return values
    coef_a[0] = vrho[0]; coef_a[1] = vsigma[0]; coef_a[2] = vsigma[1];
    coef_b[0] = vrho[1]; coef_b[1] = vsigma[2]; coef_b[2] = vsigma[1];

    // free outputs
    free(vrho); vrho = NULL;
    free(vsigma); vsigma = NULL;

    // free local outputs
    free(lenergy); lenergy = NULL;
    free(lvrho); lvrho = NULL;
    free(lvsigma); lvsigma = NULL;
    free(lv2rho2); lv2rho2 = NULL;
    free(lv2sigma2); lv2sigma2 = NULL;
    free(lv2rhosigma); lv2rhosigma = NULL;
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

/* TODO
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
*/

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
                double* sigma,
                double red,
                double cruz,
                double* lrCoef)
{
    //All outputs
    double vrho=0.0f;
    double vsigma=0.0f;
    double v2rho2=0.0f;
    double v2rhosigma=0.0f;
    double v2sigma2=0.0f;

    // Functional outputs
    double lenergy[1];
    double lvrho[1];
    double lvsigma[1];
    double lv2rho2[1];
    double lv2rhosigma[1];
    double lv2sigma2[1];

    double cc = 0.0f;

    // Exchange Calculation
    for (int ii=0; ii<nxcoef; ii++) {
        // Set zero values
        lenergy[0]  = 0.0f;
        lvrho[0]    = 0.0f;
        lv2rho2[0]  = 0.0f;
        lvsigma[0]  = 0.0f;
        lv2rhosigma[0] = 0.0f;
        lv2sigma2[0]   = 0.0f;

        switch( (&funcsId[ii])->info->family ) {
           case (XC_FAMILY_LDA):
              xc_lda(&funcsId[ii],1,rho,lenergy,lvrho,lv2rho2,
                     NULL,
                     NULL); break;
           case (XC_FAMILY_GGA):
              xc_gga(&funcsId[ii],1,rho,sigma,lenergy,lvrho,lvsigma,
                     lv2rho2,lv2rhosigma,lv2sigma2,
                     NULL,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL); break;
           default:
             printf("Unidentified Family Functional\n");
             exit(-1); break;

        } // end switch

        cc = funcsCoef[ii];
        vrho       += cc*lvrho[0];
        v2rho2     += cc*lv2rho2[0];
        vsigma     += cc*lvsigma[0];
        v2rhosigma += cc*lv2rhosigma[0];
        v2sigma2   += cc*lv2sigma2[0];
    } // end exchange

    // Correlation Calculation
    for (int ii=nxcoef; ii<ntotal_funcs; ii++) {
        // Set zero values
        lenergy[0]  = 0.0f;
        lvrho[0]    = 0.0f;
        lv2rho2[0]  = 0.0f;
        lvsigma[0]  = 0.0f;
        lv2rhosigma[0] = 0.0f;
        lv2sigma2[0]   = 0.0f;

        switch( (&funcsId[ii])->info->family ) {
           case (XC_FAMILY_LDA):
              xc_lda(&funcsId[ii],1,rho,lenergy,lvrho,lv2rho2,
                     NULL,
                     NULL); break;
           case (XC_FAMILY_GGA):
              xc_gga(&funcsId[ii],1,rho,sigma,lenergy,lvrho,lvsigma,
                     lv2rho2,lv2rhosigma,lv2sigma2,
                     NULL,NULL,NULL,NULL,
                     NULL,NULL,NULL,NULL,NULL); break;
           default:
             printf("Unidentified Family Functional\n");
             exit(-1); break;

        } // end switch

        cc = funcsCoef[ii];
        vrho       += cc*lvrho[0];
        v2rho2     += cc*lv2rho2[0];
        vsigma     += cc*lvsigma[0];
        v2rhosigma += cc*lv2rhosigma[0];
        v2sigma2   += cc*lv2sigma2[0];
    } // end correlation

   // Results
   lrCoef[0] = red * v2rho2  + 2.0f * cruz * v2rhosigma ;
   lrCoef[1] = 2.0f * red * v2rhosigma  + 4.0f * cruz * v2sigma2 ;
   lrCoef[2] = 2.0f * vsigma ;

   return;
}

template <class T, int width>
void LibxcProxy<T,width>::terms_derivs(double* rho, double* sigma,
               double &vrho, double &vsigma, double &v2rho2,
               double &v2rhosigma, double &v2sigma2, double &v3rho3,
               double &v3rho2sigma, double &v3rhosigma2, double &v3sigma3)
{
   // All outputs
   vrho = vsigma = v2rho2 = v2rhosigma = 0.0f;
   v2sigma2 = v3rho3 = v3rho2sigma = 0.0f;
   v3rhosigma2 = v3sigma3 =0.0f;

   // Local otputs libxc
   double lvrho[1], lvsigma[1], lv2rho2[1], lv2rhosigma[1], lv2sigma2[1];
   double lv3rho3[1], lv3rho2sigma[1], lv3rhosigma2[1], lv3sigma3[1], lenergy[1];
   double cc = 0.0f;

   // Exchange Calculation
   for (int ii=0; ii<nxcoef; ii++) {
       // Set zero values
       lenergy[0]  = 0.0f;
       lvrho[0]    = 0.0f;
       lv2rho2[0]  = 0.0f;
       lvsigma[0]  = 0.0f;
       lv2rhosigma[0] = 0.0f;
       lv2sigma2[0]   = 0.0f;
       lv3rho3[0]     = 0.0f;
       lv3rho2sigma[0]= 0.0f;
       lv3rhosigma2[0]= 0.0f;
       lv3sigma3[0]   = 0.0f;

       switch( (&funcsId[ii])->info->family ) {
          case (XC_FAMILY_LDA):
             xc_lda(&funcsId[ii],1,rho,lenergy,lvrho,lv2rho2,
                    lv3rho3,
                    NULL); break;
          case (XC_FAMILY_GGA):
             xc_gga(&funcsId[ii],1,rho,sigma,lenergy,lvrho,lvsigma,
                    lv2rho2,lv2rhosigma,lv2sigma2,
                    lv3rho3,lv3rho2sigma,lv3rhosigma2,lv3sigma3,
                    NULL,NULL,NULL,NULL,NULL); break;
          default:
            printf("Unidentified Family Functional\n");
            exit(-1); break;

       } // end switch

       cc = funcsCoef[ii];
       vrho       += cc*lvrho[0];
       vsigma     += cc*lvsigma[0];
       v2rho2     += cc*lv2rho2[0];
       v2rhosigma += cc*lv2rhosigma[0];
       v2sigma2   += cc*lv2sigma2[0];
       v3rho3     += cc*lv3rho3[0];
       v3rho2sigma+= cc*lv3rho2sigma[0];
       v3sigma3   += cc*lv3sigma3[0];
       v3rhosigma2+= cc*lv3rhosigma2[0];
   } // end exchange

   // Correlation Calculation
   for (int ii=nxcoef; ii<ntotal_funcs; ii++) {
       // Set zero values
       lenergy[0]  = 0.0f;
       lvrho[0]    = 0.0f;
       lv2rho2[0]  = 0.0f;
       lvsigma[0]  = 0.0f;
       lv2rhosigma[0] = 0.0f;
       lv2sigma2[0]   = 0.0f;
       lv3rho3[0]     = 0.0f;
       lv3rho2sigma[0]= 0.0f;
       lv3rhosigma2[0]= 0.0f;
       lv3sigma3[0]   = 0.0f;

       switch( (&funcsId[ii])->info->family ) {
          case (XC_FAMILY_LDA):
             xc_lda(&funcsId[ii],1,rho,lenergy,lvrho,lv2rho2,
                    lv3rho3,
                    NULL); break;
          case (XC_FAMILY_GGA):
             xc_gga(&funcsId[ii],1,rho,sigma,lenergy,lvrho,lvsigma,
                    lv2rho2,lv2rhosigma,lv2sigma2,
                    lv3rho3,lv3rho2sigma,lv3rhosigma2,lv3sigma3,
                    NULL,NULL,NULL,NULL,NULL); break;
          default:
            printf("Unidentified Family Functional\n");
            exit(-1); break;

       } // end switch

       cc = funcsCoef[ii];
       vrho       += cc*lvrho[0];
       vsigma     += cc*lvsigma[0];
       v2rho2     += cc*lv2rho2[0];
       v2rhosigma += cc*lv2rhosigma[0];
       v2sigma2   += cc*lv2sigma2[0];
       v3rho3     += cc*lv3rho3[0];
       v3rho2sigma+= cc*lv3rho2sigma[0];
       v3sigma3   += cc*lv3sigma3[0];
       v3rhosigma2+= cc*lv3rhosigma2[0];
   } // end correlation
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

   // All outputs
   double vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, v3rho3, v3rho2sigma;
   double v3rhosigma2, v3sigma3, exc;
   vrho = vsigma = 0.0f; v2rho2 = v2rhosigma = v2sigma2 = 0.0f;
   v3rho3 = v3rho2sigma = v3rhosigma2 = v3sigma3 = exc = 0.0f;

   // Functional outputs
   double lenergy[1];
   double lvrho[1];
   double lvsigma[1];
   double lv2rho2[1];
   double lv2rhosigma[1];
   double lv2sigma2[1];
   double lv3rho3[1];
   double lv3rho2sigma[1];
   double lv3rhosigma2[1];
   double lv3sigma3[1];

   double cc = 0.0f;

   // Exchange Calculation
   for (int ii=0; ii<nxcoef; ii++) {
       // Set zero values
       lenergy[0]  = 0.0f;
       lvrho[0]    = 0.0f;
       lv2rho2[0]  = 0.0f;
       lvsigma[0]  = 0.0f;
       lv2rhosigma[0] = 0.0f;
       lv2sigma2[0]   = 0.0f;
       lv3rho3[0]     = 0.0f;
       lv3rho2sigma[0]= 0.0f;
       lv3rhosigma2[0]= 0.0f;
       lv3sigma3[0]   = 0.0f;

       switch( (&funcsId[ii])->info->family ) {
          case (XC_FAMILY_LDA):
             xc_lda(&funcsId[ii],1,&dens,lenergy,lvrho,lv2rho2,
                    lv3rho3,
                    NULL); break;
          case (XC_FAMILY_GGA):
             xc_gga(&funcsId[ii],1,&dens,&sigma,lenergy,lvrho,lvsigma,
                    lv2rho2,lv2rhosigma,lv2sigma2,
                    lv3rho3,lv3rho2sigma,lv3rhosigma2,lv3sigma3,
                    NULL,NULL,NULL,NULL,NULL); break;
          default:
            printf("Unidentified Family Functional\n");
            exit(-1); break;

       } // end switch

       cc = funcsCoef[ii];
       v2rho2     += cc*lv2rho2[0];
       v2rhosigma += cc*lv2rhosigma[0];
       v2sigma2   += cc*lv2sigma2[0];
       v3rho3     += cc*lv3rho3[0];
       v3rho2sigma+= cc*lv3rho2sigma[0];
       v3rhosigma2+= cc*lv3rhosigma2[0];
       v3sigma3   += cc*lv3sigma3[0];
       
   } // end exchange

   // Correlation Calculation
   for (int ii=nxcoef; ii<ntotal_funcs; ii++) {
       // Set zero values
       lenergy[0]  = 0.0f;
       lvrho[0]    = 0.0f;
       lv2rho2[0]  = 0.0f;
       lvsigma[0]  = 0.0f;
       lv2rhosigma[0] = 0.0f;
       lv2sigma2[0]   = 0.0f;
       lv3rho3[0]     = 0.0f;
       lv3rho2sigma[0]= 0.0f;
       lv3rhosigma2[0]= 0.0f;
       lv3sigma3[0]   = 0.0f;

       switch( (&funcsId[ii])->info->family ) {
          case (XC_FAMILY_LDA):
             xc_lda(&funcsId[ii],1,&dens,lenergy,lvrho,lv2rho2,
                    lv3rho3,
                    NULL); break;
          case (XC_FAMILY_GGA):
             xc_gga(&funcsId[ii],1,&dens,&sigma,lenergy,lvrho,lvsigma,
                    lv2rho2,lv2rhosigma,lv2sigma2,
                    lv3rho3,lv3rho2sigma,lv3rhosigma2,lv3sigma3,
                    NULL,NULL,NULL,NULL,NULL); break;
          default:
            printf("Unidentified Family Functional\n");
            exit(-1); break;

       } // end switch
       cc = funcsCoef[ii];
       v2rhosigma += cc*lv2rhosigma[0];
       v2sigma2   += cc*lv2sigma2[0];
       v3rho3     += cc*lv3rho3[0];
       v3rho2sigma+= cc*lv3rho2sigma[0];
       v3rhosigma2+= cc*lv3rhosigma2[0];
       v3sigma3   += cc*lv3sigma3[0];
   } // end correlation

   double prod_rr, prod_pr;
   double term1, term2, term3, term4;

   prod_rr = txyz[0]*txyz[0] + txyz[1]*txyz[1] + txyz[2]*txyz[2];
   prod_pr = dxyz[0]*txyz[0] + dxyz[1]*txyz[1] + dxyz[2]*txyz[2];

   term1 = red * red * v3rho3;
   term2 = 4.0f * red * prod_pr * v3rho2sigma;
   term3 = 4.0f * prod_pr * prod_pr * v3rhosigma2;
   term4 = 2.0f * prod_rr * v2rhosigma;
   zcoef[0] = term1 + term2 + term3 + term4;
   zcoef[0] *= 2.0f;

   term1 = 2.0f * red * red * v3rho2sigma;
   term2 = 8.0f * red * prod_pr * v3rhosigma2;
   term3 = 8.0f * prod_pr * prod_pr * v3sigma3;
   term4 = 4.0f * prod_rr * v2sigma2;
   zcoef[1] = term1 + term2 + term3 + term4;
   zcoef[1] *= 2.0f;

   zcoef[2] = 4.0f * red * v2rhosigma + 8.0f * prod_pr * v2sigma2;
   zcoef[2] *= 2.0f;

   free(dxyz); dxyz = NULL;
   free(txyz); txyz = NULL;

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
