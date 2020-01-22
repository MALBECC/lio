#ifndef LIBXCPROXY_CUDA_H
#define LIBXCPROXY_CUDA_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include <xc.h>
#include <xc_cuda.h>
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
class LibxcProxy_cuda
{
private:

    // The libxc components
    xc_func_type_cuda funcForExchange;
    xc_func_type_cuda funcForCorrelation;

    // Functional ids
    int funcIdForExchange;
    int funcIdForCorrelation;
    int nspin;

    // Is inited
    bool inited_cuda;

    // Exchange Exact
    double fact_exchange;

    void printFunctionalInformation (xc_func_type_cuda* func);

public:
    LibxcProxy_cuda ();
    LibxcProxy_cuda (int exchangeFunctionId, int correlationFuncionalId, int nspin, double fexc);
    ~LibxcProxy_cuda ();

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

    void coefLR(double* rho,double* sigma,
                double red,double cruz,double* lrCoef); // CPU

    void coefLR(const int npoints, double* rho,double* sigma,
                double* red,double* cruz,double* lrCoef); // GPU

    void coefZv(const int npoints, double* rho, double* sigma,
                double* red, G2G::vec_type<T,WIDTH>* dxyz, G2G::vec_type<T,WIDTH>* txyz,
                double* lrCoef);

    void obtain_gpu_der(const int npoints, T* rho, T* sigma, 
                        G2G::vec_type<T,2>* vrho,
                        G2G::vec_type<T,2>* vsigma,
                        G2G::vec_type<T,2>* v2rho2,
                        G2G::vec_type<T,2>* v2rhosigma,
                        G2G::vec_type<T,2>* v2sigma2,
                        G2G::vec_type<T,2>* v3rho3,
                        G2G::vec_type<T,2>* v3rho2sigma,
                        G2G::vec_type<T,2>* v3rhosigma2,
                        G2G::vec_type<T,2>* v3sigma3);

    void doLDA (T dens,
                const G2G::vec_type<T,width>& grad,
                const G2G::vec_type<T,width>& hess1,
                const G2G::vec_type<T,width>& hess2,
                T& ex,
                T& ec,
                T& y2a);

    void init_cuda (int exId, int xcId, int nspin);
    void closeProxy_cuda ();
    void printFunctionalsInformation (int exchangeFunctionalId, int correlationFunctionalId);
};

template <class T, int width>
LibxcProxy_cuda <T, width>::LibxcProxy_cuda()
{
    funcIdForExchange = 0;
    funcIdForCorrelation = 0;
    nspin = 0;
    inited_cuda = false;
}

template <class T, int width>
LibxcProxy_cuda <T, width>::LibxcProxy_cuda (int exchangeFunctionalId, int correlationFuncionalId, int nSpin, double fexc)
{
//    printf("LibxcProxy_cuda::LibxcProxy_cuda (%u, %u, %u) \n", exchangeFunctionalId, correlationFuncionalId, nSpin);
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
    init_cuda (exchangeFunctionalId, correlationFuncionalId, nSpin);
}

template <class T, int width>
void LibxcProxy_cuda<T, width>::init_cuda (int exchangeFunctionalId, int correlationFunctionalId, int nSpin)
{
//    printf("LibxcProxy_cuda::init(%u, %u, %u)\n", exchangeFunctionalId, correlationFunctionalId, nSpin);
    funcIdForExchange = exchangeFunctionalId;
    funcIdForCorrelation = correlationFunctionalId;
    nspin = nSpin;

    //printf("libxcproxy_cuda.h XC_GGA_X_PBE_cuda XC_GGA_C_PBE_cuda %d %d\n",XC_GGA_X_PBE_cuda,XC_GGA_C_PBE_cuda);
    //printf("libxcproxy_cuda.h exID coID %d %d\n",exchangeFunctionalId,correlationFunctionalId);


    if (xc_func_init_cuda (&funcForExchange, funcIdForExchange, nspin) != 0) {
        fprintf (stderr, "cuda X Functional '%d' not found\n", funcIdForExchange);
	exit(-1);
    }

    if (xc_func_init_cuda (&funcForCorrelation, funcIdForCorrelation, nspin) != 0){
	fprintf (stderr, "cuda E Functional '%d' not found\n", funcIdForCorrelation);
	exit(-1);
    }

    inited_cuda = true;
}

template <class T, int width>
LibxcProxy_cuda <T, width>::~LibxcProxy_cuda ()
{
    //xc_func_end (&funcForExchange);
    //xc_func_end (&funcForCorrelation);
//    printf("LibxcProxy_cuda::~LibxcProxy_cuda()\n");
    closeProxy_cuda ();
}

template <class T, int width>
void LibxcProxy_cuda <T, width>::closeProxy_cuda ()
{
//    printf("LibxcProxy_cuda::closeProxy()\n");
    if (inited_cuda) {
	xc_func_end_cuda (&funcForExchange);
        xc_func_end_cuda (&funcForCorrelation);
	inited_cuda = false;
    }
}

template <class T, int width>
void LibxcProxy_cuda<T, width>::printFunctionalsInformation (int exchangeFunctionalId, int correlationFunctionalId) 
{
    if (!inited_cuda) 
    {
	init_cuda (exchangeFunctionalId, correlationFunctionalId, 1);
    }

    printFunctionalInformation (&funcForExchange);
    printFunctionalInformation (&funcForCorrelation);

    if (inited_cuda) 
    {
	closeProxy_cuda ();
    }
}

template <class T, int width>
void LibxcProxy_cuda<T, width>::printFunctionalInformation (xc_func_type_cuda* func)
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
// LibxcProxy_cuda::doGGA - CPU Version 1
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
void LibxcProxy_cuda <T, width>::doGGA(T dens,
    const G2G::vec_type<T, width> &grad,
    const G2G::vec_type<T, width> &hess1,
    const G2G::vec_type<T, width> &hess2,
    T &ex, T &ec, T &y2a)
{
    //printf("LibxcProxy_cuda::doGGA cpu simple(...) \n");

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
    double y2a_x;

    // The outputs for correlation
    double vrhoC [1];
    double vsigmaC [1];
    double v2rhoC [1];
    double v2rhosigmaC [1];
    double v2sigmaC [1];
    double y2a_c;

    // The exchange values
    xc_gga_cuda (&funcForExchange, 1,
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
    xc_gga_cuda (&funcForCorrelation, 1,
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
/*
    vrho[0] += vrhoC[0];
    vsigma[0] += vsigmaC[0];
    v2rho[0] += v2rhoC[0];
    v2rhosigma[0] += v2rhosigmaC[0];
    v2sigma[0] += v2sigmaC[0];
*/

    // Now, compute y2a value.
    y2a_x = vrho[0] - (2 * sigma[0] * v2rhosigma[0]
            + 2 * (hess1.x + hess1.y + hess1.z) * vsigma[0]
            + 4 * v2sigma[0] * (grad.x * grad.x * hess1.x + grad.y * grad.y * hess1.y + grad.z * grad.z * hess1.z + 2 * grad.x * grad.y * hess2.x + 2 * grad.x * grad.z * hess2.y + 2 * grad.y * grad.z * hess2.z));


    y2a_c = vrhoC[0] - (2 * sigma[0] * v2rhosigmaC[0]
            + 2 * (hess1.x + hess1.y + hess1.z) * vsigmaC[0]
            + 4 * v2sigmaC[0] * (grad.x * grad.x * hess1.x + grad.y * grad.y * hess1.y + grad.z * grad.z * hess1.z + 2 * grad.x * grad.y * hess2.x + 2 * grad.x * grad.z * hess2.y + 2 * grad.y * grad.z * hess2.z));

    y2a = fact_exchange * y2a_x + y2a_c;

    return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy_cuda::doGGA - CPU Version 2
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
void LibxcProxy_cuda <T, width>::doGGA(T* dens,
    const int number_of_points,
    const G2G::vec_type<T, width>* grad,
    const G2G::vec_type<T, width>* hess1,
    const G2G::vec_type<T, width>* hess2,
    T* ex,
    T* ec,
    T* y2a)
{
    //printf("LibxcProxy_cuda::doGGA cpu multiple (...) \n");

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
    xc_gga_cuda (&funcForExchange, number_of_points,
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
    xc_gga_cuda (&funcForCorrelation, number_of_points,
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
// LibxcProxy_cuda::doGGA - CPU Version 3
// Calls the XC_GGA_FXC function from libxc for multiple points.
// rho: pointer for the density array.
// sigma: contracted gradient array.
// number_of_points: the size of all the input arrays.
// v2rho2: second partial derivative of the energy per unit volume in terms of the density.
// v2rhosigma: second partial derivative of the energy per unit volume in terms of the density and sigma.
// v2sigma2: second partial derivative of the energy per unit volume in terms of sigma.
//
template <class T, int width>
void LibxcProxy_cuda <T, width>::doGGA (T* rho,
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
    xc_gga_fxc_cuda (&funcForExchange, number_of_points,
                rho,
                sigma,
                v2rho2X,
                v2rhosigmaX,
                v2sigma2X);

    // Now the correlation value.
    xc_gga_fxc_cuda (&funcForCorrelation, number_of_points,
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

#ifdef __CUDACC__

template <class T, int width>
__global__ void local_coef_gpu(double* tred, double* cruz, double* Coef,double* v2rho2,
                               double* v2rhosigma,double* v2sigma2,double* vsigma, 
                               const int npoints,double fact_ex,const bool exchange) 
{
   double fex;
   int i = blockDim.x * blockIdx.x + threadIdx.x;

   fex = fact_ex;

   if ( i < npoints ) {

     if ( exchange ) { // EXCHANGE

       Coef[i] = (2.0f * tred[i] * v2rho2[i] + 8.0f * cruz[i] * v2rhosigma[i]) * fex;
       Coef[npoints+i] = (8.0f * tred[i] * v2rhosigma[i] + 32.0f * cruz[i] * v2sigma2[i]) * fex;
       Coef[2*npoints+i] = 4.0f * vsigma[i] * fex;

     } else { // CORRELATION

       Coef[i] += 2.0f * tred[i] * v2rho2[i] + 8.0f * cruz[i] * v2rhosigma[i];
       Coef[npoints+i] += 8.0f * tred[i] * v2rhosigma[i] + 32.0f * cruz[i] * v2sigma2[i];
       Coef[2*npoints+i] += 4.0f * vsigma[i];

     } // end if exc-corr

   } // end if points
}
#endif


// GPU VERSION
template <class T, int width>
void LibxcProxy_cuda <T, width>::coefLR (const int npoints, double* rho,
               double* sigma, double* red, double* cruz, double* lrCoef) // GPU
{
#ifdef __CUDACC__

   bool full_double = (sizeof(T) == 8); // LIBXC required double variables
   int size = sizeof(double) * npoints;
   cudaError_t err = cudaSuccess;

// Outputs libxc
   double* vrho       = NULL;
   double* vsigma     = NULL;
   double* v2rho2     = NULL;
   double* v2rhosigma = NULL;
   double* v2sigma2   = NULL;
   double* energy     = NULL;

// Allocate Outputs
   err = cudaMalloc((void**)&vrho,size);
   if (err != cudaSuccess) {
      printf("Error to allocate vrho");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&vsigma,size);
   if (err != cudaSuccess) {
      printf("Error to allocate vsigma");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v2rho2,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v2rho2");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v2rhosigma,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v2rhosigma");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v2sigma2,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v2sigma2");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&energy,size);
   if (err != cudaSuccess) {
      printf("Error to allocate energy");
      exit(EXIT_FAILURE);
   }

   cudaMemset(energy,0.0f,size);
   cudaMemset(vrho,0.0f,size);
   cudaMemset(vsigma,0.0f,size);
   cudaMemset(v2rho2,0.0f,size);
   cudaMemset(v2rhosigma,0.0f,size);
   cudaMemset(v2sigma2,0.0f,size);

// Call LIBXC for Exchange
   try {
       xc_gga_cuda(&funcForExchange,npoints,
               rho,
               sigma,
               energy,
               vrho,
               vsigma,
               v2rho2,
               v2rhosigma,
               v2sigma2,
               NULL, NULL, NULL, NULL);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }

   int threadsPerBlock = 256;
   int blocksPerGrid = (npoints + threadsPerBlock - 1) / threadsPerBlock;

// Obtain coef of Exchange
   local_coef_gpu<T,width> <<<blocksPerGrid,threadsPerBlock>>>(red,cruz,
   lrCoef,v2rho2,v2rhosigma,v2sigma2,vsigma,npoints,fact_exchange,true);

   cudaMemset(energy,0.0f,size);
   cudaMemset(vrho,0.0f,size);
   cudaMemset(vsigma,0.0f,size);
   cudaMemset(v2rho2,0.0f,size);
   cudaMemset(v2rhosigma,0.0f,size);
   cudaMemset(v2sigma2,0.0f,size);

// Call LIBXC for Coerrelation
   try {
       xc_gga_cuda(&funcForCorrelation,npoints,
               rho,
               sigma,
               energy,
               vrho,
               vsigma,
               v2rho2,
               v2rhosigma,
               v2sigma2,
               NULL, NULL, NULL, NULL);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }

// Obtain coef of Coerrelation
   local_coef_gpu<T,width> <<<blocksPerGrid,threadsPerBlock>>>(red,cruz,
   lrCoef,v2rho2,v2rhosigma,v2sigma2,vsigma,npoints,0.0f,false);

/*
   double* v2rho2_cpu = (double*) malloc(3*size);
   memset(v2rho2_cpu,0.0f,3*size);
   err = cudaMemcpy(v2rho2_cpu,lrCoef,3*size,cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) cout << "caca copy" << endl;
   for (int i=0; i<3*npoints; i++)
      printf("point(%d)=\t%f\t%f\t%f\n",i,v2rho2_cpu[i],v2rho2_cpu[npoints+i],v2rho2_cpu[2*npoints+i]);
      //printf("point(%d)=\t%f\n",i,v2rho2_cpu[i]);

   exit(-1);
*/

// Free memory
   cudaFree(vrho);       vrho = NULL;
   cudaFree(vsigma);     vsigma = NULL;
   cudaFree(v2rho2);     v2rho2 = NULL;
   cudaFree(v2rhosigma); v2rhosigma = NULL;
   cudaFree(v2sigma2);   v2sigma2 = NULL;
   cudaFree(energy);     energy = NULL;

#endif
   return;
}


// CPU VERSION
template <class T, int width>
void LibxcProxy_cuda <T, width>::coefLR (double *rho,
                double* sgm,
                double red,
                double cruz,
                double* lrCoef)
{
   // The otputs for exchange
   double vrhoX[2], vsigmaX[3],v2rho2X[3],v2rhosigmaX[6],v2sigma2X[6];

   // The ouputs for correlation
   double vrhoC[2], vsigmaC[3],v2rho2C[3],v2rhosigmaC[6],v2sigma2C[6];

   // NOT Refence
   double exc = 0.0f;

   // convert to alfa and beta;
   double dens[2], sigma[3];
          dens[0] = dens[1] = *rho * 0.5f;
          sigma[0] = sigma[1] = sigma[2] = *sgm * 0.25f;

   // Exchange values
   xc_gga_cuda(&funcForExchange,1,dens,sigma,&exc,vrhoX,vsigmaX,
          v2rho2X,v2rhosigmaX,v2sigma2X,NULL,NULL,NULL,NULL);

   v2rho2X[0]     *= fact_exchange;
   v2rhosigmaX[0] *= fact_exchange;
   v2sigma2X[0]   *= fact_exchange;
   vsigmaX[0]     *= fact_exchange;
   
   // Correlation values
   xc_gga_cuda(&funcForCorrelation,1,dens,sigma,&exc,vrhoC,vsigmaC,
          v2rho2C,v2rhosigmaC,v2sigma2C,NULL,NULL,NULL,NULL);

   // Results
   double term1, term2, term3;
   term1 = red * v2rho2X[0] + 2.0f * v2rhosigmaX[0] * cruz;
   term2 = red * v2rho2C[0] + 2.0f * v2rhosigmaC[0] * cruz;
         term2 += cruz * v2rhosigmaC[1];
   term3 = red * v2rho2C[1] + 2.0f * v2rhosigmaC[2] * cruz;
         term3 += cruz * v2rhosigmaC[1];
   lrCoef[0] = term1 + term2 + term3;

   term1 = red * v2rhosigmaX[0] * 2.0f + cruz * v2sigma2X[0] * 4.0f;
   term2 = red * (2.0f * v2rhosigmaC[0] + v2rhosigmaC[1]);
   term2 += cruz * (4.0f * v2sigma2C[0] + 2.0f * v2sigma2C[1]);
   term2 += cruz * (2.0f * v2sigma2C[1] + v2sigma2C[3]);
   term3 = red * (2.0f * v2rhosigmaC[3] + v2rhosigmaC[1]);
   term3 += cruz * (4.0f*v2sigma2C[2]+2.0f*v2sigma2C[1]+2.0f*v2sigma2C[4]+v2sigma2C[3]);
   lrCoef[1] = term1 + term2 + term3;
   
   term1 = 2.0f * vsigmaX[0];
   term2 = 2.0f * vsigmaC[0];
   term3 = vsigmaC[1];
   lrCoef[2] = term1 + term2 + term3;

   return;
}

#ifdef __CUDACC__

template <class T, int width>
__global__ void Zv_exchange(const int npoints,double* td,
           G2G::vec_type<T,WIDTH>* dxyz, G2G::vec_type<T,WIDTH>* txyz,
           double* Coef,double* v2rho2, double* v2rhosigma, double* v2sigma2,
           double* v3rho3,double* v3rho2sigma,double* v3rhosigma2,double* v3sigma3,double fact_ex)

{
   double fex = fact_ex;
   int i = blockDim.x * blockIdx.x + threadIdx.x;

   double DUMNV[2],DXV[2],DYV[2],DZV[2],DUMGRV[4],DUMXX[4];
   double C[10];

   if ( i < npoints ) {
      DUMNV[0] = DUMNV[1] = td[i];

      DUMGRV[0]=DUMGRV[1]=txyz[i].x*dxyz[i].x*0.5f + txyz[i].y*dxyz[i].y*0.5f + txyz[i].z*dxyz[i].z*0.5f;
      DUMGRV[2]=DUMGRV[3]=DUMGRV[0];

      DUMXX[0]=DUMXX[1]=txyz[i].x*txyz[i].x + txyz[i].y*txyz[i].y + txyz[i].z*txyz[i].z;
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
      XDUMA  = C[0] * 4.00f * v2rhosigma[i];
      XDUMAG = C[0] * 16.0f * v2sigma2[i];
      XDUMA += C[1] * 4.00f * v3rho3[i];
      XDUMAG+= C[1] * 16.0f * v3rho2sigma[i];
      XDUMA += C[2] * 8.00f * v3rho2sigma[i];
      XDUMAG+= C[2] * 32.0f * v3rhosigma2[i];
      XDUMA += C[3] * 8.00f * v3rho2sigma[i];
      XDUMAG+= C[3] * 32.0f * v3rhosigma2[i];
      XDUMA += C[6] * 16.0f * v3rhosigma2[i];
      XDUMAG+= C[6] * 64.0f * v3sigma3[i];

      double XDUMAGEA=0.0f;
      XDUMAGEA  = 4.0f * DUMNV[0]  * 4.0f * v2rhosigma[i];
      XDUMAGEA += 8.0f * DUMGRV[1] * 8.0f * v2sigma2[i];

      Coef[i]   = XDUMA * fex;
      Coef[npoints+i] = XDUMAG * fex;
      Coef[npoints*2+i] = XDUMAGEA * fex;
   }
}
#endif

#ifdef __CUDACC__

template <class T, int width>
__global__ void Zv_coulomb(const int npoints,double* td,
           G2G::vec_type<T,WIDTH>* dxyz, G2G::vec_type<T,WIDTH>* txyz,
           double* Coef, double* v2rhosigma, double* v2sigma2,
           double* v3rho3,double* v3rho2sigma,double* v3rhosigma2,double* v3sigma3)
{
   int i = blockDim.x * blockIdx.x + threadIdx.x;
   double DUMNV[2],DXV[2],DYV[2],DZV[2],DUMGRV[4],DUMXX[4];
   double C[20];


   if ( i < npoints ) {
      DUMNV[0] = DUMNV[1] = td[i];

      DUMGRV[0]=DUMGRV[1]=txyz[i].x*dxyz[i].x*0.5f + txyz[i].y*dxyz[i].y*0.5f + txyz[i].z*dxyz[i].z*0.5f;
      DUMGRV[2]=DUMGRV[3]=DUMGRV[0];

      DUMXX[0]=DUMXX[1]=txyz[i].x*txyz[i].x + txyz[i].y*txyz[i].y + txyz[i].z*txyz[i].z;
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
      CDUMA=C[0]*v2rhosigma[i];
      CDUMAG1=C[0]*2.0f*v2sigma2[i];
      CDUMAG2=C[0]*v2sigma2[i]*2.0f;
      CDUMA=CDUMA+C[1]*v3rho3[i];
      CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma[i];
      CDUMAG2=CDUMAG2+C[1]*v3rho2sigma[i]*2.0f;
      CDUMA=CDUMA+C[2]*v3rho2sigma[i];
      CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2[i];
      CDUMAG2=CDUMAG2+C[2]*v3rhosigma2[i]*2.0f;
      CDUMA=CDUMA+C[3]*v3rho2sigma[i];
      CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2[i];
      CDUMAG2=CDUMAG2+C[3]*v3rhosigma2[i]*2.0f;
      CDUMA=CDUMA+C[4]*v3rho2sigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[4]*v3rhosigma2[i]*4.0f;
      CDUMA=CDUMA+C[5]*v3rho2sigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[5]*v3rhosigma2[i]*4.0f;
      CDUMA=CDUMA+C[6]*v3rhosigma2[i];
      CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3[i];
      CDUMAG2=CDUMAG2+C[6]*v3sigma3[i]*2.0f;
      CDUMA=CDUMA+C[7]*v3rhosigma2[i]*2.0f;
      CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3[i]*2.0f;
      CDUMAG2=CDUMAG2+C[7]*v3sigma3[i]*4.0f;
      CDUMA=CDUMA+C[8]*v3rhosigma2[i]*2.0f;
      CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3[i]*2.0f;
      CDUMAG2=CDUMAG2+C[8]*v3sigma3[i]*4.0f;
      CDUMA=CDUMA+C[9]*v3rhosigma2[i]*4.0f;
      CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3[i]*4.0f;
      CDUMAG2=CDUMAG2+C[9]*v3sigma3[i]*8.0f;

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


      CDUMA=CDUMA+C[0]*v2rhosigma[i];
      CDUMAG1=CDUMAG1+C[0]*2.0f*v2sigma2[i];
      CDUMAG2=CDUMAG2+C[0]*v2sigma2[i]*2.0f;
      CDUMA=CDUMA+C[1]*v3rho3[i];
      CDUMAG1=CDUMAG1+C[1]*2.0f*v3rho2sigma[i];
      CDUMAG2=CDUMAG2+C[1]*v3rho2sigma[i]*2.0f;
      CDUMA=CDUMA+C[2]*v3rho2sigma[i];
      CDUMAG1=CDUMAG1+C[2]*2.0f*v3rhosigma2[i];
      CDUMAG2=CDUMAG2+C[2]*v3rhosigma2[i]*2.0f;
      CDUMA=CDUMA+C[3]*v3rho2sigma[i];
      CDUMAG1=CDUMAG1+C[3]*2.0f*v3rhosigma2[i];
      CDUMAG2=CDUMAG2+C[3]*v3rhosigma2[i]*2.0f;
      CDUMA=CDUMA+C[4]*v3rho2sigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[4]*2.0f*v3rhosigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[4]*v3rhosigma2[i]*4.0f;
      CDUMA=CDUMA+C[5]*v3rho2sigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[5]*2.0f*v3rhosigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[5]*v3rhosigma2[i]*4.0f;
      CDUMA=CDUMA+C[6]*v3rhosigma2[i];
      CDUMAG1=CDUMAG1+C[6]*2.0f*v3sigma3[i];
      CDUMAG2=CDUMAG2+C[6]*v3sigma3[i]*2.0f;
      CDUMA=CDUMA+C[7]*v3rhosigma2[i]*2.0f;
      CDUMAG1=CDUMAG1+C[7]*2.0f*v3sigma3[i]*2.0f;
      CDUMAG2=CDUMAG2+C[7]*v3sigma3[i]*4.0f;
      CDUMA=CDUMA+C[8]*v3rhosigma2[i]*2.0f;
      CDUMAG1=CDUMAG1+C[8]*2.0f*v3sigma3[i]*2.0f;
      CDUMAG2=CDUMAG2+C[8]*v3sigma3[i]*4.0f;
      CDUMA=CDUMA+C[9]*v3rhosigma2[i]*4.0f;
      CDUMAG1=CDUMAG1+C[9]*2.0f*v3sigma3[i]*4.0f;
      CDUMAG2=CDUMAG2+C[9]*v3sigma3[i]*8.0f;


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


      CDUMA=CDUMA+C[10]*v2rhosigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[10]*v2sigma2[i]*4.0f;
      CDUMA=CDUMA+C[11]*v3rho3[i];
      CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma[i];
      CDUMAG2=CDUMAG2+C[11]*v3rho2sigma[i]*2.0f;
      CDUMA=CDUMA+C[12]*v3rho2sigma[i];
      CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2[i];
      CDUMAG2=CDUMAG2+C[12]*v3rhosigma2[i]*2.0f;
      CDUMA=CDUMA+C[13]*v3rho2sigma[i];
      CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2[i];
      CDUMAG2=CDUMAG2+C[13]*v3rhosigma2[i]*2.0f;
      CDUMA=CDUMA+C[14]*v3rho2sigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[14]*v3rhosigma2[i]*4.0f;
      CDUMA=CDUMA+C[15]*v3rho2sigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[15]*v3rhosigma2[i]*4.0f;
      CDUMA=CDUMA+C[16]*v3rhosigma2[i];
      CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3[i];
      CDUMAG2=CDUMAG2+C[16]*v3sigma3[i]*2.0f;
      CDUMA=CDUMA+C[17]*v3rhosigma2[i]*2.0f;
      CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3[i]*2.0f;
      CDUMAG2=CDUMAG2+C[17]*v3sigma3[i]*4.0f;
      CDUMA=CDUMA+C[18]*v3rhosigma2[i]*2.0f;
      CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3[i]*2.0f;
      CDUMAG2=CDUMAG2+C[18]*v3sigma3[i]*4.0f;
      CDUMA=CDUMA+C[19]*v3rhosigma2[i]*4.0f;
      CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3[i]*4.0f;
      CDUMAG2=CDUMAG2+C[19]*v3sigma3[i]*8.0f;


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

      CDUMA=CDUMA+C[10]*v2rhosigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[10]*2.0f*v2sigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[10]*v2sigma2[i]*4.0f;
      CDUMA=CDUMA+C[11]*v3rho3[i];
      CDUMAG1=CDUMAG1+C[11]*2.0f*v3rho2sigma[i];
      CDUMAG2=CDUMAG2+C[11]*v3rho2sigma[i]*2.0f;
      CDUMA=CDUMA+C[12]*v3rho2sigma[i];
      CDUMAG1=CDUMAG1+C[12]*2.0f*v3rhosigma2[i];
      CDUMAG2=CDUMAG2+C[12]*v3rhosigma2[i]*2.0f;
      CDUMA=CDUMA+C[13]*v3rho2sigma[i];
      CDUMAG1=CDUMAG1+C[13]*2.0f*v3rhosigma2[i];
      CDUMAG2=CDUMAG2+C[13]*v3rhosigma2[i]*2.0f;
      CDUMA=CDUMA+C[14]*v3rho2sigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[14]*2.0f*v3rhosigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[14]*v3rhosigma2[i]*4.0f;
      CDUMA=CDUMA+C[15]*v3rho2sigma[i]*2.0f;
      CDUMAG1=CDUMAG1+C[15]*2.0f*v3rhosigma2[i]*2.0f;
      CDUMAG2=CDUMAG2+C[15]*v3rhosigma2[i]*4.0f;
      CDUMA=CDUMA+C[16]*v3rhosigma2[i];
      CDUMAG1=CDUMAG1+C[16]*2.0f*v3sigma3[i];
      CDUMAG2=CDUMAG2+C[16]*v3sigma3[i]*2.0f;
      CDUMA=CDUMA+C[17]*v3rhosigma2[i]*2.0f;
      CDUMAG1=CDUMAG1+C[17]*2.0f*v3sigma3[i]*2.0f;
      CDUMAG2=CDUMAG2+C[17]*v3sigma3[i]*4.0f;
      CDUMA=CDUMA+C[18]*v3rhosigma2[i]*2.0f;
      CDUMAG1=CDUMAG1+C[18]*2.0f*v3sigma3[i]*2.0f;
      CDUMAG2=CDUMAG2+C[18]*v3sigma3[i]*4.0f;
      CDUMA=CDUMA+C[19]*v3rhosigma2[i]*4.0f;
      CDUMAG1=CDUMAG1+C[19]*2.0f*v3sigma3[i]*4.0f;
      CDUMAG2=CDUMAG2+C[19]*v3sigma3[i]*8.0f;

      double CDUMAGEA=0.0f;
      CDUMAGEA=CDUMAGEA+4.0f*DUMNV[0]*v2rhosigma[i];
      CDUMAGEA=CDUMAGEA+8.0f*DUMGRV[0]*v2sigma2[i];
      CDUMAGEA=CDUMAGEA+4.0f*DUMGRV[2]*v2sigma2[i]*2.0f;
      CDUMAGEA=CDUMAGEA+4.0f*DUMNV[1]*v2rhosigma[i];
      CDUMAGEA=CDUMAGEA+8.0f*DUMGRV[1]*v2sigma2[i];
      CDUMAGEA=CDUMAGEA+4.0f*DUMGRV[3]*v2sigma2[i]*2.0f;

      double CDUMAGEC=0.0f;
      CDUMAGEC=CDUMAGEC+2.0f*DUMNV[0]*v2rhosigma[i]*2.0f;
      CDUMAGEC=CDUMAGEC+4.0f*DUMGRV[0]*v2sigma2[i]*2.0f;
      CDUMAGEC=CDUMAGEC+2.0f*DUMGRV[2]*v2sigma2[i]*4.0f;
      CDUMAGEC=CDUMAGEC+2.0f*DUMNV[1]*v2rhosigma[i]*2.0f;
      CDUMAGEC=CDUMAGEC+4.0f*DUMGRV[1]*v2sigma2[i]*2.0f;
      CDUMAGEC=CDUMAGEC+2.0f*DUMGRV[3]*v2sigma2[i]*4.0f;

      Coef[i]   += CDUMA;
      Coef[npoints+i] += CDUMAG1+CDUMAG2;
      Coef[npoints*2+i] += CDUMAGEA+CDUMAGEC;
   }
}
#endif

// GPU VERSION OF CoefZV
template <class T, int width>
void LibxcProxy_cuda <T, width>::coefZv(const int npoints, double* rho,
               double* sigma, double* red, G2G::vec_type<T,WIDTH>* dxyz, 
               G2G::vec_type<T,WIDTH>* txyz, double* lrCoef)
{
#ifdef __CUDACC__

   bool full_double = (sizeof(T) == 8); // LIBXC required double variables
   int size = sizeof(double) * npoints;
   cudaError_t err = cudaSuccess;

// Outputs libxc
   double* vrho        = NULL;
   double* vsigma      = NULL;
   double* v2rho2      = NULL;
   double* v2rhosigma  = NULL;
   double* v2sigma2    = NULL;
   double* v3rho3      = NULL;
   double* v3rho2sigma = NULL;
   double* v3rhosigma2 = NULL;
   double* v3sigma3    = NULL;
   double* energy      = NULL;

// Allocate Outputs
   err = cudaMalloc((void**)&vrho,size);
   if (err != cudaSuccess) {
      printf("Error to allocate vrho");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&vsigma,size);
   if (err != cudaSuccess) {
      printf("Error to allocate vsigma");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v2rho2,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v2rho2");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v2rhosigma,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v2rhosigma");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v2sigma2,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v2sigma2");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v3rho3,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v3rho3");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v3rho2sigma,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v3rho2sigma");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v3rhosigma2,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v3rhosigma2");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&v3sigma3,size);
   if (err != cudaSuccess) {
      printf("Error to allocate v3sigma3");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&energy,size);
   if (err != cudaSuccess) {
      printf("Error to allocate energy");
      exit(EXIT_FAILURE);
   }

   cudaMemset(energy,0.0f,size);
   cudaMemset(vrho,0.0f,size);
   cudaMemset(vsigma,0.0f,size);
   cudaMemset(v2rho2,0.0f,size);
   cudaMemset(v2rhosigma,0.0f,size);
   cudaMemset(v2sigma2,0.0f,size);
   cudaMemset(v3rho3,0.0f,size);
   cudaMemset(v3rho2sigma,0.0f,size);
   cudaMemset(v3rhosigma2,0.0f,size);
   cudaMemset(v3sigma3,0.0f,size);

// Call LIBXC for Exchange
   try {
       xc_gga_cuda(&funcForExchange,npoints,
               rho, sigma,
               energy,
               vrho, vsigma,
               v2rho2, v2rhosigma, v2sigma2,
               v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }

   int threadsPerBlock = 256;
   int blocksPerGrid = (npoints + threadsPerBlock - 1) / threadsPerBlock;

   Zv_exchange<T,width><<<blocksPerGrid,threadsPerBlock>>>(npoints,red,
                        dxyz,txyz,lrCoef,v2rho2,v2rhosigma,v2sigma2,
                        v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3,fact_exchange);

   cudaMemset(energy,0.0f,size);
   cudaMemset(vrho,0.0f,size);
   cudaMemset(vsigma,0.0f,size);
   cudaMemset(v2rho2,0.0f,size);
   cudaMemset(v2rhosigma,0.0f,size);
   cudaMemset(v2sigma2,0.0f,size);
   cudaMemset(v3rho3,0.0f,size);
   cudaMemset(v3rho2sigma,0.0f,size);
   cudaMemset(v3rhosigma2,0.0f,size);
   cudaMemset(v3sigma3,0.0f,size);

// Call LIBXC for Coerrelation
   try {
       xc_gga_cuda(&funcForCorrelation,npoints,
               rho,
               sigma,
               energy,
               vrho,
               vsigma,
               v2rho2,
               v2rhosigma,
               v2sigma2,
               v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }

   Zv_coulomb<T,width><<<blocksPerGrid,threadsPerBlock>>>(npoints,red,
                        dxyz,txyz,lrCoef,v2rhosigma,v2sigma2,
                        v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3);

   cudaFree(vrho);        vrho = NULL;
   cudaFree(vsigma);      vsigma = NULL;
   cudaFree(v2rho2);      v2rho2 = NULL;
   cudaFree(v2rhosigma);  v2rhosigma = NULL;
   cudaFree(v2sigma2);    v2sigma2 = NULL;
   cudaFree(v3rho3);      v3rho3 = NULL;
   cudaFree(v3rho2sigma); v3rho2sigma = NULL;
   cudaFree(v3rhosigma2); v3rhosigma2 = NULL;
   cudaFree(v3sigma3);    v3sigma3 = NULL;
   cudaFree(energy);      energy = NULL;

/*
   double* coef_cpu = (double*) malloc(3*size);
   double* red_cpu = (double*) malloc(size);
   memset(coef_cpu,0.0f,3*size);
   err = cudaMemcpy(coef_cpu,lrCoef,3*size,cudaMemcpyDeviceToHost);
   err = cudaMemcpy(red_cpu,red,size,cudaMemcpyDeviceToHost);
   if (err != cudaSuccess) cout << "caca copy" << endl;
   for (int i=0; i<npoints; i++)
      printf("%d\t%f\t%f\t%f\t%f\n",i,red_cpu[i],coef_cpu[3*i],coef_cpu[3*i+1],coef_cpu[3*i+2]);
      //printf("point(%d)=\t%f\n",i,v2rho2_cpu[i]);
   exit(-1);
*/


#endif
   return;
}

#ifdef __CUDACC__

template <class T, int width>
__global__ void save_derivs(const int npoints,const double fact_ex,const int xch,
           // INPUTS //
           const T* vrho_in, const T* vsigma_in, const T* v2rho2_in, const T* v2rhosigma_in, const T* v2sigma2_in,
           const T* v3rho3_in, const T* v3rho2sigma_in, const T* v3rhosigma2_in, const T* v3sigma3_in,
           // OUTPUTS //
           G2G::vec_type<T,2>* vrho, G2G::vec_type<T,2>* vsigma, G2G::vec_type<T,2>* v2rho2,
           G2G::vec_type<T,2>* v2rhosigma, G2G::vec_type<T,2>* v2sigma2, G2G::vec_type<T,2>* v3rho3,
           G2G::vec_type<T,2>* v3rho2sigma, G2G::vec_type<T,2>* v3rhosigma2, G2G::vec_type<T,2>* v3sigma3)
{
 
   const double fex = fact_ex;
   int i = blockDim.x * blockIdx.x + threadIdx.x;

   if ( i < npoints ) {
      if ( xch == 0 ) {
         vrho[i].x = vrho_in[i] * fex; vsigma[i].x = vsigma_in[i] * fex;
         v2rho2[i].x = v2rho2_in[i] * fex; v2rhosigma[i].x = v2rhosigma_in[i] * fex; 
         v2sigma2[i].x = v2sigma2_in[i] * fex;
         v3rho3[i].x = v3rho3_in[i] * fex; v3rho2sigma[i].x = v3rho2sigma_in[i] * fex; 
         v3rhosigma2[i].x = v3rhosigma2_in[i] * fex; v3sigma3[i].x = v3sigma3_in[i] * fex;
      } else {
         vrho[i].y = vrho_in[i]; vsigma[i].y = vsigma_in[i];
         v2rho2[i].y = v2rho2_in[i]; v2rhosigma[i].y = v2rhosigma_in[i];
         v2sigma2[i].y = v2sigma2_in[i];
         v3rho3[i].y = v3rho3_in[i]; v3rho2sigma[i].y = v3rho2sigma_in[i]; 
         v3rhosigma2[i].y = v3rhosigma2_in[i]; v3sigma3[i].y = v3sigma3_in[i];
      }
   }
}

#endif


template <class T, int width>
void LibxcProxy_cuda<T,width>::obtain_gpu_der(const int npoints,
     T* rho, T* sigma,
     // outputs //
     G2G::vec_type<T,2>* vrho,
     G2G::vec_type<T,2>* vsigma,
     G2G::vec_type<T,2>* v2rho2,
     G2G::vec_type<T,2>* v2rhosigma,
     G2G::vec_type<T,2>* v2sigma2,
     G2G::vec_type<T,2>* v3rho3,
     G2G::vec_type<T,2>* v3rho2sigma,
     G2G::vec_type<T,2>* v3rhosigma2,
     G2G::vec_type<T,2>* v3sigma3)
{
#ifdef __CUDACC__
   int size = sizeof(double) * npoints;
   cudaError_t err = cudaSuccess;

// Local libxc Variables
   double* lvrho        = NULL;
   double* lvsigma      = NULL;
   double* lv2rho2      = NULL;
   double* lv2rhosigma  = NULL;
   double* lv2sigma2    = NULL;
   double* lv3rho3      = NULL;
   double* lv3rho2sigma = NULL;
   double* lv3rhosigma2 = NULL;
   double* lv3sigma3    = NULL;
   double* energy       = NULL;

// Allocate Outputs
   err = cudaMalloc((void**)&lvrho,size);
   if (err != cudaSuccess) {
      printf("Error to allocate lvrho");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&lvsigma,size);
   if (err != cudaSuccess) {
      printf("Error to allocate lvsigma");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&lv2rho2,size);
   if (err != cudaSuccess) {
      printf("Error to allocate lv2rho2");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&lv2rhosigma,size);
   if (err != cudaSuccess) {
      printf("Error to allocate lv2rhosigma");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&lv2sigma2,size);
   if (err != cudaSuccess) {
      printf("Error to allocate vl2sigma2");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&lv3rho3,size);
   if (err != cudaSuccess) {
      printf("Error to allocate lv3rho3");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&lv3rho2sigma,size);
   if (err != cudaSuccess) {
      printf("Error to allocate vl3rho2sigma");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&lv3rhosigma2,size);
   if (err != cudaSuccess) {
      printf("Error to allocate lv3rhosigma2");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&lv3sigma3,size);
   if (err != cudaSuccess) {
      printf("Error to allocate lv3sigma3");
      exit(EXIT_FAILURE);
   }

   err = cudaMalloc((void**)&energy,size);
   if (err != cudaSuccess) {
      printf("Error to allocate energy");
      exit(EXIT_FAILURE);
   }

   cudaMemset(energy,0.0f,size);
   cudaMemset(lvrho,0.0f,size);
   cudaMemset(lvsigma,0.0f,size);
   cudaMemset(lv2rho2,0.0f,size);
   cudaMemset(lv2rhosigma,0.0f,size);
   cudaMemset(lv2sigma2,0.0f,size);
   cudaMemset(lv3rho3,0.0f,size);
   cudaMemset(lv3rho2sigma,0.0f,size);
   cudaMemset(lv3rhosigma2,0.0f,size);
   cudaMemset(lv3sigma3,0.0f,size);

// Call LIBXC for Exchange
   try {
       xc_gga_cuda(&funcForExchange,npoints,
               rho, sigma,
               energy,
               lvrho, lvsigma,
               lv2rho2, lv2rhosigma, lv2sigma2,
               lv3rho3, lv3rho2sigma, lv3rhosigma2, lv3sigma3);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }

   int threadsPerBlock = 256;
   int blocksPerGrid = (npoints + threadsPerBlock - 1) / threadsPerBlock;

   save_derivs<T,width><<<blocksPerGrid,threadsPerBlock>>>(npoints,fact_exchange,0,
              // INPUTS //
              lvrho, lvsigma, lv2rho2, lv2rhosigma, lv2sigma2, lv3rho3,
              lv3rho2sigma, lv3rhosigma2, lv3sigma3,
              // OUTPUTS //
              vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, v3rho3,
              v3rho2sigma, v3rhosigma2, v3sigma3);

   cudaMemset(energy,0.0f,size);
   cudaMemset(lvrho,0.0f,size);
   cudaMemset(lvsigma,0.0f,size);
   cudaMemset(lv2rho2,0.0f,size);
   cudaMemset(lv2rhosigma,0.0f,size);
   cudaMemset(lv2sigma2,0.0f,size);
   cudaMemset(lv3rho3,0.0f,size);
   cudaMemset(lv3rho2sigma,0.0f,size);
   cudaMemset(lv3rhosigma2,0.0f,size);
   cudaMemset(lv3sigma3,0.0f,size);

// Call LIBXC for Coerrelation
   try {
       xc_gga_cuda(&funcForCorrelation,npoints,
               rho, sigma, energy,
               lvrho, lvsigma, lv2rho2, lv2rhosigma, lv2sigma2, 
               lv3rho3, lv3rho2sigma, lv3rhosigma2, lv3sigma3);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }

   save_derivs<T,width><<<blocksPerGrid,threadsPerBlock>>>(npoints,0.0f,1,
              // INPUTS //
              lvrho, lvsigma, lv2rho2, lv2rhosigma, lv2sigma2, lv3rho3,
              lv3rho2sigma, lv3rhosigma2, lv3sigma3,
              // OUTPUTS //
              vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, v3rho3,
              v3rho2sigma, v3rhosigma2, v3sigma3);

   cudaFree(lvrho);        lvrho = NULL;
   cudaFree(lvsigma);      lvsigma = NULL;
   cudaFree(lv2rho2);      lv2rho2 = NULL;
   cudaFree(lv2rhosigma);  lv2rhosigma = NULL;
   cudaFree(lv2sigma2);    lv2sigma2 = NULL;
   cudaFree(lv3rho3);      lv3rho3 = NULL;
   cudaFree(lv3rho2sigma); lv3rho2sigma = NULL;
   cudaFree(lv3rhosigma2); lv3rhosigma2 = NULL;
   cudaFree(lv3sigma3);    lv3sigma3 = NULL;
   cudaFree(energy);       energy = NULL;
#endif

}


template <class T, int width>
void LibxcProxy_cuda <T, width>::doGGA (T rho,
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
    xc_gga_fxc_cuda (&funcForExchange,1,
                &dens,
                &sig,
                &v2rho2X,
                &v2rhosigmaX,
                &v2sigma2X);

    // Correlation values
    xc_gga_fxc_cuda (&funcForCorrelation,1,
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
// LibxcProxy_cuda::joinResults
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
/*
	vrho[i] += vrhoC[i];
        vsigma[i] += vsigmaC[i];
        v2rho[i] += v2rhoC[i];
        v2rhosigma[i] += v2rhosigmaC[i];
        v2sigma[i] += v2sigmaC[i];
*/
	vrho[i]       = fEE * vrho[i] + vrhoC[i];
        vsigma[i]     = fEE * vsigma[i] + vsigmaC[i];
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
__forceinline__ __global__ void convertFloatToDouble(const float* input, double* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	output[i] = (double)input[i];
    }
}

__forceinline__ __global__ void convertDoubleToFloat(const double* input, float* output, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements)
    {
	output[i] = (float)input[i];
    }
}

__forceinline__ __global__ void convertFloatToDouble(const G2G::vec_type<float,4>* input, G2G::vec_type<double,4>* output, int numElements)
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
__forceinline__ __global__ void convertDoubleToFloat(const double* input, double* output, int numElements)
{
    return;
}

// Esto es para engañar al compilador pq el FLAG FULL_DOUBLE
// a veces permite que T sea double y se rompe todo.
__forceinline__ __global__ void convertFloatToDouble(const double* input, double* output, int numElements)
{
    return;
}

// Esto es para engañar al compilador pq el FLAG FULL_DOUBLE
// a veces permite que T sea double y se rompe todo.
__forceinline__ __global__ void convertFloatToDouble(const G2G::vec_type<double,4>* input, G2G::vec_type<double,4>* output, int numElements)
{
    return;
}

__forceinline__ __global__ void convertDoubleToFloat(const G2G::vec_type<double,4>* input, G2G::vec_type<float,4>* output, int numElements)
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// LibxcProxy_cuda::doGGA - GPU version
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
void LibxcProxy_cuda <T, width>::doGGA(T* dens,
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
    //printf("doGGA - GPU, exact %f\n",fact_exchange);
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
        fprintf(stderr, "Failed to allocate device rho si!\n");
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
        xc_gga_cuda (&funcForExchange, number_of_points,
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
        xc_gga_cuda (&funcForCorrelation, number_of_points,
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
// LibxcProxy_cuda::doLDA - GPU version
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
void LibxcProxy_cuda <T, width>::doLDA(T dens, const G2G::vec_type<T, width> &grad, const G2G::vec_type<T, width> &hess1, const G2G::vec_type<T, width> &hess2, T &ex, T &ec, T &y2a)
{
    //TODO: not implemented yet!
    return;
}

#endif // LIBXCPROXY_CUDA_H
