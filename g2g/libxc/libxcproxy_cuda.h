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
    LibxcProxy_cuda (int* func_id,double* func_coef, int nx_coef, int nc_coef,
                                int nsr_id, double screen, int nSpin);

    void init_cuda (int, int, int);

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
LibxcProxy_cuda <T, width>::LibxcProxy_cuda (int* func_id,double* func_coef, int nx_coef, int nc_coef,
                                   int nsr_id, double screen, int nSpin)
{
// Checkeamos que el funcional sea pbe o pbe0, ya que libxc en GPU solo funciona con estos 
// funcionales
    int ntot = nx_coef + nc_coef;
    if ( ntot != 2 ) {
       cout << "Libxc gpu implementation only works with pbe or pbe0" << endl;
       cout << "If you want other functinal, please compile LIO in cpu with libxc=1" << endl;
       exit(-1);
    }
    int other_func = 0;
    if ( func_id[0] != 101 ) other_func = 1;
    if ( func_id[1] != 130 ) other_func = 1;

    if ( other_func == 1 ) {
       cout << "Libxc gpu implementation only works with pbe or pbe0" << endl;
       cout << "If you want other functinal, please compile LIO in cpu with libxc=1" << endl;
       exit(-1);
    }
    fact_exchange = func_coef[0];
    int func_x, func_c;
    func_x = func_id[0] + 1000;
    func_c = func_id[1] + 1000;
    init_cuda (func_x, func_c, nSpin);
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
   // printf("LibxcProxy_cuda::closeProxy()\n");
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
__global__ void local_coef_gpu(
           double* tred, double* cruz, const int npoints, double fact_ex, // INPUTS
           double* Coef, // OUTPUT
           double* xv2rho2, double* xv2rhosigma,double* xv2sigma2,double* xvsigma, // Exchange
           double* cv2rho2, double* cv2rhosigma,double* cv2sigma2,double* cvsigma) // Correlation
{
   double fex = fact_ex;
   int i = blockDim.x * blockIdx.x + threadIdx.x;

   if ( i < npoints ) {
      Coef[i] = tred[i] * (fex*xv2rho2[i]+cv2rho2[i]) + 2.0f * cruz[i] * (fex*xv2rhosigma[i]+cv2rhosigma[i]);
      Coef[npoints+i] = 2.0f * tred[i] * (fex*xv2rhosigma[i]+cv2rhosigma[i]) + 4.0f * cruz[i] * (fex*xv2sigma2[i]+cv2sigma2[i]);
      Coef[2*npoints+i] = 2.0f * (fex*xvsigma[i]+cvsigma[i]);
   } // end if points
}
#endif


// GPU VERSION
template <class T, int width>
void LibxcProxy_cuda <T, width>::coefLR (const int npoints, double* rho,
               double* sigma, double* red, double* cruz, double* lrCoef) // GPU
{
#ifdef __CUDACC__

   int size = sizeof(double) * npoints;
   cudaError_t err = cudaSuccess;

// Outputs for each functional: Exchange and Correlation
   double* xvrho       = NULL;
   double* xvsigma     = NULL;
   double* xv2rho2     = NULL;
   double* xv2rhosigma = NULL;
   double* xv2sigma2   = NULL;
   double* xenergy     = NULL;
   double* cvrho       = NULL;
   double* cvsigma     = NULL;
   double* cv2rho2     = NULL;
   double* cv2rhosigma = NULL;
   double* cv2sigma2   = NULL;
   double* cenergy     = NULL;

// Allocate Exchange Outputs
   err = cudaMalloc((void**)&xvrho,size);
   err = cudaMalloc((void**)&xvsigma,size);
   err = cudaMalloc((void**)&xv2rho2,size);
   err = cudaMalloc((void**)&xv2rhosigma,size);
   err = cudaMalloc((void**)&xv2sigma2,size);
   err = cudaMalloc((void**)&xenergy,size);

// Allocate Correlation Outputs
   err = cudaMalloc((void**)&cvrho,size);
   err = cudaMalloc((void**)&cvsigma,size);
   err = cudaMalloc((void**)&cv2rho2,size);
   err = cudaMalloc((void**)&cv2rhosigma,size);
   err = cudaMalloc((void**)&cv2sigma2,size);
   err = cudaMalloc((void**)&cenergy,size);

   // Set zero exchange outputs
   cudaMemset(xvrho,0.0f,size);
   cudaMemset(xvsigma,0.0f,size);
   cudaMemset(xv2rho2,0.0f,size);
   cudaMemset(xv2rhosigma,0.0f,size);
   cudaMemset(xv2sigma2,0.0f,size);
   cudaMemset(xenergy,0.0f,size);

   // Set zero correlation outputs
   cudaMemset(cvrho,0.0f,size);
   cudaMemset(cvsigma,0.0f,size);
   cudaMemset(cv2rho2,0.0f,size);
   cudaMemset(cv2rhosigma,0.0f,size);
   cudaMemset(cv2sigma2,0.0f,size);
   cudaMemset(cenergy,0.0f,size);

// Call LIBXC for Exchange
   try {
       xc_gga_cuda(&funcForExchange,npoints,
               rho,
               sigma,
               xenergy,
               xvrho,
               xvsigma,
               xv2rho2,
               xv2rhosigma,
               xv2sigma2,
               NULL, NULL, NULL, NULL);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }

// Call LIBXC for Correlation
   try {
       xc_gga_cuda(&funcForCorrelation,npoints,
               rho,
               sigma,
               cenergy,
               cvrho,
               cvsigma,
               cv2rho2,
               cv2rhosigma,
               cv2sigma2,
               NULL, NULL, NULL, NULL);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }

   int threadsPerBlock = 256;
   int blocksPerGrid = (npoints + threadsPerBlock - 1) / threadsPerBlock;

// Obtain coef LR
   local_coef_gpu<T,width> <<<blocksPerGrid,threadsPerBlock>>>(
                 red,cruz,npoints,fact_exchange, // INPUTS
                 lrCoef, // OUTPUT
                 xv2rho2,xv2rhosigma,xv2sigma2,xvsigma, // Exchange
                 cv2rho2,cv2rhosigma,cv2sigma2,cvsigma); // Correlation

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
   cudaFree(xvrho);       xvrho = NULL;
   cudaFree(xvsigma);     xvsigma = NULL;
   cudaFree(xv2rho2);     xv2rho2 = NULL;
   cudaFree(xv2rhosigma); xv2rhosigma = NULL;
   cudaFree(xv2sigma2);   xv2sigma2 = NULL;
   cudaFree(xenergy);     xenergy = NULL;
   cudaFree(cvrho);       cvrho = NULL;
   cudaFree(cvsigma);     cvsigma = NULL;
   cudaFree(cv2rho2);     cv2rho2 = NULL;
   cudaFree(cv2rhosigma); cv2rhosigma = NULL;
   cudaFree(cv2sigma2);   cv2sigma2 = NULL;
   cudaFree(cenergy);     cenergy = NULL;

#endif
   return;
}

#ifdef __CUDACC__

template <class T, int width>
__global__ void Zv_factor(
                const int npoints,double* td, G2G::vec_type<T,WIDTH>* dxyz, G2G::vec_type<T,WIDTH>* txyz, // inputs var.

                double* xv3rho3,double* xv3rho2sigma,double* xv3rhosigma2,
                double* xv2rhosigma, double* xv3sigma3, double* xv2sigma2, // inputs exchange

                double* cv3rho3,double* cv3rho2sigma,double* cv3rhosigma2,
                double* cv2rhosigma, double* cv3sigma3, double* cv2sigma2, // inputs correlation

                double fact_ex, double* Coef) // output
{
   double fex = fact_ex;
   int ii = blockDim.x * blockIdx.x + threadIdx.x;
   bool valid_thread  =  ( ii < npoints );

   double v3rho3, v3rho2sigma, v3rhosigma2, v2rhosigma, v3sigma3, v2sigma2;
   double red, prod_rr, prod_pr;
   double term1, term2, term3, term4;
  
   if ( ii < npoints ) {
      // generate inputs
      red     = td[ii];
      prod_rr = txyz[ii].x*txyz[ii].x + txyz[ii].y*txyz[ii].y + txyz[ii].z*txyz[ii].z;
      prod_pr = dxyz[ii].x*txyz[ii].x + dxyz[ii].y*txyz[ii].y + dxyz[ii].z*txyz[ii].z;

      // total functional derivative
      v3rho3      = fex * xv3rho3[ii] + cv3rho3[ii];
      v3rho2sigma = fex * xv3rho2sigma[ii] + cv3rho2sigma[ii];
      v3rhosigma2 = fex * xv3rhosigma2[ii] + cv3rhosigma2[ii];
      v2rhosigma  = fex * xv2rhosigma[ii] + cv2rhosigma[ii];
      v3sigma3    = fex * xv3sigma3[ii] + cv3sigma3[ii];
      v2sigma2    = fex * xv2sigma2[ii] + cv2sigma2[ii];

      // Coef 0
      term1 = red * red * v3rho3;
      term2 = 4.0f * red * prod_pr * v3rho2sigma;
      term3 = 4.0f * prod_pr * prod_pr * v3rhosigma2;
      term4 = 2.0f * prod_rr * v2rhosigma;
      Coef[ii] = (term1 + term2 + term3 + term4) * 2.0f;

      // Coef 1
      term1 = 2.0f * red * red * v3rho2sigma;
      term2 = 8.0f * red * prod_pr * v3rhosigma2;
      term3 = 8.0f * prod_pr * prod_pr * v3sigma3;
      term4 = 4.0f * prod_rr * v2sigma2;
      Coef[npoints+ii] = (term1 + term2 + term3 + term4) * 2.0f;

      // Coef 2
      term1 = 4.0f * red * v2rhosigma + 8.0f * prod_pr * v2sigma2;
      Coef[npoints*2+ii] = term1 * 2.0f;
   } // end valid thread
}
#endif

// GPU VERSION OF CoefZV
template <class T, int width>
void LibxcProxy_cuda <T, width>::coefZv(const int npoints, double* rho,
               double* sigma, double* red, G2G::vec_type<T,WIDTH>* dxyz, 
               G2G::vec_type<T,WIDTH>* txyz, double* lrCoef)
{
#ifdef __CUDACC__

   int size = sizeof(double) * npoints;
   cudaError_t err = cudaSuccess;

// Outputs exchange libxc
   double* xv3rho3      = NULL;
   double* xv3rho2sigma = NULL;
   double* xv3rhosigma2 = NULL;
   double* xv2rhosigma  = NULL;
   double* xv3sigma3    = NULL;
   double* xv2sigma2    = NULL;

// Outputs correlation libxc
   double* cv3rho3      = NULL;
   double* cv3rho2sigma = NULL;
   double* cv3rhosigma2 = NULL;
   double* cv2rhosigma  = NULL;
   double* cv3sigma3    = NULL;
   double* cv2sigma2    = NULL;

// Unusefull Variables
   double* uvrho        = NULL;
   double* uvsigma      = NULL;
   double* uv2rho2      = NULL;
   double* uenergy      = NULL;

// Allocate Exchange outputs
   err = cudaMalloc((void**)&xv3rho3,size);
   err = cudaMalloc((void**)&xv3rho2sigma,size);
   err = cudaMalloc((void**)&xv3rhosigma2,size);
   err = cudaMalloc((void**)&xv2rhosigma,size);
   err = cudaMalloc((void**)&xv3sigma3,size);
   err = cudaMalloc((void**)&xv2sigma2,size);

// Allocate Correlation outputs
   err = cudaMalloc((void**)&cv3rho3,size);
   err = cudaMalloc((void**)&cv3rho2sigma,size);
   err = cudaMalloc((void**)&cv3rhosigma2,size);
   err = cudaMalloc((void**)&cv2rhosigma,size);
   err = cudaMalloc((void**)&cv3sigma3,size);
   err = cudaMalloc((void**)&cv2sigma2,size);

// Allocate Unusefull outputs
   err = cudaMalloc((void**)&uvrho,size);
   err = cudaMalloc((void**)&uvsigma,size);
   err = cudaMalloc((void**)&uv2rho2,size);
   err = cudaMalloc((void**)&uenergy,size);

// Set Zero Exchange ouputs
   cudaMemset(xv3rho3,0.0f,size);
   cudaMemset(xv3rho2sigma,0.0f,size);
   cudaMemset(xv3rhosigma2,0.0f,size);
   cudaMemset(xv2rhosigma,0.0f,size);
   cudaMemset(xv3sigma3,0.0f,size);
   cudaMemset(xv2sigma2,0.0f,size);

// Set Zero Correlation ouputs
   cudaMemset(cv3rho3,0.0f,size);
   cudaMemset(cv3rho2sigma,0.0f,size);
   cudaMemset(cv3rhosigma2,0.0f,size);
   cudaMemset(cv2rhosigma,0.0f,size);
   cudaMemset(cv2sigma2,0.0f,size);
   cudaMemset(cv3sigma3,0.0f,size);

// Set Zero Unusefull outputs
   cudaMemset(uenergy,0.0f,size);
   cudaMemset(uvrho,0.0f,size);
   cudaMemset(uvsigma,0.0f,size);
   cudaMemset(uv2rho2,0.0f,size);

// Call LIBXC for Exchange
   try {
       xc_gga_cuda(&funcForExchange,npoints,
               rho, sigma,
               uenergy,
               uvrho, uvsigma,
               uv2rho2, xv2rhosigma, xv2sigma2,
               xv3rho3, xv3rho2sigma, xv3rhosigma2, xv3sigma3);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }

// Call LIBXC for Coerrelation
   try {
       xc_gga_cuda(&funcForCorrelation,npoints,
               rho, sigma,
               uenergy,
               uvrho, uvsigma,
               uv2rho2, cv2rhosigma, cv2sigma2,
               cv3rho3, cv3rho2sigma, cv3rhosigma2, cv3sigma3);
   } catch (int exception) {
       fprintf (stderr, "Exception ocurred calling xc_gga for Exchange '%d' \n", exception);
       return;
   }
   // Free unusefull outputs
   cudaFree(uvrho);   uvrho = NULL;
   cudaFree(uvsigma); uvsigma = NULL;
   cudaFree(uv2rho2); uv2rho2 = NULL;
   cudaFree(uenergy); uenergy = NULL;

   int threadsPerBlock = 256;
   int blocksPerGrid = (npoints + threadsPerBlock - 1) / threadsPerBlock;
   Zv_factor<T,width><<<blocksPerGrid,threadsPerBlock>>>(
        npoints,red,dxyz,txyz,
        xv3rho3,xv3rho2sigma,xv3rhosigma2,xv2rhosigma,xv3sigma3,xv2sigma2,
        cv3rho3,cv3rho2sigma,cv3rhosigma2,cv2rhosigma,cv3sigma3,cv2sigma2,
        fact_exchange,lrCoef);

// Free all memory
   cudaFree(xv3rho3);      xv3rho3 = NULL;
   cudaFree(xv3rho2sigma); xv3rho2sigma = NULL;
   cudaFree(xv3rhosigma2); xv3rhosigma2 = NULL;
   cudaFree(xv2rhosigma);  xv2rhosigma = NULL;
   cudaFree(xv2sigma2);    xv2sigma2 = NULL;
   cudaFree(xv3sigma3);    xv3sigma3 = NULL;
   cudaFree(cv3rho3);      cv3rho3 = NULL;
   cudaFree(cv3rho2sigma); cv3rho2sigma = NULL;
   cudaFree(cv3rhosigma2); cv3rhosigma2 = NULL;
   cudaFree(cv2rhosigma);  cv2rhosigma = NULL;
   cudaFree(cv2sigma2);    cv2sigma2 = NULL;
   cudaFree(cv3sigma3);    cv3sigma3 = NULL;
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
