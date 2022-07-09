#include <iostream>
#include <omp.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"

#include <stdio.h>
#include <string.h>

using namespace G2G;
extern Partition partition;

void calc_gradients(double* dens, double* diff, double* trad,
                    double sigma, double* dfac, double* pfac,
                    double* tfac, int calc_fxc);


extern "C" void g2g_calcgradxc_(double* P,double* V, double* F, int& met)
{
   partition.solveForcesExc(P,V,F,met);
}

namespace G2G {

void Partition::solveForcesExc(double* P,double* V,double* F,int met)
{
    std::vector< HostMatrix<double> > forces;
    forces.resize(G2G::cpu_threads + G2G::gpu_threads);

#pragma omp parallel for num_threads(cpu_threads + gpu_threads) schedule(static)
    for(uint i=0;i<work.size();i++) {
#if GPU_KERNELS
       bool gpu_thread = false;
       if (i >= cpu_threads) {
          gpu_thread = true;
          cudaSetDevice(i - cpu_threads);
       }
#endif

      forces[i].resize(fortran_vars.max_atoms, 3); // max_atoms = atoms (atomos qm)
      forces[i].zero();

      for(uint j=0;j<work[i].size();j++) {
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind - cubes.size()]->solve_for_exc(P,V,forces[i],met);
         } else {
           cubes[ind]->solve_for_exc(P,V,forces[i],met);
         }
#if GPU_KERNELS
         if (gpu_thread) {
           cudaDeviceSynchronize();
         }
#endif
      }
   }

   for(int k=0; k<forces.size(); k++) {
      for(int i=0; i<fortran_vars.atoms; i++) {
         F[i*3]   += forces[k](i,0);
         F[i*3+1] += forces[k](i,1);
         F[i*3+2] += forces[k](i,2);
      }
   }
   std::vector<HostMatrix<double> >().swap(forces);
   fflush(stdout);
}

template<class scalar_type> void PointGroupCPU<scalar_type>::
               solve_for_exc(double*P,double*V,HostMatrix<double>& F,int MET)
{
   const uint group_m = this->total_functions();
   const int npoints = this->points.size();
   bool lda = false;
   bool compute_forces = true;

   compute_functions(compute_forces,!lda);

   int M = fortran_vars.m;
   int natom = fortran_vars.atoms;
   int local_atoms = this->total_nucleii();
   double* dens = (double*) malloc(4*sizeof(double));
   double* trad = (double*) malloc(4*sizeof(double));
   double* diff = (double*) malloc(4*sizeof(double));
   double* dfac = (double*) malloc(4*sizeof(double));
   double* pfac = (double*) malloc(4*sizeof(double));
   double* tfac = (double*) malloc(4*sizeof(double));

// Force for this group
   int dim = 3*local_atoms;
   double* sForce = (double*) malloc(dim*sizeof(double));
   memset(sForce,0.0f,dim*sizeof(double));

// Reduced Matrix for this group
   HostMatrix<scalar_type> rmm_input(group_m,group_m);
   HostMatrix<scalar_type> Pred(group_m,group_m);
   HostMatrix<scalar_type> Vred(group_m,group_m);
   HostMatrix<scalar_type> Pbig(M*(M+1)/2);
   HostMatrix<scalar_type> Vbig(M*(M+1)/2);

   get_rmm_input(rmm_input);
   int index = 0;
   for(int row=0; row<M; row++) {
     Pbig(index) = P[row*M+row];
     Vbig(index) = V[row*M+row];
     index += 1;
     for(int col=row+1; col<M; col++) {
       Pbig(index) = P[row*M+col] + P[col*M+row];
       Vbig(index) = V[row*M+col] + V[col*M+row];
       index += 1;
     }
   }
   get_tred_input(Pred,Pbig);
   get_tred_input(Vred,Vbig);
   Pbig.deallocate();
   Vbig.deallocate();

   // Gradients for nuclei
   HostMatrix<scalar_type> ddx, ddy, ddz;
   ddx.resize(local_atoms, 1);
   ddy.resize(local_atoms, 1);
   ddz.resize(local_atoms, 1);

   HostMatrix<scalar_type> groundD(4);
   HostMatrix<scalar_type> pdiffD(4);
   HostMatrix<scalar_type> transD(4);

   for(int point=0;point<npoints;point++) {
    double pd, pdx, pdy, pdz; pd = pdx = pdy = pdz = 0.0f;
    double pp, ppx, ppy, ppz; pp = ppx = ppy = ppz = 0.0f;
    double pt, ptx, pty, ptz; pt = ptx = pty = ptz = 0.0f;
    // functions and gradients values
    const scalar_type* fv = function_values.row(point);
    const scalar_type* gfx = gX.row(point);
    const scalar_type* gfy = gY.row(point);
    const scalar_type* gfz = gZ.row(point);
    // Hessianos de las funciones bases
    const scalar_type* hxx = hPX.row(point); //XX
    const scalar_type* hyy = hPY.row(point); //YY
    const scalar_type* hzz = hPZ.row(point); //ZZ
    const scalar_type* hxy = hIX.row(point); //XY
    const scalar_type* hxz = hIY.row(point); //XZ
    const scalar_type* hyz = hIZ.row(point); //YZ

    // Calculate GS and transition densities and derivatives in the point
    #define recalc_params \
    function_values.row(point), gX.row(point), gY.row(point), gZ.row(point), \
    rmm_input, Pred, Vred, point, \
    groundD, pdiffD, transD
    recalc_densGS3(recalc_params);
    #undef recalc_params

    // Copying outputs
    pd = groundD(0); pdx = groundD(1); pdy = groundD(2); pdz = groundD(3);
    pp = pdiffD(0);  ppx = pdiffD(1);  ppy = pdiffD(2);  ppz = pdiffD(3);
    pt = transD(0);  ptx = transD(1);  pty = transD(2);  ptz = transD(3);

    double sigma = (pdx * pdx) + (pdy * pdy) + (pdz * pdz);
    dens[0] = pd; dens[1] = pdx; dens[2] = pdy; dens[3] = pdz;
    diff[0] = pp; diff[1] = ppx; diff[2] = ppy; diff[3] = ppz;
    trad[0] = pt; trad[1] = ptx; trad[2] = pty; trad[3] = ptz;

    // Terms Calculate
    calc_gradients(dens,diff,trad,sigma,dfac,pfac,tfac,MET);

    // FORCES CALCULATE
    double DJII, PJII, VJII;
    double temp[4];
    // Gradients for the nuclei
    ddx.zero();
    ddy.zero();
    ddz.zero();

    const scalar_type wp = this->points[point].weight;

    for (int i = 0, ii = 0; i < this->total_functions_simple(); i++) {// cantidad de bases, s+p+d (sin tener en cuenta las 3p y 5d
       uint nuc = this->func2local_nuc(ii);// da a que nucleo LOCAL pertenece la base ii
       uint inc_i = this->small_function_type(i);// da cuantas funciones tiene ese tipo: s->1, p->3, d->5
       double grdx, grdy, grdz; grdx = grdy = grdz = 0.0f;

       for (uint k = 0; k < inc_i; k++, ii++) {
          for (uint j = 0; j < group_m; j++) { 
             double factor = (ii == j ? 2.0f : 1.0f);
             DJII = rmm_input(j,ii) * factor * 0.5f;
             PJII = Pred(j,ii) * factor;
             VJII = Vred(j,ii) * factor;
             temp[0] = 2.0f * (dfac[0]*DJII + pfac[0]*PJII + tfac[0]*VJII);
             temp[1] = 2.0f * (dfac[1]*DJII + pfac[1]*PJII + tfac[1]*VJII);
             temp[2] = 2.0f * (dfac[2]*DJII + pfac[2]*PJII + tfac[2]*VJII);
             temp[3] = 2.0f * (dfac[3]*DJII + pfac[3]*PJII + tfac[3]*VJII);

             grdx += temp[0]*gfx[ii]*fv[j];
             grdx += temp[1]*( hxx[ii]*fv[j] + gfx[ii]*gfx[j] );
             grdx += temp[2]*( hxy[ii]*fv[j] + gfx[ii]*gfy[j] );
             grdx += temp[3]*( hxz[ii]*fv[j] + gfx[ii]*gfz[j] );

             grdy += temp[0]*gfy[ii]*fv[j];
             grdy += temp[1]*( hxy[ii]*fv[j] + gfy[ii]*gfx[j] );
             grdy += temp[2]*( hyy[ii]*fv[j] + gfy[ii]*gfy[j] );
             grdy += temp[3]*( hyz[ii]*fv[j] + gfy[ii]*gfz[j] );

             grdz += temp[0]*gfz[ii]*fv[j];
             grdz += temp[1]*( hxz[ii]*fv[j] + gfz[ii]*gfx[j] );
             grdz += temp[2]*( hyz[ii]*fv[j] + gfz[ii]*gfy[j] );
             grdz += temp[3]*( hzz[ii]*fv[j] + gfz[ii]*gfz[j] );
          }
       }
       ddx(nuc) += grdx;
       ddy(nuc) += grdy;
       ddz(nuc) += grdz;
    }
    for (int i = 0; i < local_atoms; i++) {
       sForce[i*3] -= ddx(i) * wp;
       sForce[i*3+1] -= ddy(i) * wp;
       sForce[i*3+2] -= ddz(i) * wp;
    }
   } // END POINTS

   // Acumulate forces for this group
   for (int i = 0; i < local_atoms; i++) {
      uint global_atom = this->local2global_nuc[i]; // da el indice del atomo LOCAL al GLOBAL
      F(global_atom,0) += sForce[i*3];
      F(global_atom,1) += sForce[i*3+1];
      F(global_atom,2) += sForce[i*3+2];
   }

   // Free Memory
   free(dens); dens = NULL;
   free(diff); diff = NULL;
   free(trad); trad = NULL;
   free(dfac); dfac = NULL;
   free(pfac); pfac = NULL;
   free(tfac); tfac = NULL;
   free(sForce); sForce = NULL;
   Pred.deallocate(); Vred.deallocate();
   rmm_input.deallocate();
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}


    









