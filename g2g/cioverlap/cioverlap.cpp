/*
    Copyright (C) 2008 Jiri Pittner <jiri.pittner@jh-inst.cas.cz> or <jiri@pittnerovi.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//this program computes overlaps of arbitrary wave functions expressed in Slater determinants or GUGA CSFs
//UHF version now available
//NO particular order of determinants is assumed - all dets are listed in Slaterfile and CI coefficients
//in matching order are in the input

#undef oldversion

#include <iostream>
#include <sys/types.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
//#include "gelfand.h"
#include "drt.h"
#include "readintegrals.h"

#include "options.h"

//LA
#include "vec.h"
#include "mat.h"
#include "smat.h"
#include "sparsemat.h"
#include "nonclass.h"
#include "qsort.h"

#include "cis_slatergen.h"

#undef DEBUG
#undef DEBUG2
#undef DEBUG3

#include <cstdlib>
#include <cmath>

template<class MAT>
struct twooverlaps
	{
	NRMat<typename LA_traits<MAT>::elementtype> alpha;
	NRMat<typename LA_traits<MAT>::elementtype> beta;
	int rowparity; //parity of permutation to get rows in alpha,beta order
	int colparity; //parity of permutation to get columns in alpha,beta order
	};

template<class MAT, class INDEX>
const twooverlaps<MAT>  slatersubmatrices(const MAT (&a)[2], const int nelec, const INDEX rows, const INDEX cols, const int ncore, const int nactive, const bool inversedorbitals)
{
twooverlaps<MAT> r;

const int ncoreshift=0; //even with frozen core, the numbering of MOs is NOT shifted

int nalpha=0;
int nbeta=0;
for(int j=0; j<nelec; ++j) if(rows[j]<0) ++nbeta; else ++nalpha;

r.alpha.resize(ncore+nalpha,ncore+nalpha);
r.beta.resize(ncore+nbeta,ncore+nbeta);
r.alpha.clear();
r.beta.clear();

//count the number of betas all alphas must transpose with when moved to the left in conservative order
r.rowparity=ncore*(ncore-1)/2; 
r.colparity=ncore*(ncore-1)/2;
int rowbeta=ncore;
int colbeta=ncore; 
for(int j=0; j<nelec; ++j) 
	{
	if(rows[j]<0) ++rowbeta; else r.rowparity += rowbeta;
	if(cols[j]<0) ++colbeta; else r.colparity += colbeta;
	}


//core-core block
for(int i=0; i<ncore; ++i)
        for(int j=0; j<ncore; ++j)
                {
                r.alpha(i,j) = a[0](i,j);
                r.beta(i,j) = a[1](i,j);
                }

//core-active and active-core
for(int i=0; i<ncore; ++i)
	{
	int jac=0;
	int jbc=0;
	int jar=0;
	int jbr=0;
        for(int j=0; j<nelec; ++j)
		{
		if(inversedorbitals)
                        {
                        if(cols[j]<0) r.beta(i,ncore+ jbc++) = a[1](i,nactive-abs(cols[j])+ncoreshift);
                        else          r.alpha(i,ncore+ jac++) = a[0](i,nactive-abs(cols[j])+ncoreshift);

                        if(rows[j]<0) r.beta(ncore+ jbr++,i) = a[1](nactive-abs(rows[j])+ncoreshift,i);
                        else          r.alpha(ncore+ jar++,i) = a[0](nactive-abs(rows[j])+ncoreshift,i);
                        }
                else
                        {
			if(cols[j]<0) r.beta(i,ncore+ jbc++) = a[1](i,abs(cols[j])+ncoreshift-1); 
			else	      r.alpha(i,ncore+ jac++) = a[0](i,abs(cols[j])+ncoreshift-1); 

			if(rows[j]<0) r.beta(ncore+ jbr++,i) = a[1](abs(rows[j])+ncoreshift-1,i);
			else	      r.alpha(ncore+ jar++,i) = a[0](abs(rows[j])+ncoreshift-1,i);
                        }
		
		}
	}


//active-active block
if(inversedorbitals)
{
int ia=0; int ib=0;
for(int i=0; i<nelec; ++i)
        {
        int ja=0; int jb=0;
        for(int j=0; j<nelec; ++j)
                {
		cout <<"here2\n";
                if(rows[i]<0 && cols[j]<0) r.beta(ncore+ib,ncore+jb) = a[1](nactive-abs(rows[i])+ncoreshift,nactive-abs(cols[j])+ncoreshift);
                if(rows[i]>0 && cols[j]>0) r.alpha(ncore+ia,ncore+ja) = a[0](nactive-abs(rows[i])+ncoreshift,nactive-abs(cols[j])+ncoreshift);
                if(cols[j]<0) ++jb; else ++ja;
                }
        if(rows[i]<0) ++ib; else ++ia;
        }
}
else
{
int ia=0; int ib=0;
for(int i=0; i<nelec; ++i)
	{
	int ja=0; int jb=0;
        for(int j=0; j<nelec; ++j)
		{
		
		if(rows[i]<0 && cols[j]<0) r.beta(ncore+ib,ncore+jb) = a[1](abs(rows[i])+ncoreshift-1,abs(cols[j])+ncoreshift-1);
		if(rows[i]>0 && cols[j]>0) r.alpha(ncore+ia,ncore+ja) = a[0](abs(rows[i])+ncoreshift-1,abs(cols[j])+ncoreshift-1);
		if(cols[j]<0) ++jb; else ++ja;
		}
	if(rows[i]<0) ++ib; else ++ia;
	}
}

return r;
}
	

#if defined(oldversion) || defined(DEBUG3)
template<class MAT, class INDEX>
const NRMat<typename LA_traits<MAT>::elementtype> slatersubmatrix(const MAT (&a)[2], const int nelec, const INDEX rows, const INDEX cols, const int ncore, const int nactive, const bool inversedorbitals)
{
const int ncoreshift=0; //even with frozen core, the numbering of MOs is NOT shifted

int ncore2=2*ncore;
NRMat<typename LA_traits<MAT>::elementtype> r(nelec+ncore2,nelec+ncore2);

//core-core block
for(int i=0; i<ncore2; ++i)
	for(int j=0; j<ncore2; ++j)
		{
		r(i,j) = (i^j)&1 ? 0 : a[i&1](i/2,j/2);
		}

//core-active and active-core
for(int i=0; i<ncore2; ++i)
	for(int j=0; j<nelec; ++j)
		{
		int si= i&1? -1:1;
		if(inversedorbitals)
			{
			r(i,ncore2+j) = si*cols[j]<0 ? 0. : a[i&1](i/2,nactive-abs(cols[j])+ncoreshift);
                        r(ncore2+j,i) = si*rows[j]<0 ? 0. : a[i&1](nactive-abs(rows[j])+ncoreshift,i/2);
			}
		else
			{
			r(i,ncore2+j) = si*cols[j]<0 ? 0. : a[i&1](i/2,abs(cols[j])+ncoreshift-1);
			r(ncore2+j,i) = si*rows[j]<0 ? 0. : a[i&1](abs(rows[j])+ncoreshift-1,i/2);
			}
		}

//active-active block
if(inversedorbitals)
for(int i=0; i<nelec; ++i)
        for(int j=0; j<nelec; ++j)
                r(ncore2+i,ncore2+j) = rows[i]*cols[j]<0?0.:a[rows[i]<0?1:0](nactive-abs(rows[i])+ncoreshift,nactive-abs(cols[j])+ncoreshift);
else
for(int i=0; i<nelec; ++i)
        for(int j=0; j<nelec; ++j)
                r(ncore2+i,ncore2+j) = rows[i]*cols[j]<0?0.:a[rows[i]<0?1:0](abs(rows[i])+ncoreshift-1,abs(cols[j])+ncoreshift-1);
return r;
}
#endif


struct part_eliminated
	{
	NRMat<REAL> submat;
	NRSMat<REAL> coeffs;
	REAL det;
	};

static part_eliminated *precompa=NULL, *precompb=NULL;

void prepare_elimination(part_eliminated **precomp, const  NRMat<REAL> &a, int inactive)
{
 *precomp  = new part_eliminated;
(*precomp)->submat = a.submatrix(0,inactive-1,0,inactive-1);
(*precomp)->coeffs.resize(inactive); 
(*precomp)->coeffs.clear();

cout <<"Constant part of overlap matrices\n"<<(*precomp)->submat;

//no pivoting necessary, assumed inactive orbitals do not change character and diagonal elements are about 0.99
for(int i=0; i<inactive-1; ++i)
	{
	REAL p=(*precomp)->submat(i,i);
	for(int j=i+1; j<inactive; ++j)
		{
		REAL q= -(*precomp)->submat(j,i)/p;
		(*precomp)->coeffs(i,j) = q;
		xaxpy(inactive,q,(*precomp)->submat[i],1,(*precomp)->submat[j],1);
		}
	}
REAL d=1.;
for(int i=0; i<inactive; ++i) d *= (*precomp)->submat(i,i);
(*precomp)->det=d;

cout <<"Pre-eliminated ovelap submatrix\n"<<(*precomp)->submat;
cout.flush();
cout <<"Elimination coefficients\n"<<(*precomp)->coeffs;
cout.flush();
}


REAL finish_elimination(const part_eliminated *precomp,  NRMat<REAL> &a, int inactive)
{
//store in the matrix the precomputed part
a.storesubmatrix(0,0,precomp->submat);

//perform the precomputed linear combinations on the right columns
for(int i=0; i<inactive-1; ++i)
        for(int j=i+1; j<inactive; ++j)
                xaxpy(a.ncols()-inactive,precomp->coeffs(i,j),a[i]+inactive,1,a[j]+inactive,1);


//finish elimination
//should a pivoting be done for accuracy (in principle an orbital swap could lead to division by a near-zero) ?
for(int i=0; i<a.nrows()-1; ++i)
	{
	REAL p=a(i,i);
	int jlow= i+1;
	if(inactive > jlow) jlow= inactive;
        for(int j=jlow; j<a.nrows(); ++j)
		{
		REAL q= -a(j,i)/p;
                xaxpy(a.ncols(),q,a[i],1,a[j],1);
		}
	}

#ifdef DEBUG
cout << "eliminated matrix\n"<<a;
#endif

//compute determinant
REAL d=precomp->det;
for(int i=inactive; i<a.nrows(); ++i) d *= a(i,i);
#ifdef DEBUG
cout << "determinant = "<<d<<endl;
#endif
return d;
}



void perform_overlap(lexindex i, lexindex j, unsigned long *skipcount,
	bool (&needpermute)[2],int nelec,const slaterbasis &slaters, const NRVec<int> (&bra2ketperm)[2],
	bool permuteslvectors, const NRVec<SPMatindex> &slperm,REAL screeningthr, 
	const NRMat<REAL> &brasl,const NRMat<REAL> &ketsl,const NRMat<REAL> (&Smo)[2],int ncore,int mo,
	bool inversedorbitals,unsigned long *computedcount,SparseMat<REAL> &Ssl, int inactive,
	NRMat<int> *screening_mask)
{
if(screening_mask)
	{
	if((*screening_mask).nrows() != brasl.ncols() || (*screening_mask).ncols() != ketsl.ncols()) laerror("error in screening_mask dimensions");
	}
			lexindex k;
			if(needpermute[0]||needpermute[1]) //only general need in principle
				{
				bool waspermuted=false;
				slaterdet permuted(nelec);
				for(int e=0; e<nelec; ++e) 
					{
					spinorbindex ee=slaters(j,e);
					permuted[e]=bra2ketperm[ee<0?1:0][abs(ee)];
					if(ee<0) permuted[e]= -permuted[e];
					if(permuted[e]!=ee) waspermuted=true;
					}
				
				if(waspermuted)
					{
					permuted.sort();
					k=permuted.find(slaters,slaterorder); 
#ifdef DEBUG
					cout << "original determinant "<<slaters.row(j)<<endl;
					cout << "permuted determinant "<<j<<" "<<k<<" : "<<permuted<<endl;
					cout <<" ci size = "<<slaters.nrows()<<endl;
#endif
					if(k == (lexindex)-1) 
						{
#ifdef DEBUG
						cout <<"WARNING: permuted slater determinant outside CI, computed overlap will be inaccurate\n";
#endif
						++ (*skipcount);
						return;
						}
					}
				else
					k=j;
				}
			else
				k=j;

			lexindex itrue = permuteslvectors ? slperm[i] : i ;
			lexindex ktrue = permuteslvectors ? slperm[k] : k ;

			//implement screening
			if(screeningthr>0)
				{
				REAL s=0.;
				for(int a=0; a<brasl.ncols(); ++a) 
					for(int b=0; b<ketsl.ncols(); ++b) 
					    {
					    if(!screening_mask || (*screening_mask)(a,b) )
						{
						s += abs(brasl[itrue][a] * ketsl[ktrue][b]);
						}
					    }
				if(s<screeningthr)
					{
					++ (*skipcount);
					return;
					}
				}


{

#ifdef oldversion
NRMat<REAL> overlaps=slatersubmatrix(Smo,nelec,slaters[i],slaters[k],ncore,mo,inversedorbitals);
REAL dd=determinant_destroy(overlaps);
#else
#ifdef DEBUG3
NRMat<REAL> overlaps=slatersubmatrix(Smo,nelec,slaters[i],slaters[k],ncore,mo,inversedorbitals);
#ifdef DEBUG2
if(i==0 && k==13)
{
cout <<"Ssl "<<i<<" "<<k<<" (j="<<j<<") : " <<itrue<<" "<<ktrue<<" = "<<endl;
cout <<"Overlap submatrix original algorithm " <<overlaps<<endl;
}
#endif
REAL ddorig=determinant_destroy(overlaps);
cout <<"original_result = "<<ddorig<<endl;
#endif


twooverlaps<NRMat<REAL> > overlaps2 = slatersubmatrices(Smo,nelec,slaters[i],slaters[k],ncore,mo,inversedorbitals);

#ifdef DEBUG2
if(i==0 && k==13)
{
cout <<"Ssl "<<i<<" "<<k<<" (j="<<j<<") : " <<itrue<<" "<<ktrue<<" = "<<endl;
cout <<"i= ";for(int ii=0;ii<nelec;++ii) cout <<slaters[i][ii]<<" "; cout <<endl;
cout <<"k= ";for(int ii=0;ii<nelec;++ii) cout <<slaters[k][ii]<<" "; cout <<endl;
cout <<"Overlap submatrices alpha, beta:\n"<<overlaps2.alpha<<endl<<overlaps2.beta<<endl;
}
#endif
REAL da,db;
if(inactive>1) //simplified computation from a precomputed intermediate for efficiency
	{
	if(!precompa) prepare_elimination(&precompa,overlaps2.alpha,inactive);
	if(!precompb) prepare_elimination(&precompb,overlaps2.beta,inactive);
	da = finish_elimination(precompa,overlaps2.alpha,inactive);
	db = finish_elimination(precompb,overlaps2.beta,inactive);
	}
else //full determinant computation
	{
	da=determinant_destroy(overlaps2.alpha);
	db=determinant_destroy(overlaps2.beta);
	}

REAL dd=da*db; if((overlaps2.rowparity+overlaps2.colparity)&1) dd= -dd;
#ifdef DEBUG
cout <<"da db rowparity colparity = "<<da<<" "<<db<<" "<<overlaps2.rowparity<<" "<<overlaps2.colparity<<endl;
cout <<"full_result = "<<dd<<endl;
#ifdef DEBUG3
if(fabs(ddorig-dd) > 1e-13) cout <<"INTERNAL ERROR\n";
#endif
#endif
#endif //oldversion
++ *computedcount;
Ssl.add(itrue,ktrue,dd);
}
}


//suitable for small cis, we do not care much about memory needs, just try not to really waste
//all paths go to this routine

NRMat<REAL> CIoverlap_slater(
	bool uhf,
	const NRMat<REAL> &brasl, //slaterCIsize x n
	const NRMat<REAL> &ketsl, //slaterCIsize x m
	const NRMat<REAL> &Sraw, //2*AO x 2*AO (generated by artificial doubled geometry input)
	const NRMat<REAL> (&braLCAO)[2], //AO x MO
        const NRMat<REAL> (&ketLCAO)[2], //AO x MO
	int slaterfile, //file descriptor for slater basis 
	int excitlistfile, //file descriptor for list of not neglected excitations
        const int nelec,
        const int mo, //active MOs
        const int ncore=0,
        const int ndiscarded=0,
	bool make_excitlist=true,
	const int truncation=4, //neglect overlap of Slater dets differing in more than ... spinorbitals (take into account arbitrary MO swapping)
	const bool inversedorbitals=true, //lowest MOs are at the top of DRT when true
	const int slaterorder=1,
        const bool reorderslbasis=true, //the slater basis was not sorted before
	const int slaterpermfile=-1, //file handle for permutation of slater overlap matrix (slater CI vectors will remain same, slater basis is permuted)
	REAL screeningthr=0,
	NRMat<int> *activerange=NULL,
	int inactive=0, //for performance improving precompute a constant part of the determinants 
	NRMat<int> *screening_mask = NULL
	) 
{

int motot=mo+ncore+ndiscarded;
cout <<"Number of frozen+active+discarded MOs = "<<motot<<endl;
int ao=braLCAO[0].nrows();
if(ao!= braLCAO[1].nrows() || ao!= ketLCAO[0].nrows() ||  ao!= ketLCAO[1].nrows()) laerror("inconsistent row dimension of LCAO matrices");
if(braLCAO[0].ncols() != braLCAO[1].ncols() || braLCAO[0].ncols() != ketLCAO[0].ncols() ||  braLCAO[0].ncols() != ketLCAO[1].ncols() ) laerror("inconsistent column dimension of LCAO matrices");

lexindex sl=(lexindex)brasl.nrows();


if(activerange)
	{
	if((*activerange)(0,0)<=0||(*activerange)(0,1)<=0||(*activerange)(1,0)<=0||(*activerange)(1,1)<=0||
		(*activerange)(0,0)>mo||(*activerange)(0,1)>mo||(*activerange)(1,0)>mo||(*activerange)(1,1)>mo) laerror("illegal active orbital range specified");
	}

//consistency checks
if(ao!=ketLCAO[0].nrows()
|| ao!=ketLCAO[1].nrows()
|| ao!=braLCAO[0].nrows()
|| ao!=braLCAO[1].nrows()
|| 2*ao!=Sraw.nrows()
|| motot!=braLCAO[0].ncols()
|| motot!=ketLCAO[0].ncols()
|| motot!=braLCAO[1].ncols()
|| motot!=ketLCAO[1].ncols()
|| sl!=(lexindex)ketsl.nrows()
  ) 
	{
	cout <<"ao " <<ao <<" ketLCAO[0].nrows() "<<ketLCAO[0].nrows()<< endl;
	cout <<"ao " <<ao <<" ketLCAO[1].nrows() "<<ketLCAO[1].nrows()<< endl;
	cout <<"ao " <<ao <<" braLCAO[0].nrows() "<<braLCAO[0].nrows()<< endl;
	cout <<"ao " <<ao <<" braLCAO[1].nrows() "<<braLCAO[1].nrows()<< endl;
	cout <<"2*ao "<<2*ao<<" Sraw.nrows() "<<Sraw.nrows()<<endl;
	cout <<"motot "<<motot<<" braLCAO[0].ncols() "<<braLCAO[0].ncols()<<" ketLCAO[0].ncols() "<<ketLCAO[0].ncols()<<endl;
	cout <<"motot "<<motot<<" braLCAO[1].ncols() "<<braLCAO[1].ncols()<<" ketLCAO[1].ncols() "<<ketLCAO[1].ncols()<<endl;
	cout <<"sl "<<sl<<" ketsl.nrows() "<<ketsl.nrows()<<endl;
	laerror("inconsistent dimensions in CIoverlap");
	}

//calculate one-particle overlaps in the MO basis
NRMat<REAL> Smo[2];
Smo[0].resize(motot,motot);
Smo[1].resize(motot,motot);
{
NRMat<REAL> Sao12=Sraw.submatrix(0,ao-1,ao,2*ao-1);
cout <<"test Sao12\n"<<Sao12;
NRMat<REAL> tmp;
for(int spin=0; spin<2; ++spin)
	{
	tmp=Sao12*ketLCAO[spin];
	Smo[spin].gemm(0.,braLCAO[spin],'t',tmp,'n',1.); 
	}
}

/* for test symmetrize Smo and test symmetry of Ssl
Smo[0]= (Smo[0]+Smo[0].transpose())*.5;
Smo[1]= (Smo[1]+Smo[1].transpose())*.5;
*/

for(int spin=0; spin<2; ++spin) cout <<"test Smo["<<spin<<"]\n"<<Smo[spin];

//#define UNIT_SMO
#ifdef UNIT_SMO
//try to put here a fake diagonal Smo 
for(int spin=0; spin<2; ++spin)
			{
			for(int kk=0; kk<Smo[spin].nrows(); ++kk)
				for(int ll=0; ll<Smo[spin].ncols(); ++ll)
					{							
					if(kk!=ll) Smo[spin](kk,ll)=0; else
						{
						if(Smo[spin](kk,ll)<0) Smo[spin](kk,ll)= -1; else Smo[spin](kk,ll)=1;
						}
					}
			}
			cout <<"Adjusted Smo["<<spin<<"]\n"<<Smo[spin];
#endif


//find out largest overlap mapping between bra and ket MOs for later need
//this should preferably be based not on Smo but on overlap of Loewdin-normalized LCAO at the two geometries
//thus it will yield identity for rotated-translated molecule
bool needpermute[2];
needpermute[0]=false;
needpermute[1]=false;
NRVec<int> bra2ketperm[2];
NRVec<int> bra2ket[2];
bra2ket[0].resize(motot);
bra2ket[1].resize(motot);
{
NRMat<REAL> Sao11=realsqrt(Sraw.submatrix(0,ao-1,0,ao-1));
NRMat<REAL> Sao22=realsqrt(Sraw.submatrix(ao,2*ao-1,ao,2*ao-1));
for(int spin=0; spin<2; ++spin)
{
NRMat<REAL> braL = Sao11 * braLCAO[spin];
NRMat<REAL> ketL = Sao22 * ketLCAO[spin];
NRMat<REAL> Sloewdin(motot,motot);
Sloewdin.gemm(0.,braL,'t',ketL,'n',1.);
cout <<"test Sloewdin["<<spin<<"]\n"<<Sloewdin;

NRVec<int> counts(motot);
NRVec<int> mophase(motot);
counts=0;
mophase=0;
for(int i=0; i<motot; ++i)
	{
	int jmaxforbid = -1;
	findmax:
	int jmax= -1;
	int pmax= 1;
	REAL smax= -1.;
	for(int j=0; j<motot; ++j)
		if(abs(Sloewdin(i,j)) > smax && j != jmaxforbid)
			{
			jmax=j;
			smax=abs(Sloewdin(i,j));
			pmax = Sloewdin(i,j)>0. ?1 : -1;
			}
	bra2ket[spin][i]=jmax;
	mophase[i] = pmax;
	if(counts[jmax] && jmaxforbid== -1) //this already had a big overlap and we were not in the second search already
		{
		//find second biggest one
		jmaxforbid=jmax;
		goto findmax;
		}
	else counts[jmax]++;
	}

cout <<"test bra2ket["<<spin<<"], counts, mophase\n"<<bra2ket[spin]<<counts<<mophase;
for(int i=0; i<motot; ++i) if(counts[i]!=1) 
	{
	cout << "cannot determine MO permutation, orbital rotations do not represent simple orbital swaps\n";
	for(int j=0; j<motot; ++j) bra2ket[spin][j]=j;
	goto skippermute;
	}

cout <<"Smo["<<spin<<"] (permuted) diagonal elements\n";
for(int i=0; i<motot; ++i) cout << i<<" : "<<Smo[spin](i,bra2ket[spin][i])<<endl;

for(int i=ncore; i<motot-ndiscarded; ++i) if(bra2ket[spin][i]!=i) 
	{
	needpermute[spin]=true;
	if(bra2ket[spin][i]<ncore||bra2ket[spin][i]>=motot-ndiscarded) {cout <<"Warning: MO swap outside of active space occured for spin "<<spin<<"\n";}
	}

skippermute:

//adjust the bra2ket permutation to the CI-active orbitals
bra2ketperm[spin].resize(mo+1); //index from 1 over active
	{
	bra2ketperm[spin][0]= -1; //dummy
	if(inversedorbitals) for(int i=0; i<mo; ++i) 
			{
			bra2ketperm[spin][i+1]=mo-(bra2ket[spin][ncore+mo-i-1]-ncore); 
			}
	else for(int i=0; i<mo; ++i) bra2ketperm[spin][i+1]=bra2ket[spin][ncore+i]-ncore+1;
	}
cout <<"bra2ketperm["<<spin<<"] "<<bra2ketperm[spin]<<endl;
cout <<"needpermute["<<spin<<"] "<<needpermute[spin]<<endl;

}//spin
}//local scope

cout <<"application of the permutation has been switched off, since typically the determinants from permuted orbitals are outside of the CI space anyway\n";
needpermute[0]=false;
needpermute[1]=false;

//read the slater basis into memory
cout << "Slater basis dimensions are "<<sl<<" "<<nelec<<endl;
cout <<"(total number of electrons including the ones in frozen orbitals is "<<nelec+2*ncore<<")"<<endl;
slaterbasis slaters(sl,nelec);
bool permuteslvectors=false;
{
struct stat buf;
if(! fstat(slaterpermfile,&buf)) if(buf.st_size>0) permuteslvectors=true;
}
if(reorderslbasis && !permuteslvectors) //was not reordered previously
	{
	//read the input basis
	lseek(slaterfile,0,SEEK_SET);
	slaters.get(slaterfile,false);
	//check that the whole slaterfile has been exhausted by the read to prevent mismatch in the number of electrons
	{
	unsigned char dummy;
	int r = read(slaterfile,&dummy,1);
	if (r>0) laerror("slaterfile is larger than expected, check the number of electrons in cioverlap input");
	if (r<0) laerror("unexpected error when testing slaterfile eof");
	}
	//check for illegal content (orbital number 0)	
	if(slaters.checkzero()) laerror("malformed slaterfile encountered, perhaps cipc/mcpc was compiled without --assume byterecl");
	//sort
	NRVec<SPMatindex> slperm(sl);
	for(lexindex i=0; i<sl; ++i) slperm[i]=i;
	slsort_sldetbase=&slaters;
        slsort_permbase=&slperm[0];
        slsort_nelectrons=nelec;
        genqsort(0,(int)sl-1,slaterorder>0 ? slsort_slbascmp : slsort_slbascmp2 ,slsort_slbasswap);
	//save the new basis
        lseek(slaterfile,0,SEEK_SET);
        for(lexindex i=0; i<sl; ++i)
                if(slsort_nelectrons*(int)sizeof(spinorbindex) != write(slaterfile, slaters[slperm[i]] ,slsort_nelectrons*sizeof(spinorbindex))) laerror("write error in CIoverlap_slater");

	//save the permutation, to be used in contraction with CI coefficients
	lseek(slaterpermfile,0,SEEK_SET);
        slperm.put(slaterpermfile,true);
	permuteslvectors=true;
	}
lseek(slaterfile,0,SEEK_SET);
slaters.get(slaterfile,false);
{
unsigned char dummy;
int r = read(slaterfile,&dummy,1);
if (r>0) laerror("slaterfile is larger than expected, check the number of electrons in cioverlap input");
if (r<0) laerror("unexpected error when testing slaterfile eof");
}

//check for illegal content (orbital number 0)
if(slaters.checkzero()) laerror("malformed slaterfile encountered, perhaps cipc was compiled without --assume byterecl");
#ifdef DEBUG
cout <<"slater basis "<<slaters<<endl;
#endif

//precompute and store in a file a list giving all determinants to consider from given one in the Ssl computation
//Note: symmetry of non-neglected Ssl entries is not in general present, if the bra2ket permutation has cycles longer than two

if(truncation<0) 
	{
	cout <<"Omitting excitlistfile generation\n";
	make_excitlist=false;
	}

if(make_excitlist)
	{
	clock_t t0=clock();
	cout <<"Generating excitation list for CI size " <<sl<<endl;
	lseek(excitlistfile,0,SEEK_SET);
	lexindex maxlen=0;
	if(sizeof(lexindex)!=write(excitlistfile,&maxlen,sizeof(lexindex))) laerror("write error");
	NRVec<lexindex> allowedlist(sl);
	NRVec<unsigned char> n_act_elec((unsigned char)0,sl);
	if(activerange) //count electrons in active orbitals for each determinant
	    for(lexindex i=0; i<sl; ++i)
		{
		for(int e=0; e<nelec; ++e) 
			{
			int imo=slaters[i][e];
			if(imo>0 && imo>=(*activerange)(0,0) && imo<=(*activerange)(0,1)
				|| imo<0 && imo>=(*activerange)(1,0) && imo<=(*activerange)(1,1)) ++n_act_elec[i];
			}
		}
	for(lexindex i=0; i<sl; ++i) //prepare the list for each slater det.
		{
		NRVec<char> notocc_i(2*(mo+ncore)+1);  //is in principle a bitvector but the "unpacked form" should be more efficient
		for(int o=0; o<2*(mo+ncore)+1; ++o) notocc_i[o]=1;
		for(int e=0; e<nelec; ++e) notocc_i[mo+ncore+slaters[i][e]]=0;
		//setup an active space, where a difference in occupation does not count
		//not to degrade performance it is also handeled via notocc_i flags
		if(activerange) 
			{
			for(int o=(*activerange)(0,0); o<=(*activerange)(0,1); ++o) notocc_i[mo+ncore+o]=0;
			for(int o=(*activerange)(1,0); o<=(*activerange)(1,1); ++o) notocc_i[mo+ncore+o]=0;
			}
		//
		allowedlist[0]=i;
		lexindex ntot=1;
        	for(lexindex j=i+1; j<sl; ++j) 
			{
			//if(slater_excitation(nelec,slaters[i],slaters[j],false,mo)>truncation) goto skipthis ;
			//replaced by more optimized code
			int include_it = truncation;
			int activediff = (int)n_act_elec[j] - (int)n_act_elec[i];
			spinorbindex *sj = slaters[j];
			if(abs(activediff) > truncation) goto skipthis;
			if(activediff>0) include_it -= activediff; //account for electrons excited into active orbitals
			if(include_it<0) goto skipthis;
			for(int e=0; e<nelec; ++e)
				{
				if(notocc_i[mo+ncore+sj[e]]) --include_it;
				if(include_it<0) goto skipthis;
				}	
			allowedlist[ntot++]=j;
			skipthis:;
			}
		cout <<"For "<<i<< " "<<ntot <<endl;
		//write the list for i-th slater
		if(ntot>maxlen) maxlen=ntot;
		if(sizeof(lexindex)!=write(excitlistfile,&ntot,sizeof(lexindex))) laerror("write error");
		if((lexindex)ntot*sizeof(lexindex)!=(lexindex)write(excitlistfile,&allowedlist[0],ntot*sizeof(lexindex))) laerror("write error");
		}
	lseek(excitlistfile,0,SEEK_SET);
	if(sizeof(lexindex)!=write(excitlistfile,&maxlen,sizeof(lexindex))) laerror("write error");
	cout <<"make_excitlist cpu time "<<(1.*clock()-t0)/CLOCKS_PER_SEC<<endl;
	}


//prepare many-particle overlap matrix in the Slater basis
NRVec<SPMatindex> slperm;
if(permuteslvectors)
	{
	lseek(slaterpermfile,0,SEEK_SET);
	slperm.get(slaterpermfile,true);
	if(sl!= (lexindex)slperm.size()) laerror("inconsistent slaterpermfile");
	}

SparseMat<REAL> Ssl(sl,sl);
unsigned long skipcount=0;
unsigned long computedcount=0;
clock_t t0=clock();
if(truncation>=0) lseek(excitlistfile,0,SEEK_SET);
lexindex maxlen=0;
if(truncation>=0)
	{
	if(sizeof(lexindex)!=read(excitlistfile,&maxlen,sizeof(lexindex))) laerror("excitlist read error 1");
	if(maxlen==0) laerror("unfinished excitlistfile encountered");
	}
else
	maxlen=1;
NRVec<lexindex> excitlist(maxlen);
cout<<"Number of determinants in the wavefunctions = "<<sl<<endl;
cout<<"Done 0%\n"; cout.flush();
int percent=0;
for(lexindex i=0; i<sl; ++i)
	{
#ifdef DEBUG
cout<< "Processing determinant "<<i<<endl;
#endif
	if((i*100)/sl > percent)
		{
		percent = (i*100)/sl;
		cout<<" "<<percent<<"%"; cout.flush();
		}
	lexindex ntot;
	if(truncation>=0)
		{
		if(sizeof(lexindex)!=read(excitlistfile,&ntot,sizeof(lexindex))) laerror("excitlist read error 2");
		if(ntot*sizeof(lexindex)!=(lexindex)read(excitlistfile,&excitlist[0],ntot*sizeof(lexindex))) laerror("excitlist read error 3");
		}
	for(lexindex jj= (truncation>=0? 0 : i) ; jj<(truncation>=0? ntot :sl); ++jj)	
			{
			lexindex jjx;
			jjx = truncation>=0 ? excitlist[jj] : jj;
#ifdef DEBUG
cout<< "Processing its peer "<<jj<<" "<<jjx<<endl;
#endif
			perform_overlap(i,jjx,&skipcount,needpermute,nelec,slaters,bra2ketperm,permuteslvectors,slperm,screeningthr,brasl,ketsl,Smo,ncore,mo,inversedorbitals,&computedcount,Ssl,inactive,screening_mask);
			if(i!=jjx) perform_overlap(jjx,i,&skipcount,needpermute,nelec,slaters,bra2ketperm,permuteslvectors,slperm,screeningthr,brasl,ketsl,Smo,ncore,mo,inversedorbitals,&computedcount,Ssl,inactive,screening_mask);
			}
	}
cout<< " 100%\noverlap calculation CPU time "<<(1.*clock()-t0)/CLOCKS_PER_SEC<<endl;
cout.flush();
slaters.resize(0,0);// free memory

#ifdef DEBUG
NRMat<REAL>Ssltest(Ssl);
/*
for(lexindex i=0; i<sl; ++i) 
    {
	cout <<"for "<<i<<endl;
	for(lexindex j=0; j<sl; ++j)
	{
	if(abs(Ssltest(i,j))>.01) cout <<"Ssl "<<i<<" "<<j<<" = "<<Ssltest(i,j)<<endl;
	//cout <<"Ssl asymetry "<<i<<" "<<j<<" "<<Ssltest(i,j)-Ssltest(j,i) <<" from "<<Ssltest(i,j)<<endl;
	}
	cout <<endl;
    }
*/
cout <<"Ssl diagonal elements\n";
for(lexindex i=0; i<sl; ++i) cout <<"det "<<i<<": "<<Ssltest(i,i)<<endl;
#endif

cout <<"Number of computed overlap determinant pairs = "<<computedcount<<endl;
if(screeningthr>0) cout <<"Number of skipped determinant pairs at screening threshold "<<screeningthr<<" = "<<skipcount<<endl;


//contract bras and kets with the overlaps in Slater basis 
NRMat<REAL> tmp=Ssl*ketsl;
Ssl.resize(0,0);//free memory
NRMat<REAL> result(brasl.ncols(),ketsl.ncols());
result.gemm(0.,brasl,'t',tmp,'n',1.);

//zero out inaccurate elements according to screening_mask
if(screening_mask)
	{
	int i,j;
	for(i=0; i<result.nrows(); ++i) for(j=0; j<result.ncols(); ++j)
		{
		if((*screening_mask)(i,j) == 0) result(i,j)=0;
		}
	}

return result;
}


NRMat<int> *read_screening_mask(char *screening_mask_file,int nbras, int nkets)
{
if(!screening_mask_file) return NULL;
NRMat<int> *r = new NRMat<int>(0,nbras,nkets);
FILE *f = fopen(screening_mask_file,"r");
if(!f) laerror("cannot read screening_mask_file");
char line[1024];
fgets(line, 1024,f);
int x1,x2,x3,x4;
while(4==fscanf(f,"%d %d %d %d",&x1,&x2,&x3,&x4))
	{
	(*r)(x2-1,x4-1)=1;
	(*r)(x4-1,x2-1)=1;
	}
fclose(f);
for(int i=0; i<min(nbras,nkets); ++i) (*r)[i][i]=1;
cout <<"Screening mask\n"<< *r<<endl;
return r;
}


NRMat<REAL> CIoverlap(
	bool uhf,
	const NRMat<REAL> &bras, //CIsize x n
	const NRMat<REAL> &kets, //CIsize x m
	const NRMat<REAL> &Sraw, //2*AO x 2*AO (generated by artificial doubled geometry input)
	const NRMat<REAL> (&braLCAO)[2], //AO x MO
	const NRMat<REAL> (&ketLCAO)[2], //AO x MO
	const drt &d,
	int slaterfile, //file descriptor for slater basis 
	int citrfile, //file descriptor for CItrafo matrix
	int excitlistfile, //file descriptor for list of not neglected excitations
        const int ncore=0,
        const int ndiscarded=0,
	bool generate_trafo=true,
	bool make_excitlist=true,
	const int truncation=4, //neglect overlap of Slater dets differing in more than ... spinorbitals (take into account arbitrary MO swapping)
	const bool inversedorbitals=true, //lowest MOs are at the top of DRT when true
	const bool inverselex=false,
	const int slaterorder=1,
	REAL screeningthr=0,
	NRMat<int> *activerange=NULL,
	int inactive=0,
	NRMat<int> *screening_mask = NULL
	) 
{
//generate and save, or read from files Slater basis and transformation matrix
//this computation can be done only once and remains constant
lexindex ci=d.cisize();
NRMat<REAL> brasl, ketsl;
SparseMat<REAL> citrafo;
if(generate_trafo)
	{
	citrafo=CItrafo(slaterfile,d,slaterorder,inverselex);
	citrafo.put(citrfile,true);
	}
else
	{
	lseek(citrfile,0,SEEK_SET);
	citrafo.get(citrfile,true);
	if((lexindex)citrafo.ncols()!=ci) laerror("inconsistent dimension of citrafo");
	}

//transform CI vectors to the Slater basis
brasl=citrafo*bras;
ketsl=citrafo*kets;
citrafo.resize(0,0); //free memory

return CIoverlap_slater(uhf,brasl,ketsl,Sraw,braLCAO,ketLCAO,slaterfile,excitlistfile,d.electrons(),d.orbnum(),ncore,ndiscarded,make_excitlist,truncation,inversedorbitals,slaterorder,false,-1,screeningthr,activerange,inactive,screening_mask);
}

NRMat<REAL> docislateroverlap(bool uhf, REAL screeningthr, int excitrank,NRMat<int> *activerange=NULL,int inactive=0, char *screening_mask_file=NULL)
{

//NOTE: ncore are only such core orbitals which are not included in slaterfile
//NOTE: slaterfile numbering is expected to start from ncore+1
//NOTE  inactive orbitals by -i option are NOT counted as core here, they are only approximation to make computation faster
//NOTE: discarded orbitals have to be included so that all dimensions fit
//NOTE: nelec should NOT include the 2*ncore frozen electrons
//
int nbas,ncore,ndisc,nactive,nelec;
//cin >>nbas >>ncore >>ndisc >>nelec; //NOTE: nelec here does not have to distinguish alpha and beta
nbas=92; ncore=0; ndisc=0; nelec=42;
nactive=nbas-ncore-ndisc;
cout << "Number of frozen core orbitals (not in slaterfile) = "<<ncore<<endl;
cout << "Number of discarded virtual orbitals = "<<ndisc<<endl;
cout << "Number of not frozen electrons = " <<nelec<<endl;

if(inactive*2>=nelec) laerror("more inactive doubly occ orbitals than nelec/2");

bool include_occ=1, number_occ=1;
int nocc[2], nvirt[2], inactive_occ=0, inactive_virt=0;
nocc[0]=nocc[1]=21; nvirt[0]=nvirt[1]=71;
do_cis_slater(uhf,ncore,nocc,nvirt,inactive_occ,inactive_virt,include_occ,number_occ);

int slaterfile=open("slaterfile",O_RDWR|O_LARGEFILE,0777); if(slaterfile<0) {perror("cannot open slaterfile");laerror("IO error");}
int excitlistfile=open("excitlistfile",O_CREAT|O_LARGEFILE|O_RDWR,0777); if(excitlistfile<0) {perror("cannot open excitlistfile");laerror("IO error");}
int slaterpermfile=open("slaterpermfile",O_CREAT|O_RDWR|O_LARGEFILE,0777); if(slaterpermfile<0) {perror("cannot open slaterpermfile");laerror("IO error");}

NRMat<REAL> bras,kets;
int f;
f=open("eivectors1",O_RDONLY|O_LARGEFILE); if(f<0) {perror("cannot open eivectors1");laerror("IO error");}
bras.get(f,true,true);
close(f);
f=open("eivectors2",O_RDONLY|O_LARGEFILE); if(f<0) {perror("cannot open eivectors2");laerror("IO error");}
kets.get(f,true,true);
close(f);

#ifdef DEBUG
cout <<"test eivectors1 (bras) "<<bras<<endl;
cout <<"test eivectors2 (kets) "<<kets<<endl;
#else
cout <<"Dominant contributions to bras\n";
for(int k=0; k<bras.ncols(); ++k)
	{
	cout <<"bra state "<<k<<":\n";
	for(lexindex l=0; l<bras.nrows(); ++l)
		if(abs(bras(l,k))>0.1) cout <<"det "<<l<<" "<<bras(l,k)<<endl;
	}
cout <<"Dominant contributions to kets\n";
for(int k=0; k<kets.ncols(); ++k)
        {
        cout <<"ket state "<<k<<":\n";
        for(lexindex l=0; l<kets.nrows(); ++l)
                if(abs(kets(l,k))>0.1) cout <<"det "<<l<<" "<<kets(l,k)<<endl;
        }
#endif

NRMat<int> *screening_mask = read_screening_mask(screening_mask_file,bras.ncols(), kets.ncols());

cout <<"Coefficient-overlap bra-ket:\n";
cout << bras.transpose() * kets;

NRMat<REAL> Sraw;
Sraw.resize(nbas*2,nbas*2);
ifstream rf; rf.open("Srawfile.txt");
REAL value=0.;
for (int ii=0; ii<nbas*2; ii++)
   for (int jj=0; jj<nbas*2; jj++) {
        rf >> value; 
        Sraw(ii,jj) = value;
   }
rf.close();

NRMat<REAL> braLCAO[2],ketLCAO[2];
braLCAO[0].resize(nbas,nbas); braLCAO[1].resize(nbas,nbas);
ketLCAO[0].resize(nbas,nbas); ketLCAO[1].resize(nbas,nbas);
ifstream bf; bf.open("braLCAOfile.txt");
ifstream kf; kf.open("ketLCAOfile.txt");
REAL v1=0.0, v2=0.0;
for (int ii=0; ii<nbas; ii++)
   for (int jj=0; jj<nbas; jj++) {
        bf >> v1;
        kf >> v2;
        braLCAO[0](ii,jj)=braLCAO[1](ii,jj)=v1;
        ketLCAO[0](ii,jj)=ketLCAO[1](ii,jj)=v2;
   }
bf.close();
kf.close();

/*
cin >> braLCAO[0];
if(uhf) cin >> braLCAO[1]; else  braLCAO[1] = braLCAO[0];
cin >> ketLCAO[0];
if(uhf) cin >> ketLCAO[1]; else ketLCAO[1] = ketLCAO[0];
*/


bool make_excitlist=true;
{
struct stat64 buf;
if(fstat64(excitlistfile,&buf)) laerror("cannot stat excitlistfile");
if(buf.st_size>0) make_excitlist=false;
}


NRMat<REAL> SCI = CIoverlap_slater(uhf, bras, kets, Sraw, braLCAO, ketLCAO, slaterfile, excitlistfile, nelec, nactive, ncore, ndisc, make_excitlist, excitrank,inverseorb,slaterorder,false /*we do not need lexical order now*/,slaterpermfile,screeningthr,activerange,inactive,screening_mask);

close(slaterfile);
close(slaterpermfile);
close(excitlistfile);

return SCI;
}

void cioverlap( )
{
clock_t clock0 = clock();

cout.setf(ios::fixed);
//cout.setf(ios::scientific);
    cout.precision(12);

bool uhf=0;
/*
if(argc>1  && !strcmp(argv[1],"-U"))
        {
	uhf=1;
        argc -=1; argv+=1;
        }
*/


char *screening_mask_file=NULL;
screening_mask_file="transmomin";
/*
if(argc>2  && !strcmp(argv[1],"-s"))
        {
	screening_mask_file=strdup(argv[2]);
	argc -=2; argv+=2;
        }
*/



bool alignrows=1;
bool useoldphase=1;
/*
if(argc>1  && !strcmp(argv[1],"-a"))
        {
	alignrows=1;
	useoldphase=1;
	--argc; ++argv;
	}
*/

bool aligncolumns=0;
/*
if(argc>1  && !strcmp(argv[1],"-b"))
        {
        aligncolumns=1;
	useoldphase=1;
        --argc; ++argv;
        }

if(argc>1  && !strcmp(argv[1],"-o"))
        {
        useoldphase=0;
        --argc; ++argv;
        }
*/

REAL screeningthr=5e-4; 
/*
if(argc>2 && !strcmp(argv[1],"-t"))
	{
	sscanf(argv[2],"%lf",&screeningthr);
	argc -=2; argv+=2;
	}
*/

int excitrank=-1;
/*

if(argc>2 && !strcmp(argv[1],"-e"))
        {
        sscanf(argv[2],"%d",&excitrank);
        argc -=2; argv+=2;
        }
*/

int inactive=0;
/*
if(argc>2 && !strcmp(argv[1],"-i"))
        {
        sscanf(argv[2],"%d",&inactive);
        argc -=2; argv+=2;
        }
*/

NRMat<int> *activerange=NULL;

/*
if(argc>(uhf?5:3) && !strcmp(argv[1],"-A"))
        {
	activerange = new NRMat<int>(2,2);
	if(uhf)
		{
		sscanf(argv[2],"%d",&(*activerange)(0,0));
                sscanf(argv[3],"%d",&(*activerange)(0,1));
		sscanf(argv[2],"%d",&(*activerange)(1,0));
                sscanf(argv[3],"%d",&(*activerange)(1,1));
		argc -=5; argv+=5;
		}
	else
		{
        	sscanf(argv[2],"%d",&(*activerange)(0,0));
        	sscanf(argv[3],"%d",&(*activerange)(0,1));
		(*activerange)(1,0)= -(*activerange)(0,1);
		(*activerange)(1,1)= -(*activerange)(0,0);
        	argc -=3; argv+=3;
		}
        }
*/

cout << "Number of (not frozen) inactive orbitals = "<<inactive<<endl;
cout << "Excitation rank threshold = "<<excitrank<<endl;
cout << "Screening threshold = "<<screeningthr<<endl;

//int type=0; if(argc>1) sscanf(argv[1],"%d",&type);
NRMat<REAL> SCI;

//if(type) SCI = docigugaoverlap(uhf,screeningthr,excitrank,activerange,inactive,screening_mask_file);
SCI = docislateroverlap(uhf,screeningthr,excitrank,activerange,inactive,screening_mask_file);

cout << "CI overlap before alignment\n";
cout << SCI;

//find maximum overlap , make it positive; store the resulting phases, use previous phases

SCI.copyonwrite();
if(aligncolumns) SCI.transposeme();
if(alignrows||aligncolumns)
	{
	NRVec<REAL> phases;
	phases.resize(SCI.nrows());

//apply old phases if available 
if(useoldphase)
{
int f=open("phases.old",O_RDONLY); 
if(f>=0)
	{
	phases.get(f,true,false);
	close(f);
	cout <<"Old phases:\n"<<phases<<endl;
	if(phases.size() != SCI.ncols() || phases.size() != SCI.nrows()) laerror("size mismatch in phases.old");
	for(int j=0; j<SCI.ncols(); ++j) for(int i=0; i<SCI.nrows(); ++i) SCI(i,j) *= phases[j];
	cout <<"CI overlap patched by old phases\n"<< (aligncolumns? SCI.transpose():SCI);
	}

}

	for(int i=0; i<SCI.nrows(); ++i)
		{
		REAL maxoverlap = -1;
		int maxoverlapindex = -1;
		for(int j=0; j<SCI.ncols(); ++j) if(abs(SCI(i,j)) > maxoverlap) 
			{
			maxoverlapindex=j;
			maxoverlap=abs(SCI(i,j));
			}

		if(i!=maxoverlapindex) cout <<"Note: state character exchanged: "<<i<<" -> "<<maxoverlapindex<<endl;
	
		if(SCI(i,i)<0)
			{	
	 		phases[i]= -1;
			for(int j=0; j<SCI.ncols(); ++j) SCI(i,j) *= -1;
			}
		else
			phases[i]= 1;
		}
	int p=open("phases",O_WRONLY|O_CREAT,0777);
	if(p<0) laerror("cannot write-open phases");
	phases.put(p,true);
	close(p);
	cout <<"Phases aligned:\n"<<phases<<endl;
	}

if(aligncolumns) SCI.transposeme();

cout << "CPU time "<< (clock()-clock0)/CLOCKS_PER_SEC<<endl;

cout << "CI overlap matrix\n";
cout << SCI;
}


