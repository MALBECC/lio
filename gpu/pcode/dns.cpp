// Local Density Functionals
// Input: cantidad de funciones s, p, d ('cant_funcs'), M ('m'), Nuc ('nuc')
// Output: densidad para el punto actual ('density')
void calc_density(double3& density, int3 cant_funcs, double m, array<double3> distancias_punto_atomos, array<double> nuc, array<int> contracciones,
									double3 posicion_punto, array<double3> posiciones_atomos)
{
//									SUBROUTINE DNS(Dens,F,Xi,ds,NORM,Nuc,ncont,nshell,a,c,r,
//     >               M,M18,NCO,RMM)
    
	/*
      implicit real*8 (a-h,o-z)
      logical NORM
      INCLUDE 'param'
      dimension c(ng,nl),a(ng,nl),Nuc(ng),ncont(ng)
      dimension r(nt,3),nshell(0:3),Xi(3)
      dimension ds(ntq),F(M),W(ng),RMM(*)
      dimension indx(ng)
c
      common /Ll/ Ll(3)
      common /Nc/ Ndens

	/*
	 * now we should evaluate all same loops as the ones used for
	 * 1 electron matrix elements, but doing only products
	 * then, the particular density functional wanted is calculated
	 */
	double fc = 1;
	if (normalize) fc = 1 / sqrt(3);
	
	density = 0;
	double density_local = 0;
	int funcs_s = cant_funcs.x;
	int funcs_p = cant_funcs.y;
	int funcs_d = cant_funcs.z;
	double m2 = 2 * m;

	// DUDA: ng?
	array<double> W(ng);
	array<double> F(m);
	
	// basis functions evaluated at r are calculated	
	for each i 1,m {
		W(i) = 0;
		F(i) = 0;
	}

	// funciones s	
	for each func 1,funcs_s {
		double dist = distancias_punto_atomos(nuc(func)); // DUDA: nuc?
		
		for each contraccion 1,contracciones(func) {
		  double rexp = factor_a(func, contraccion) * dist;

			if (rexp > 30) continue;
      F(func) += exp(-rexp) * factor_c(func,contraccion);
		}
	}
	
	// funciones p
	for each func (ns + 1),(ns + np) by 3 {
		int atomo_nucleo = nuc(func);
		double dist = distancias_punto_atomos(atomo_nucleo);
		
		for each contraccion 1,contracciones(func) {
      double rexp = factor_a(func, contraccion) * dist;
      if (rexp > 30) continue;
			
      double t = exp(-rexp) * factor_c(func, contraccion)
				
      for each dim 1,3 {
				int ii = func + dim - 1;
				F(ii) += t * (posicion_punto(dim) - (posiciones_atomos(atomo_nucleo))(dim));
			}
		}
	}
	
	// funciones d
	for each func (ns + np + 1),m by 6 {
		int atomo_nucleo = nuc(func);
		double dist = distancias_punto_atomos(atomo_nucleo);
		
		for each contraccion 1,contracciones(func) {
			double rexp = factor_a(func, contraccion) * dist;
			if (rexp > 30) continue;
			
			double t = exp(-rexp) * factor_c(func, contraccion);
			
			for each dim1 1,3 {
				double t1 = (posicion_punto(dim1) - (posiciones_atomos(atomo_nucleo))(dim1));
				
				for each dim2 1,dim1 {
					double t2 = (posicion_punto(dim2) - (posiciones_atomos(atomo_nucleo))(dim2));
					if (dim1 == dim2) t2 *= fc;
					
					double term = t * t1 * t2;
					int l12 = ???; // l12=Ll(l1)+l2
					int ii = func + l12 - 1;
					F(ii) += term;
				}				
			}			
		}		
	}
	
	//
	// calculation of vector W : density matrix scalar F
	//
	if (Ndens == 1) {
		int k = 0;
		for each j 1,m {
			if (F(j) == 0) { k += m - j + 1; continue; }
			
			for each i j,m {
				k++;
				W(j) = W(j) + RMM(k) * F(i);
			}
		}
		
		for each i 1,m {
			dens += F(i) * W(i);	
		}
	}
	else {
		float kk = M18 - 1;
		for each jj 1,NCO {
			for each ii 1,m {
				kk++;
				W(jj) += RMM(kk) * F(ii);
			}
		}
		
		// construction of sparse F vector
		// slower (at least for moderately small molecules)
		
		k = 0;
		for each i 1,m {
			if (F(i) != 0) {
				k++;
				indx(k) = i;
			}			
		}
		
		kk = M18 - 1;
		
		if (k < m / 2) {
			for each jj 1,NCO {
				W(jj) = 0;
				for each ii 1,k {
					ik = indx(ii);
					W(jj) += F(ik) * RMM(kk + ik);
				}
				kk += m;
			}
		}
		else {
			for each jj 1,NCO {
				for each ii 1,m {
					kk++;
					W(jj) += RMM(kk) * F(ii);					
				}				
			}			
		}
		
		for each ii 1,NCO {
			DENS += W(ii) * W(ii);
		}
		
		DENS *= 2;
	}
}
