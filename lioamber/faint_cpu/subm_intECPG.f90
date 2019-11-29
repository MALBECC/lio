!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INTECPG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates ECP gradients.                                                     !
!                                                                              !
! EXTERNAL INPUT: system information.                                          !
!   · natom: number of QM atoms.                                               !
!!   · ntatom: total number of atoms (QM+MM)                                    !
!!   · r(ntatom,3): atoms' coordinates.                                         !
!!   · d(natom,natom): distances between QM atoms.                              !
!!   · Iz(natom): nuclear charge for each QM atom.                              !
!!   · rho(M,M): density matrix.                                                !
!!                                                                              !
! INTERNAL INPUT: basis set information.                                       !
!   · M: number of basis functions (without contractions)                      !
!   · ncont(M): number of contractions per function.                           !
!!   · a(M,nl): basis function exponents.                                       !
!!   · c(M,nl): basis function coefficients.                                    !
!!   · nshell(0:3): number of basis functions per shell (s,p,d).                !
!   · Nuc(M): atomic index corresponding to function i.                        !
!                                                                              !
! OUTPUTS:                                                                     !
!   · ff(natom,3): ECP gradients (= -forces)                                    !
!                                                                              !
! Original : N Foglia 2019                                                     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module subm_intECPG
contains
subroutine intECPG(ff,rho,natom)
   use basis_data   , only: Nuc,M, ncont, a
   use ECP_mod, ONLY : Lxyz, ecptypes, IzECP, Cnorm, pi, ZlistECP, distx, disty, distz,VAAB,VBAC !ultmos 2 para test
   use subm_intECP   , only: AAB_LOCAL, AAB_SEMILOCAL, ABC_LOCAL, ABC_SEMILOCAL
   implicit none
   integer         , intent(in)  :: natom
   double precision, intent(out) :: ff(natom,3)
   double precision, intent(in)  :: rho(:) !why this havent got dimention ??? Nick

   ! Auxiliary variables
   integer :: i, j, k !number of basis set function
   integer :: ii, ji !number of contraction
   integer :: kecp !atoms with ECP
   integer :: lxi,lxj,lyi,lyj,lzi,lzj !l?$  potencia de la parte angular de la base
   double precision :: AAB, dAABp, dAABn
   double precision :: ABC, dABCpl, dABCnl, dABCpr, dABCnr
   double precision :: acum

   double precision :: dABC, dAAB !derivative of ABC term
   double precision :: acuml, acumr, Distcoef, dxi, dyi, dzi, dxj, dyj, dzj, dx, dy, dz
   double precision, dimension(M,M,natom) :: dHcore !para test, luego pasar a dimension 3*natot
   double precision, dimension(M,M,2,3) :: dHcore_AAB
   double precision, dimension(M,M,4) :: dHcore_AAB_temp
   double precision, dimension(M,M) :: Hcore2 !duplicated Hcore for test
   double precision, dimension(M,M) :: Hcore
   double precision :: factor
   integer :: pos
   double precision :: T1,T2,T3 !just for dbug

   dHcore=0.d0
   Hcore=0.d0
   dHcore_AAB=0.d0
   Hcore2=0.d0

   write(*,*) "TEST FFECP"
   ff=0.d0

   do i = 1, M
   do j = 1, M !cambiar luego a 1,i
      lxi=Lxyz(i,1)
      lxj=Lxyz(j,1)
      lyi=Lxyz(i,2)
      lyj=Lxyz(j,2)
      lzi=Lxyz(i,3)
      lzj=Lxyz(j,3)

	write(*,*) "doing ", i,j
!d<B|A|A> and d<A|A|B>
      if (nuc(i) .NE. nuc(j) ) THEN 
      write(*,*) "i,j", i,j
	 DO kecp=1, ecptypes ! barre atomos con ecp
	    IF (IzECP(nuc(i)) .EQ. ZlistECP(kecp) .OR. IzECP(nuc(j)) .EQ. ZlistECP(kecp)) THEN !solo calcula si el atomo tiene ECP
               dx=distx(nuc(i),nuc(j))
               dy=disty(nuc(i),nuc(j))
               dz=distz(nuc(i),nuc(j))
               Distcoef=(dx**2.d0 + dy**2.d0 + dz**2.d0)
               lxi=Lxyz(i,1)
               lxj=Lxyz(j,1)
               lyi=Lxyz(i,2)
               lyj=Lxyz(j,2)
               lzi=Lxyz(i,3)
               lzj=Lxyz(j,3) !exponentes de la parte angular

               IF (IzECP(nuc(i)) .EQ. ZlistECP(kecp)) THEN !calculo para ECP en i d<A|A|B>
                  DO ji=1, ncont(j) !barre contracciones de la base j
	             acuml=0.d0
                     acumr=0.d0
		     acum=0.d0
		     dAABp=0.d0
		     dAABn=0.d0
		     AAB=0.d0
		     DO ii=1, ncont(i) !ii barre contracciones de las funcion de base i


!local term
			AAB=  AAB_LOCAL(i,j,kecp,ii,ji,lxj  ,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)
			dAABp=AAB_LOCAL(i,j,kecp,ii,ji,lxj+1,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)
!S-local term
!			AAB=AAB+    AAB_SEMILOCAL(i,j,ii,ji,kecp,lxi,lyi,lzi,lxj  ,lyj,lzj,dx,dy,dz)
!			dAABp=dAABp+AAB_SEMILOCAL(i,j,ii,ji,kecp,lxi,lyi,lzi,lxj+1,lyj,lzj,dx,dy,dz)

			acum=acum+AAB*Cnorm(i,ii)
			acumr=acumr+dAABp*Cnorm(i,ii)*2.d0*a(j,ji)


			dAABp=0.d0
			dAABn=0.d0
			AAB=0.d0

			if (lxj.gt.0) then !p+ case
			   dAABn=AAB_LOCAL(i,j,kecp,ii,ji,lxj-1,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)
!			   dAABn=dAABn + AAB_SEMILOCAL(i,j,ii,ji,kecp,lxi,lyi,lzi,lxj-1,lyj,lzj,dx,dy,dz)
			   acuml=acuml - dAABn*Cnorm(i,ii)*lxj
			   dAABn=0.d0
			end if
		     END DO

!matriz de derivadas
			dHcore(i,j,nuc(j))= dHcore(i,j,nuc(j)) + (acumr+acuml)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))/0.529177D0
			dHcore(i,j,nuc(i))= dHcore(i,j,nuc(i)) -  (acumr+acuml)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))/0.529177D0
			Hcore(i,j)= Hcore(i,j) + acum*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))
			write(*,*) "-1-1-1", i,j,acum*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))
			acumr=0.d0
			acuml=0.d0
			acum=0.d0
!test write
			if (Hcore(i,j) .ne. 0.d0) then
			   if (dHcore(i,j,nuc(j)) .ne. 0.d0) then
			      write(*,*) "Ncheck iij", i,j, Hcore(i,j),dHcore(i,j,nuc(j)) !, -0.006999749669d0*dx*exp(-0.1241349985d0*dx**2)
			   end if
			end if

		     END DO
		  END IF

	       IF (IzECP(nuc(j)) .EQ. ZlistECP(kecp)) THEN !calculo para ECP en j d<B|A|A>
		  DO ii=1, ncont(i) ! barre contracciones de las funcion de base i
		     acuml=0.d0
		     acumr=0.d0
		     acum=0.d0
		     dAABp=0.d0
		     dAABn=0.d0
		     AAB=0.d0
		     DO ji=1, ncont(j) !barre contracciones de las funcion de base j

			AAB=  AAB_LOCAL(j,i,kecp,ji,ii,lxi  ,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)
			dAABp=AAB_LOCAL(j,i,kecp,ji,ii,lxi+1,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)

!			AAB=AAB+    AAB_SEMILOCAL(j,i,ji,ii,kecp,lxj,lyj,lzj,lxi  ,lyi,lzi,-dx,-dy,-dz)
!			dAABp=dAABp+AAB_SEMILOCAL(j,i,ji,ii,kecp,lxj,lyj,lzj,lxi+1,lyi,lzi,-dx,-dy,-dz)

			acum=acum+AAB*Cnorm(j,ji)
			acuml=acuml+dAABp*Cnorm(j,ji)*2.d0*a(i,ii)!/0.529177D0 !multiplica por el coeficiente de la base j


			dAABp=0.d0
			dAABn=0.d0
			AAB=0.d0


			if (lxi.gt.0) then !p+case
			   dAABn=AAB_LOCAL(j,i,kecp,ji,ii,lxi-1,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)
!			   dAABn=dAABn+ AAB_SEMILOCAL(j,i,ji,ii,kecp,lxj,lyj,lzj,lxi-1,lyi,lzi,-dx,-dy,-dz)
			   acumr=acumr - dAABn*Cnorm(j,ji)*lxi
			   dAABn=0.d0
			end if


		     END DO

		     dHcore(i,j,nuc(i))=dHcore(i,j,nuc(i)) + (acumr+acuml)*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii))/0.529177D0
		     dHcore(i,j,nuc(j))=dHcore(i,j,nuc(j)) - (acumr+acuml)*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii))/0.529177D0
		     Hcore(i,j)=Hcore(i,j)+ acum*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii))

                        acumr=0.d0
                        acuml=0.d0
                        acum=0.d0

!a		     if (Hcore(i,j).ne.0.d0) then
!		       if (dHcore(i,j,nuc(j)).ne.0.d0) then
!			  write(*,*) "Ncheckijj", i,j, Hcore(i,j),dHcore(i,j,nuc(j))
!		       endif
!		     endif

		  END DO
	       END IF

	    END IF
	 END DO

	else !<A|B|A> and <A|B|C> derivatives

	   DO k=1, natom !barre todos los nucleoas del sistema
	      if (nuc(i) .NE. k .AND. nuc(j) .NE. k) THEN !solo calcula si las 2 funciones de base NO corresponden al atomo con el ECP
		do kecp=1, ecptypes !barre atomos con ecp
		   if (IzECP(k) .EQ. ZlistECP(kecp).and..false.) THEN !solo calcula si el nucleo tiene ecp
		      dxi=-distx(nuc(i),k)
		      dyi=-disty(nuc(i),k)
		      dzi=-distz(nuc(i),k)
		      dxj=-distx(nuc(j),k)
		      dyj=-disty(nuc(j),k)
		      dzj=-distz(nuc(j),k) ! distancias al nucleo con ecp
		      lxi=Lxyz(i,1)
		      lxj=Lxyz(j,1)
		      lyi=Lxyz(i,2)
		      lyj=Lxyz(j,2)
		      lzi=Lxyz(i,3)
		      lzj=Lxyz(j,3)

		      acum=0.d0
		      acuml=0.d0
		      acumr=0.d0

		      ABC=0.d0
		      dABCpl=0.d0
		      dABCpr=0.d0
		
		      DO ii=1, ncont(i) !barre contracciones de la base i
			 DO ji=1, ncont(j) !barre contracciones de la base j
			    Distcoef=a(i,ii)*(dxi**2.d0 + dyi**2.d0 + dzi**2.d0) + a(j,ji)*(dxj**2.d0 + dyj**2.d0 + dzj**2.d0)

			    ABC=ABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj) 
			    dABCpl=ABC_LOCAL(i,j,ii,ji,k,lxi+1,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
			    dABCpr=ABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj+1,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)

			    ABC=ABC +     4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi,  lyi,lzi,lxj,  lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
			    dABCpl=dABCpl+4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi+1,lyi,lzi,lxj,  lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
			    dABCpr=dABCpr+4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi,  lyi,lzi,lxj+1,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)

			    acum=acum +     ABC *Cnorm(j,ji)*exp(-Distcoef)!*Cnorm(i,ii)
			    ABC=0.d0

			    acuml=acuml + dABCpl*Cnorm(j,ji)*exp(-Distcoef)*2.d0*a(i,ii)!*Cnorm(i,ii)
			    dABCpl=0.d0

			    acumr=acumr + dABCpr*Cnorm(j,ji)*exp(-Distcoef)*2.d0*a(j,ji)!*Cnorm(i,ii) !agregue Cnorm para test, sacar
			    dABCpr=0.d0


			    if (lxi.gt.0) then !p+ case
			      dABCnl=0.d0
			      dABCnl=ABC_LOCAL(i,j,ii,ji,k,lxi-1,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
			      dABCnl=dABCnl+4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi-1,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
			      acuml=acuml - dABCnl*Cnorm(j,ji)*exp(-Distcoef)*lxi
			      dABCnl=0.d0
			    end if

			    if (lxj.gt.0) then !p+ case
			      dABCnr=0.d0
			      dABCnr=ABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj-1,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
			      dABCnr=dABCnr+4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj-1,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
			      acumr=acumr - dABCnr*Cnorm(j,ji)*exp(-Distcoef)*lxj
			      dABCnr=0.d0
			    end if

			 END DO

!			 Hcore(i,j) = Hcore(i,j) + acum*Cnorm(i,ii)*4.d0*pi
!			 dHcore(i,j,nuc(i))=dHcore(i,j,nuc(i)) + acuml*Cnorm(i,ii)*4.d0*pi/0.529177D0
!			 dHcore(i,j,nuc(j))=dHcore(i,j,nuc(j)) + acumr*Cnorm(i,ii)*4.d0*pi/0.529177D0 
!			 dHcore(i,j,k)=dHcore(i,j,k) -   (acuml+acumr)*Cnorm(i,ii)*4.d0*pi/0.529177D0

!			 T1=acuml*Cnorm(i,ii)*4.d0*pi/0.529177D0
!			 T2=acumr*Cnorm(i,ii)*4.d0*pi/0.529177D0
!			 T3=-(acuml+acumr)*Cnorm(i,ii)*4.d0*pi/0.529177D0

!                          write(*,*) "Ncheckikj", i,k,j, Hcore(i,j),dHcore(i,j,nuc(i)),dHcore(i,j,nuc(j))

			 acum=0.d0
			 acuml=0.d0
			 acumr=0.d0
		      END DO
		   END IF
		END DO
	     END IF
	   END DO
	end if
     end do
  end do


!simetry test
   do i=1,M
     do j=1,M
       do k=1, natom
	     write(123456,*) "error in dHcore(i,j,k)",i,j,k, dHcore(i,j,k)-dHcore(j,i,k),dHcore(i,j,k)
       end do
     end do
   end do

!forces
   do i=1,M
     do j=1,i
	factor=1.d0
        pos=i+(1-j)*(j-2*M)/2
        do k=1, natom
	  ff(k,1)= ff(k,1) + factor*dHcore(i,j,k)*rho(pos)
        end do
     end do
   end do



!#######################################################################3
!New loops making derivativs at same time than fock elements


   do i = 1, M
   do j = 1, M !cambiar luego a 1,i
      lxi=Lxyz(i,1)
      lxj=Lxyz(j,1)
      lyi=Lxyz(i,2)
      lyj=Lxyz(j,2)
      lzi=Lxyz(i,3)
      lzj=Lxyz(j,3)

!        write(*,*) "doing ", i,j
!d<B|A|A> and d<A|A|B>
      if (nuc(i) .NE. nuc(j) ) THEN
!      write(*,*) "i,j", i,j
         DO kecp=1, ecptypes ! barre atomos con ecp
            IF (IzECP(nuc(i)) .EQ. ZlistECP(kecp) .OR. IzECP(nuc(j)) .EQ. ZlistECP(kecp)) THEN !solo calcula si el atomo tiene ECP
               dx=distx(nuc(i),nuc(j))
               dy=disty(nuc(i),nuc(j))
               dz=distz(nuc(i),nuc(j))
               Distcoef=(dx**2.d0 + dy**2.d0 + dz**2.d0)
               lxi=Lxyz(i,1)
               lxj=Lxyz(j,1)
               lyi=Lxyz(i,2)
               lyj=Lxyz(j,2)
               lzi=Lxyz(i,3)
               lzj=Lxyz(j,3) !exponentes de la parte angular
		
               IF (IzECP(nuc(i)) .EQ. ZlistECP(kecp)) THEN !calculo para ECP en i d<A|A|B>
                  DO ji=1, ncont(j) !barre contracciones de la base j
                     dHcore_AAB_temp=0.d0
                     DO ii=1, ncont(i) !ii barre contracciones de las funcion de base i
			dHcore_AAB_temp(i,j,1:4)=dHcore_AAB_temp(i,j,1:4)+dAAB_LOCAL(i,j,kecp,ii,ji,lxj,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)*Cnorm(i,ii)
!			dHcore_AAB_temp(i,j,1:4)=dHcore_AAB_temp(i,j,1:4)+dAAB_SEMILOCAL(i,j,ii,ji,kecp,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)*Cnorm(i,ii)
                     END DO
!matriz de derivadasa
			Hcore2(i,j)= Hcore2(i,j) + dHcore_AAB_temp(i,j,1)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))
			write(*,*) "000",i,j,dHcore_AAB_temp(i,j,1)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))
			dHcore_AAB(i,j,1,1:3)=dHcore_AAB(i,j,1,1:3)-dHcore_AAB_temp(i,j,2:4)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))/0.529177D0 !di/d
			dHcore_AAB(i,j,2,1:3)=dHcore_AAB(i,j,2,1:3)+dHcore_AAB_temp(i,j,2:4)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))/0.529177D0 !dj/d
			dHcore_AAB_temp=0.d0
                  END DO
               END IF


               IF (IzECP(nuc(j)) .EQ. ZlistECP(kecp)) THEN !calculo para ECP en j d<B|A|A>
                  DO ii=1, ncont(i) ! barre contracciones de las funcion de base i
		     dHcore_AAB_temp=0.d0
                     DO ji=1, ncont(j) !barre contracciones de las funcion de base j
			dHcore_AAB_temp(i,j,1:4)=dHcore_AAB_temp(i,j,1:4)+ dAAB_LOCAL(j,i,kecp,ji,ii,lxi,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)*Cnorm(j,ji)
!			dHcore_AAB_temp(i,j,1:4)=dHcore_AAB_temp(i,j,1:4)+ dAAB_SEMILOCAL(j,i,ji,ii,kecp,lxj,lyj,lzj,lxi,lyi,lzi,-dx,-dy,-dz)*Cnorm(j,ji)
                     END DO
!matriz de derivadasa
                        Hcore2(i,j)= Hcore2(i,j) + dHcore_AAB_temp(i,j,1)*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii))
			dHcore_AAB(i,j,1,1:3)=dHcore_AAB(i,j,1,1:3)+dHcore_AAB_temp(i,j,2:4)*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii))/0.529177D0
			dHcore_AAB(i,j,2,1:3)=dHcore_AAB(i,j,2,1:3)-dHcore_AAB_temp(i,j,2:4)*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii))/0.529177D0
			dHcore_AAB_temp=0.d0
                  END DO
               END IF

            END IF
         END DO

	else !<A|B|A> and <A|B|C> derivatives
        end if
     end do
  end do


   do i=1,M
     do j=1,M
        if (j.le.i) then
          pos=i+(1-j)*(j-2*M)/2

	write(9876,*) "i,j,H",i,j,Hcore(i,j), VAAB(pos)+VBAC(pos), dHcore(i,j,:),dHcore_AAB(i,j,1,1),dHcore_AAB(i,j,2,1)

	  if( (Hcore(i,j).ne.0.d0) .and. (Hcore2(i,j).ne.0.d0)) then
		if ((dHcore_AAB(i,j,1,1).ne.0.d0) .and.(dHcore_AAB(i,j,2,1).ne.0.d0)) &
            write(9877,*) "i,j,H",i,j,Hcore(i,j), Hcore2(i,j),dHcore_AAB(i,j,1,1),dHcore_AAB(i,j,2,1)

		if ((dHcore_AAB(i,j,1,2).ne.0.d0) .and. (dHcore_AAB(i,j,2,2).ne.0.d0)) &
            write(9878,*) "i,j,H",i,j,Hcore(i,j), Hcore2(i,j),dHcore_AAB(i,j,1,2),dHcore_AAB(i,j,2,2)

		if ((dHcore_AAB(i,j,1,3).ne.0.d0) .and. (dHcore_AAB(i,j,2,3).ne.0.d0)) &
            write(9879,*) "i,j,H",i,j,Hcore(i,j), Hcore2(i,j),dHcore_AAB(i,j,1,3),dHcore_AAB(i,j,2,3)
	  end if
        endif

        do k=1, natom
          write(5300,*) i,j,k,pos,dHcore(i,j,k),rho(pos)
        end do
     end do
   end do





   return;
END SUBROUTINE intECPG


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!WIP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION dAAB_LOCAL(i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz)
!Calcula el termino local del pseudopotencial centrado en i [<xi|Vi(LM)|xj>] y sus derivadas 
!los coef de la base se multiplican en la rutina que llama a esta
!i,j funcion de base
!ii,ji numero de contraccion de la funcion de base
!k atomo con ECP
USE basis_data, ONLY : a
USE ECP_mod, ONLY :nECP,bECP, aECP, ZlistECP, Lmax, expnumbersECP, Qnl, angularint,Fulltimer_ECP,tlocal,tQ1
use subm_intECP   , only: OMEGA1, comb, qtype1n
! Vl= Σ aECP * r^nECP * exp(-bECP r^2)
! ZlistECP(k) carga del atomo k con ECP
! Lmax(Z) maximo momento angular del pseudopotencial para el atomo con carga nuclear Z
! expnumbersECP(Z,l) cantidad de terminos del ECP para el atomo con carga nuclear Z y l del ECP
!              ͚ 
! Qnl(n,l) = ʃ Ml(k*r)*r^n * exp(-cr^2) dr
!            ̊ 

! Fulltimer_ECP activa los timers para int. radiales
! semilocal,tQ1 auxiliares para el calculo de tiempos

! angularint(i,j,k)= ʃ(x/r)^i (y/r)^j (z/r)^k dΩ
! tlocal,tQ1 auxiliares para timers
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi
   DOUBLE PRECISION, DIMENSION(4) :: dAAB_LOCAL !(<A|A|B>,<A|A|dB/dxB>,<A|A|dB/dyB>,<A|A|dB/dzB>
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz
   INTEGER :: z,l,Lmaxbase
! Z carga nuclear
! l maximo valor del momento angular del pseudopotencial
   INTEGER :: w,lxi,lyi,lzi,lambda !Auxiliares
   INTEGER :: limitx, limity, limitz
   INTEGER :: DERIV
   DOUBLE PRECISION :: FACTOR
   DOUBLE PRECISION :: Ccoef, Kmod,distcoefx, distcoefy,distcoefz, integral
   DOUBLE PRECISION, DIMENSION (4) :: acum
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   DOUBLE PRECISION :: t1,t2,t1q,t2q !auxiliares para timers

   limitx=-1
   limity=-1
   limitz=-1
   if(lx.eq.0) limitx=0
   if(ly.eq.0) limity=0
   if(lz.eq.0) limitz=0

   IF (Fulltimer_ECP) CALL cpu_time ( t1 )

   Z=ZlistECP(k)
   L=Lmax(Z)
   Lmaxbase=lx+ly+lz+kxi+kyi+kzi+1
   dAAB_LOCAL=0.d0
   acum=0.d0
   integral=0.d0
   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod= 2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)
!angular integrals should be calculated here

   DO w =1, expnumbersECP(z,l) !barre todos los terminos del Lmaximo
      Qnl=0.d0
      Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
      IF (Fulltimer_ECP) CALL cpu_time ( t1q )
      CALL Qtype1N(Kmod,Ccoef,Lmaxbase,necp(Z,l,w)+Lmaxbase,necp(Z,l,w)+kxi+kyi+kzi) !calcula integrales radiales

      IF (Fulltimer_ECP) THEN
         CALL cpu_time ( t2q )
         tQ1=tQ1 +t2q-t1q
      END IF

!energy
      DO lxi=0,lx !barre potencias por expansion del binomio de Newton (x - dx)^lx
	distcoefx=dx**(lx-lxi)
         DO lyi=0,ly
	  distcoefy=dy**(ly-lyi)
          DO lzi=0,lz
	    distcoefz=dz**(lz-lzi)
	    DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
		integral=integral + OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi) * &
		Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda)
	    END DO
		acum(1)= acum(1) + integral*distcoefx * distcoefy * distcoefz *comb(lx,lxi) *comb(ly,lyi) * comb(lz,lzi)
	    integral=0.d0
            END DO
         END DO
      END DO
!d/dx 
      do DERIV=1,limitx,-2
      DO lxi=0,lx+DERIV !barre potencias por expansion del binomio de Newton (x - dx)^lx
        distcoefx=dx**(lx+DERIV-lxi)
         DO lyi=0,ly
          distcoefy=dy**(ly-lyi)
          DO lzi=0,lz
            distcoefz=dz**(lz-lzi)
            DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
                integral=integral + OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi) * &
                Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda)
            END DO
		FACTOR=-dble(lx)
                if (DERIV.eq.1) FACTOR=2.d0*a(j,ji)
		acum(2)= acum(2) + integral*distcoefx * distcoefy * distcoefz *comb(lx+DERIV,lxi) *comb(ly,lyi) * comb(lz,lzi)*FACTOR
            integral=0.d0
            END DO
         END DO
      END DO
      enddo

!d/dy
      do DERIV=1,limity,-2
      DO lxi=0,lx !barre potencias por expansion del binomio de Newton (x - dx)^lx
        distcoefx=dx**(lx-lxi)
         DO lyi=0,ly+DERIV
          distcoefy=dy**(ly+DERIV-lyi)
          DO lzi=0,lz
            distcoefz=dz**(lz-lzi)
            DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
                integral=integral + OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi) * &
                Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda)
            END DO
                FACTOR=-dble(ly)
                if (DERIV.eq.1) FACTOR=2.d0*a(j,ji)
		acum(3)= acum(3) + integral*distcoefx * distcoefy * distcoefz *comb(lx,lxi) *comb(ly+DERIV,lyi) * comb(lz,lzi)*FACTOR
            integral=0.d0
            END DO
         END DO
      END DO
      end do
!d/dz
      do DERIV=1,limitz,-2
      DO lxi=0,lx !barre potencias por expansion del binomio de Newton (x - dx)^lx
        distcoefx=dx**(lx-lxi)
         DO lyi=0,ly
          distcoefy=dy**(ly-lyi)
          DO lzi=0,lz+DERIV
            distcoefz=dz**(lz-lzi)
            DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
                integral=integral + OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi) * &
                Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda)
            END DO
                FACTOR=-dble(lz)
                if (DERIV.eq.1) FACTOR=2.d0*a(j,ji)
		acum(4)= acum(4) + integral*distcoefx * distcoefy * distcoefz *comb(lx,lxi) *comb(ly,lyi) * comb(lz+DERIV,lzi)*FACTOR
            integral=0.d0
            END DO
         END DO
      END DO
      enddo

     dAAB_LOCAL=dAAB_LOCAL+aECP(z,L,w)*acum 
     acum=0.d0
   END DO

   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tlocal=tlocal+t2-t1
   END IF

   RETURN
END FUNCTION dAAB_LOCAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end module subm_intECPG

