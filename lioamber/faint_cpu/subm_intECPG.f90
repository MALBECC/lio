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
subroutine intECPG(natom)
   use basis_data   , only: Nuc,M, ncont, a
   use ECP_mod, ONLY : Lxyz, ecptypes, IzECP, Cnorm, pi, ZlistECP, distx, disty, distz,ECPatoms_order,ECPatoms, &
   dHcore_AAB, dHcore_ABC,VAAB,VBAC !ultmos 2 para test
   use subm_intECP   , only: AAB_LOCAL, AAB_SEMILOCAL, ABC_LOCAL, ABC_SEMILOCAL
   implicit none
   integer         , intent(in)  :: natom
!   double precision, intent(in)  :: rho(:) !why this havent got dimention ??? Nick

   ! Auxiliary variables
   integer :: i, j, k !number of basis set function
   integer :: ii, ji !number of contraction
   integer :: kecp !atoms with ECP
   integer :: lxi,lxj,lyi,lyj,lzi,lzj !l?$  potencia de la parte angular de la base
   double precision :: AAB, dAABp, dAABn
   double precision :: ABC, dABCpl, dABCnl, dABCpr, dABCnr
   double precision :: acum

   double precision :: acuml, acumr, Distcoef, dxi, dyi, dzi, dxj, dyj, dzj, dx, dy, dz
   double precision, dimension(M,M,natom) :: dHcore !para test, luego pasar a dimension 3*natot
   double precision, dimension(4) :: dHcore_AAB_temp !se tiene q bajar a dimension 4
   double precision, dimension(7) :: dHcore_ABC_temp
   double precision, dimension(M,M) :: Hcore2 !duplicated Hcore for test
   double precision, dimension(M,M) :: Hcore
   integer :: pos
   double precision, dimension(natom) :: F_i !just for test
   integer :: l

   dHcore=0.d0
   Hcore=0.d0
   dHcore_AAB=0.d0
   dHcore_ABC=0.d0
   dHcore_ABC_temp=0.d0
   Hcore2=0.d0

	if (.false.) then !test only remove this part at the end
   do i = 1, M
   do j = 1, M !cambiar luego a 1,i
      lxi=Lxyz(i,1)
      lxj=Lxyz(j,1)
      lyi=Lxyz(i,2)
      lyj=Lxyz(j,2)
      lzi=Lxyz(i,3)
      lzj=Lxyz(j,3)

!d<B|A|A> and d<A|A|B>
      if (nuc(i) .NE. nuc(j) ) THEN 
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
                        AAB=  AAB_LOCAL(i,j,kecp,ii,ji,lxj  ,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)!local term
                        dAABp=AAB_LOCAL(i,j,kecp,ii,ji,lxj+1,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)!local term
                        AAB=AAB+    AAB_SEMILOCAL(i,j,ii,ji,kecp,lxi,lyi,lzi,lxj  ,lyj,lzj,dx,dy,dz)!S-local term
                        dAABp=dAABp+AAB_SEMILOCAL(i,j,ii,ji,kecp,lxi,lyi,lzi,lxj+1,lyj,lzj,dx,dy,dz)!S-local term

                        acum=acum+AAB*Cnorm(i,ii)
                        acumr=acumr+dAABp*Cnorm(i,ii)*2.d0*a(j,ji)

                        dAABp=0.d0
                        dAABn=0.d0
                        AAB=0.d0

                        if (lxj.gt.0) then !p+ case
                           dAABn=AAB_LOCAL(i,j,kecp,ii,ji,lxj-1,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)
                           dAABn=dAABn + AAB_SEMILOCAL(i,j,ii,ji,kecp,lxi,lyi,lzi,lxj-1,lyj,lzj,dx,dy,dz)
                           acuml=acuml - dAABn*Cnorm(i,ii)*lxj
                           dAABn=0.d0
                        end if
                     END DO

!matriz de derivadas
   dHcore(i,j,nuc(j))= dHcore(i,j,nuc(j)) + (acumr+acuml)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))/0.529177D0
   dHcore(i,j,nuc(i))= dHcore(i,j,nuc(i)) -  (acumr+acuml)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))/0.529177D0
   Hcore(i,j)= Hcore(i,j) + acum*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))
   acumr=0.d0
   acuml=0.d0
   acum=0.d0

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

                        AAB=AAB+    AAB_SEMILOCAL(j,i,ji,ii,kecp,lxj,lyj,lzj,lxi  ,lyi,lzi,-dx,-dy,-dz)
                        dAABp=dAABp+AAB_SEMILOCAL(j,i,ji,ii,kecp,lxj,lyj,lzj,lxi+1,lyi,lzi,-dx,-dy,-dz)

                        acum=acum+AAB*Cnorm(j,ji)
                        acuml=acuml+dAABp*Cnorm(j,ji)*2.d0*a(i,ii)!/0.529177D0 !multiplica por el coeficiente de la base j


                        dAABp=0.d0
                        dAABn=0.d0
                        AAB=0.d0

                        if (lxi.gt.0) then !p+case
                           dAABn=AAB_LOCAL(j,i,kecp,ji,ii,lxi-1,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)
                           dAABn=dAABn+ AAB_SEMILOCAL(j,i,ji,ii,kecp,lxj,lyj,lzj,lxi-1,lyi,lzi,-dx,-dy,-dz)
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
                  END DO
               END IF
            END IF
         END DO

      else !<A|B|A> and <A|B|C> derivatives

         DO k=1, natom !barre todos los nucleoas del sistema
            if (nuc(i) .NE. k .AND. nuc(j) .NE. k) THEN !solo calcula si las 2 funciones de base NO corresponden al atomo con el ECP
               do kecp=1, ecptypes !barre atomos con ecp
                  if (IzECP(k) .EQ. ZlistECP(kecp)) THEN !solo calcula si el nucleo tiene ecp
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

   ABC=ABC + 4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
   dABCpl=dABCpl+4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi+1,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
   dABCpr=dABCpr+4.d0*pi*ABC_SEMILOCAL(i,j,ii,ji,k,lxi,  lyi,lzi,lxj+1,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)

   acum=acum + ABC *Cnorm(j,ji)*exp(-Distcoef)!*Cnorm(i,ii)
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

 Hcore(i,j) = Hcore(i,j) + acum*Cnorm(i,ii)*4.d0*pi
 dHcore(i,j,nuc(i))=dHcore(i,j,nuc(i)) + acuml*Cnorm(i,ii)*4.d0*pi/0.529177D0
 dHcore(i,j,nuc(j))=dHcore(i,j,nuc(j)) + acumr*Cnorm(i,ii)*4.d0*pi/0.529177D0 
 dHcore(i,j,k)=dHcore(i,j,k) -   (acuml+acumr)*Cnorm(i,ii)*4.d0*pi/0.529177D0

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

end if !if false
!#######################################################################3!#######################################################################3
!#######################################################################3!#######################################################################3
!#######################################################################3!#######################################################################3
!#######################################################################3!#######################################################################3
!#######################################################################3!#######################################################################3
!#######################################################################3!#######################################################################3
!New loops making derivativs at same time than fock elements


   do i = 1, M
   do j = 1, M !cambiar luego a 1,i
      lxi=Lxyz(i,1)
      lxj=Lxyz(j,1)
      lyi=Lxyz(i,2)
      lyj=Lxyz(j,2)
      lzi=Lxyz(i,3)
      lzj=Lxyz(j,3)

      if (nuc(i) .NE. nuc(j) ) THEN !d<B|A|A>/dx,y,z and d<A|A|B>/dx,y,z

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
   dHcore_AAB_temp(1:4)=dHcore_AAB_temp(1:4)+dAAB_LOCAL(i,j,kecp,ii,ji,lxj,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)*Cnorm(i,ii)
   dHcore_AAB_temp(1:4)=dHcore_AAB_temp(1:4)+dAAB_SEMILOCAL(i,j,ii,ji,kecp,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)*Cnorm(i,ii)
                     END DO

   Hcore2(i,j)= Hcore2(i,j) + dHcore_AAB_temp(1)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))
   dHcore_AAB(i,j,1,1:3)=dHcore_AAB(i,j,1,1:3)-dHcore_AAB_temp(2:4)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))/0.529177D0 !di/d
   dHcore_AAB(i,j,2,1:3)=dHcore_AAB(i,j,2,1:3)+dHcore_AAB_temp(2:4)*Cnorm(j,ji)*4.d0*pi*exp(-Distcoef*a(j,ji))/0.529177D0 !dj/d

   dHcore_AAB_temp=0.d0
                  END DO
               END IF


               IF (IzECP(nuc(j)) .EQ. ZlistECP(kecp)) THEN !calculo para ECP en j d<B|A|A>
                  DO ii=1, ncont(i) ! barre contracciones de las funcion de base i
                     dHcore_AAB_temp=0.d0

                     DO ji=1, ncont(j) !barre contracciones de las funcion de base j
   dHcore_AAB_temp(1:4)=dHcore_AAB_temp(1:4)+ dAAB_LOCAL(j,i,kecp,ji,ii,lxi,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)*Cnorm(j,ji)
   dHcore_AAB_temp(1:4)=dHcore_AAB_temp(1:4)+ dAAB_SEMILOCAL(j,i,ji,ii,kecp,lxj,lyj,lzj,lxi,lyi,lzi,-dx,-dy,-dz)*Cnorm(j,ji)
                     END DO

   Hcore2(i,j)= Hcore2(i,j) + dHcore_AAB_temp(1)*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii))
   dHcore_AAB(i,j,1,1:3)=dHcore_AAB(i,j,1,1:3)+dHcore_AAB_temp(2:4)*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii))/0.529177D0
   dHcore_AAB(i,j,2,1:3)=dHcore_AAB(i,j,2,1:3)-dHcore_AAB_temp(2:4)*Cnorm(i,ii)*4.d0*pi*exp(-Distcoef*a(i,ii))/0.529177D0

   dHcore_AAB_temp=0.d0
                  END DO
               END IF

            END IF
         END DO
      else !<A|B|A> and <A|B|C> derivatives


         DO k=1, natom !barre todos los nucleoas del sistema
            if (nuc(i) .NE. k .AND. nuc(j) .NE. k) THEN !solo calcula si las 2 funciones de base NO corresponden al atomo con el ECP
               do kecp=1, ecptypes !barre atomos con ecp
                  if (IzECP(k) .EQ. ZlistECP(kecp)) THEN !solo calcula si el nucleo tiene ecp
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
                     dHcore_ABC_temp=0.d0

                     DO ii=1, ncont(i) !barre contracciones de la base i
                        dHcore_ABC_temp=0.d0
                        DO ji=1, ncont(j) !barre contracciones de la base j
                           Distcoef=a(i,ii)*(dxi**2.d0 + dyi**2.d0 + dzi**2.d0) + a(j,ji)*(dxj**2.d0 + dyj**2.d0 + dzj**2.d0)

   dHcore_ABC_temp=dHcore_ABC_temp+ &
   dABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)*Cnorm(j,ji)*exp(-Distcoef)
   dHcore_ABC_temp=dHcore_ABC_temp+ &
   4.d0*pi*dABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)*Cnorm(j,ji)*exp(-Distcoef)
                         END DO

   Hcore2(i,j) = Hcore2(i,j) + dHcore_ABC_temp(1)*Cnorm(i,ii)*4.d0*pi
   dHcore_ABC(i,j,1,1:3)=dHcore_ABC(i,j,1,1:3)+dHcore_ABC_temp(2:4)*Cnorm(i,ii)*4.d0*pi/0.529177D0
   dHcore_ABC(i,j,2,1:3)=dHcore_ABC(i,j,2,1:3)+dHcore_ABC_temp(5:7)*Cnorm(i,ii)*4.d0*pi/0.529177D0
   dHcore_ABC(i,j,2+ECPatoms_order(k),1:3)=dHcore_ABC(i,j,2+ECPatoms_order(k),1:3)- &
   dHcore_ABC_temp(2:4)*Cnorm(i,ii)*4.d0*pi/0.529177D0-dHcore_ABC_temp(5:7)*Cnorm(i,ii)*4.d0*pi/0.529177D0
   dHcore_ABC_temp=0.d0

                      END DO
                  END IF
               END DO
            END IF
         END DO
      end if
   end do
   end do


   do i=1,M
     do j=1,M
        if (j.le.i) then
          pos=i+(1-j)*(j-2*M)/2
          write(9876,*) "i,j,H",i,j,Hcore(i,j), VAAB(pos)+VBAC(pos),Hcore2(i,j)
        endif
     end do
   end do

!NN CHECK

   do i=1,M
      do j=1,M
         do l=1,3
             if (dHcore_AAB(i,j,1,l).ne.dHcore_AAB(i,j,1,l)) then
                write(*,*) "NAN en: dHcore_AAB", i,j,"1",l
                stop
             elseif (dHcore_AAB(i,j,2,l).ne.dHcore_AAB(i,j,2,l)) then
                write(*,*) "NAN en: dHcore_AAB", i,j,"2",l
                stop
             end if
          end do
       end do
    end do

   do i=1,M
      do j=1,M
         do k=1,2+ECPatoms
            do l=1,3
               if (dHcore_ABC(i,j,k,l).ne.dHcore_ABC(i,j,k,l)) then
                write(*,*) "NAN en: dHcore_ABC", i,j,k,l
                stop
               end if
            end do
          end do
       end do
    end do


!este loop computa las fuerzas multiplicando por rho
!   ff=0.d0
   do i=1,M
      do j=1,M
         do l=1,3
            F_i=0.d0
            do k=1,natom
               if (k.eq.nuc(i)) F_i(k)=F_i(k)+(dHcore_AAB(i,j,1,l)+dHcore_ABC(i,j,1,l))*0.529177D0
               if (k.eq.nuc(j)) F_i(k)=F_i(k)+(dHcore_AAB(i,j,2,l)+dHcore_ABC(i,j,2,l))*0.529177D0
               do kecp=1, ecptypes !barre atomos con ecp
                  if (IzECP(k) .EQ. ZlistECP(kecp)) THEN !solo calcula si el nucleo tiene ecp
                     F_i(k)=F_i(k)+dHcore_ABC(i,j,2+ECPatoms_order(k),l)*0.529177D0
                  end if
               end do

!               if (j.le.i) then
!                  pos=i+(1-j)*(j-2*M)/2
!                  ff(k,l)=ff(k,l)+F_i(k)*rho(pos)
!               end if
            enddo
            if (j.le.i) then
               pos=i+(1-j)*(j-2*M)/2
               if(Hcore2(i,j).ne.0.d0) then
                  if(l.eq.1) write(10000,*) "i,j,H",i,j,Hcore2(i,j), VAAB(pos)+VBAC(pos), F_i
               end  if
            end if
            if (Hcore(i,j).ne.0.d0 ) write(9900+l,*) "i,j,H",i,j,Hcore(i,j),F_i
         end do
      end do
   enddo
   return;
END SUBROUTINE intECPG

SUBROUTINE ECP_gradients(ff,rho,natom)
   use basis_data   , only: M, Nuc
   use ECP_mod, ONLY : dHcore_AAB, dHcore_ABC, ecptypes, IzECP, ZlistECP, ECPatoms_order
   IMPLICIT NONE
   integer         , intent(in)  :: natom
   double precision, intent(out) :: ff(natom,3)
   double precision, intent(in)  :: rho(:) !why this havent got dimention ??? Nick
   integer :: i,j,l,k, kecp, pos
   double precision, dimension(natom) :: F_i
   ff=0.d0
   do i=1,M
      do j=1,M
         do l=1,3
            F_i=0.d0
            do k=1,natom
               if (k.eq.nuc(i)) F_i(k)=F_i(k)+(dHcore_AAB(i,j,1,l)+dHcore_ABC(i,j,1,l))*0.529177D0
               if (k.eq.nuc(j)) F_i(k)=F_i(k)+(dHcore_AAB(i,j,2,l)+dHcore_ABC(i,j,2,l))*0.529177D0
               do kecp=1, ecptypes !barre atomos con ecp
                  if (IzECP(k) .EQ. ZlistECP(kecp)) THEN !solo calcula si el nucleo tiene ecp
                     F_i(k)=F_i(k)+dHcore_ABC(i,j,2+ECPatoms_order(k),l)*0.529177D0
                  end if
               end do

               if (j.le.i) then
                  pos=i+(1-j)*(j-2*M)/2
                  ff(k,l)=ff(k,l)+F_i(k)*rho(pos)
               end if
            enddo
         end do
      end do
   enddo
   return
END SUBROUTINE ECP_gradients

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION dAAB_LOCAL(i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz)
!Calcula el termino local del pseudopotencial centrado en i [<xi|Vi(LM)|xj>] y sus derivadas 
!los coef de la base se multiplican en la rutina que llama a esta
!i,j funcion de base
!ii,ji numero de contraccion de la funcion de base
!k atomo con ECP
USE basis_data, ONLY : a
USE ECP_mod, ONLY :nECP,bECP, aECP, ZlistECP, Lmax, expnumbersECP, Qnl, Fulltimer_ECP,tlocal,tQ1
use subm_intECP   , only: OMEGA1, comb, qtype1n, Anal_radial_int
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
   DOUBLE PRECISION :: Ccoef, Kmod,distcoefx, distcoefy,distcoefz, integral
   DOUBLE PRECISION, DIMENSION (4) :: acum
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   DOUBLE PRECISION :: t1,t2,t1q,t2q !auxiliares para timers

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
      CALL Anal_radial_int(1)

      IF (Fulltimer_ECP) THEN
         CALL cpu_time ( t2q )
         tQ1=tQ1 +t2q-t1q
      END IF

!Fock
      acum(1)= acum(1) + AAB_LOCAL_loops(i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz,w)
!d/dx
      acum(2)= acum(2) + AAB_LOCAL_loops(i,j,k,ii,ji,lx+1,ly,lz,kxi,kyi,kzi,dx,dy,dz,w)*2.d0*a(j,ji)
      if (lx .gt.0) acum(2)= acum(2) - AAB_LOCAL_loops(i,j,k,ii,ji,lx-1,ly,lz,kxi,kyi,kzi,dx,dy,dz,w)*dble(lx)
!d/dy
      acum(3)= acum(3) + AAB_LOCAL_loops(i,j,k,ii,ji,lx,ly+1,lz,kxi,kyi,kzi,dx,dy,dz,w)*2.d0*a(j,ji)
      if (ly .gt.0) acum(3)= acum(3) - AAB_LOCAL_loops(i,j,k,ii,ji,lx,ly-1,lz,kxi,kyi,kzi,dx,dy,dz,w)*dble(ly)
!d/dz
      acum(4)= acum(4) + AAB_LOCAL_loops(i,j,k,ii,ji,lx,ly,lz+1,kxi,kyi,kzi,dx,dy,dz,w)*2.d0*a(j,ji)
      if (lz .gt.0) acum(4)= acum(4) - AAB_LOCAL_loops(i,j,k,ii,ji,lx,ly,lz-1,kxi,kyi,kzi,dx,dy,dz,w)*dble(lz)

      dAAB_LOCAL=dAAB_LOCAL+aECP(z,L,w)*acum 
      acum=0.d0
   END DO

   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tlocal=tlocal+t2-t1
   END IF

   RETURN
END FUNCTION dAAB_LOCAL





DOUBLE PRECISION FUNCTION AAB_LOCAL_loops(i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz,w)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : Qnl,necp, ZlistECP,Lmax
   use subm_intECP   , only: comb, OMEGA1
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi,w
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz
   INTEGER :: z,l
! Z carga nuclear
! l maximo valor del momento angular del pseudopotencial
   INTEGER :: lxi,lyi,lzi,lambda !Auxiliares
   DOUBLE PRECISION :: Kmod,distcoefx, distcoefy,distcoefz, integral
   DOUBLE PRECISION :: acum
   DOUBLE PRECISION,DIMENSION (3) :: Kvector

   Z=ZlistECP(k)
   L=Lmax(Z)
   AAB_LOCAL_loops=0.d0
   acum=0.d0
   integral=0.d0
   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod= 2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)
!angular integrals should be calculated here

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
            acum= acum + integral*distcoefx * distcoefy * distcoefz *comb(lx,lxi) *comb(ly,lyi) * comb(lz,lzi)
            integral=0.d0
         END DO
      END DO
   END DO
   AAB_LOCAL_loops=acum
END FUNCTION AAB_LOCAL_loops




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION dAAB_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)
! calcula el termino semi-local del pseudopotencial centrado en i

!  l                                   
!  Σ[<xi|lm> Vi(l-LM) <lm|xj>] 
! m=-l                               

! los coef de la base se multiplican en la rutina que llama a esta

   USE basis_data, ONLY : a !a(i,ni) exponente de la funcion de base i, contrccion ni
   USE ECP_mod, ONLY :ZlistECP,Lmax,aECP,nECP,bECP, expnumbersECP,Qnl,Fulltimer_ECP,tsemilocal,tQ1
   use subm_intECP   , only: comb, OMEGA2, Aintegral, Qtype1N, Anal_radial_int

! ZlistECP(k) carga del atomo k con ECP
! Lmax(Z) maximo momento angular del pseudopotencial para el atomo con carga nuclear Z
! Vl= Σ aECP * r^nECP * exp(-bECP r^2)
! expnumbersECP(Z,l) terminos del ECP para el atomo con carga nuclear Z y momento angular l del ECP
!              ͚ 
! Qnl(n,l) = ʃ Ml(k*r)*r^n * exp(-cr^2) dr
!            ̊ 
! Fulltimer_ECP activa los timers para int. radiales
! tsemilocal,tQ1 auxiliares para el calculo de tiempos

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj
! i,j funciones de la base
! ii,ji numero de contraccion de la funcion
! k atomo con ecp
! lx, ly, lz; i,j exponente de la parte angular de la base x^lx y^ly z^lz
   DOUBLE PRECISION, DIMENSION(4) :: dAAB_SEMILOCAL !(<A|A|B>,<A|A|dB/dxB>,<A|A|dB/dyB>,<A|A|dB/dzB>
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz ! dx, dy, dz distancia entre nucleos
   DOUBLE PRECISION, DIMENSION(3) :: Kvector

   INTEGER :: l,m, term, lx,ly,lz, lambda,lmaxbase !auxiliades para ciclos
   INTEGER :: Z,n !Z= carga del nucleo
   DOUBLE PRECISION :: acumang, acumint, AABx, AABy, AABz, Kmod,Ccoef, auxdistx,auxdisty,auxdistz
!auxiliares
   INTEGER :: lambmin !minimo valor de lambda para integral angular no nula

   DOUBLE PRECISION :: t1,t2,t1q, t2q !auxiliares para timers

   INTEGER :: limitx, limity, limitz

   IF (Fulltimer_ECP) CALL cpu_time ( t1 )
   dAAB_SEMILOCAL=0.d0
   Z=ZlistECP(k)
   n=lxi+lxj+lyi+lyj+lzi+lzj
   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod=2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)
   lmaxbase=lxj+lyj+lzj
   AABx=0.d0
   AABy=0.d0
   AABz=0.d0
   acumint=0.d0
   acumang=0.d0

!aca deberian guardarse las integrales angulares

   DO l = 0 , Lmax(z)-1 !barre todos los l de la parte no local
      DO term=1, expnumbersECP(z,l) !barre contracciones del ECP para el atomo con carga z y l del ecp
         Ccoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
         Qnl=0.d0
         IF (Fulltimer_ECP) CALL cpu_time ( t1q )
         CALL Qtype1N(Kmod,Ccoef,lmaxbase+l+1,necp(Z,l,term)+n+1,necp(Z,l,term)) !calcula integrales radiales
!              ͚ 
! Qnl(n,l) = ʃ Ml(k*r)*r^n * exp(-cr^2) dr
!  

         CALL Anal_radial_int(1)

         IF (Fulltimer_ECP) THEN
            CALL cpu_time ( t2q )
            tQ1=tQ1 +t2q-t1q
         END IF

!FOCK
         dAAB_SEMILOCAL(1)=AAB_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz,l,term)
!d/dx 
         dAAB_SEMILOCAL(2)=dAAB_SEMILOCAL(2)+&
         AAB_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj+1,lyj,lzj,dx,dy,dz,l,term)*2.d0*a(j,ji)
         IF(lxj.gt.0) dAAB_SEMILOCAL(2)=dAAB_SEMILOCAL(2) &
         -AAB_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj-1,lyj,lzj,dx,dy,dz,l,term)*dble(lxj)
!d/dy 
         dAAB_SEMILOCAL(3)=dAAB_SEMILOCAL(3)+&
         AAB_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj+1,lzj,dx,dy,dz,l,term)*2.d0*a(j,ji)
         IF(lyj.gt.0) dAAB_SEMILOCAL(3)=dAAB_SEMILOCAL(3) &
         -AAB_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj-1,lzj,dx,dy,dz,l,term)*dble(lyj)
!d/dz 
         dAAB_SEMILOCAL(4)=dAAB_SEMILOCAL(4)& 
         +AAB_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj+1,dx,dy,dz,l,term)*2.d0*a(j,ji)
         IF(lzj.gt.0) dAAB_SEMILOCAL(4)=dAAB_SEMILOCAL(4)&
         -AAB_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj-1,dx,dy,dz,l,term)*dble(lzj)
      END DO
   END DO


   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tsemilocal=tsemilocal+t2-t1
   END IF
   RETURN
END FUNCTION dAAB_SEMILOCAL




DOUBLE PRECISION FUNCTION AAB_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz,l,term)
   USE basis_data, ONLY : a !a(i,ni) exponente de la funcion de base i, contrccion ni
   USE ECP_mod, ONLY :ZlistECP,Lmax,aECP,nECP,bECP, expnumbersECP,Qnl,Fulltimer_ECP,tsemilocal,tQ1
   use subm_intECP   , only: comb, OMEGA2, Aintegral, Qtype1N, Anal_radial_int

! ZlistECP(k) carga del atomo k con ECP
! Lmax(Z) maximo momento angular del pseudopotencial para el atomo con carga nuclear Z
! Vl= Σ aECP * r^nECP * exp(-bECP r^2)
! expnumbersECP(Z,l) terminos del ECP para el atomo con carga nuclear Z y momento angular l del ECP
!              ͚ 
! Qnl(n,l) = ʃ Ml(k*r)*r^n * exp(-cr^2) dr
!            ̊ 
! Fulltimer_ECP activa los timers para int. radiales
! tsemilocal,tQ1 auxiliares para el calculo de tiempos

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,l,term
! i,j funciones de la base
! ii,ji numero de contraccion de la funcion
! k atomo con ecp
! lx, ly, lz; i,j exponente de la parte angular de la base x^lx y^ly z^lz
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz ! dx, dy, dz distancia entre nucleos
   DOUBLE PRECISION, DIMENSION(3) :: Kvector

   INTEGER :: m, lx,ly,lz, lambda,lmaxbase !auxiliades para ciclos
   INTEGER :: Z,n !Z= carga del nucleo
   DOUBLE PRECISION :: acumang, acumint, AABx, AABy, AABz, Kmod,Ccoef, auxdistx,auxdisty,auxdistz
!auxiliares
   INTEGER :: lambmin !minimo valor de lambda para integral angular no nula

   DOUBLE PRECISION :: t1,t2,t1q, t2q !auxiliares para timers
   AAB_SEMILOCAL_loops=0.d0
   Z=ZlistECP(k)

   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod=2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)
   AABx=0.d0
   AABy=0.d0
   AABz=0.d0
   acumint=0.d0
   acumang=0.d0


   DO lx=0,lxj !barre potencias por expansion del binomio de Newton (x - dx)^lxj
      auxdistx=dx**(lxj-lx)
      IF (auxdistx .NE. 0.d0) THEN !si el factor de distancia es 0 deja de calcular
         DO ly=0,lyj
            auxdisty=dy**(lyj-ly)
            IF (auxdisty .NE. 0.d0) THEN
               DO lz=0,lzj
                  auxdistz=dz**(lzj-lz)
                  IF (auxdistz .NE. 0.d0) THEN
                     acumint=0.d0
                     lambmin=0

   IF (l-lxj-lyj-lzj .GT. 0) lambmin=l-lxj-lyj-lzj !minimo valor de lambda para integral angular no nula
   DO lambda=lxj+lyj+lzj+l,lambmin,-1
      acumang=0.d0
      DO m=-l,l
         acumang=acumang+Aintegral(l,m,lxi,lyi,lzi)*OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
      END DO
      acumint=acumint+acumang*Qnl(necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda)*aECP(z,L,term)
   END DO
   AABz=AABz+acumint * auxdistz* comb(lzj,lz)
   acumint=0.d0

                  END IF
               END DO
               AABy=AABy+AABz * auxdisty * comb(lyj,ly)
               AABz=0.d0
            END IF
         END DO
         AABx=AABx+AABy * auxdistx * comb(lxj,lx)
         AABy=0.d0
      END IF
   END DO
   AAB_SEMILOCAL_loops=AABx
END FUNCTION AAB_SEMILOCAL_loops



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!



FUNCTION dABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY :expnumbersECP, Qnl,bECP,IzECP,Fulltimer_ECP,tlocal,tQ1,Lmax,necp
   use subm_intECP   , only: OMEGA1, Q0, comb, Qtype1N, Anal_radial_int

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji !terminos de la base
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dx1,dy1,dz1,dx2,dy2,dz2 !distancias del centro con ecp a cada nucleo
   INTEGER, INTENT(IN) :: k !numero de atomo
   INTEGER :: L, Z !L maximo del ecp, Z carga nuclear sin modificar
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   DOUBLE PRECISION :: Kmod,Ccoef, integral,acum
   INTEGER :: lmaxbase,w
   DOUBLE PRECISION :: t1,t2,t1q,t2q ! auxiliares para timers
   DOUBLE PRECISION, DIMENSION(7) :: dABC_LOCAL !(<A|B|C>,<dA/dxA|B|C>,<dA/dyA|B|C>,<dA/dzA|B|C>,<A|B|dC/dxC>,<A|B|dC/dyC>,<A|B|dC/dzC>

   IF (Fulltimer_ECP) CALL cpu_time ( t1 )


   dABC_LOCAL=0.d0
   z=IzECP(k)
   L=Lmax(Z)
   lmaxbase=lxi+lyi+lzi+lxj+lyj+lzj
   Kvector=(/a(i,ii)*dx1+a(j,ji)*dx2,a(i,ii)*dy1+a(j,ji)*dy2,a(i,ii)*dz1+a(j,ji)*dz2/)
   Kvector=-2.d0*Kvector
   Kmod=sqrt(Kvector(1)**2.d0+Kvector(2)**2.d0+Kvector(3)**2.d0)
   integral=0.d0
   acum=0.d0
!aca deberia meter en memoria las integralas angulares

   DO w =1, expnumbersECP(z,l) !barre terminos del ECP para el atomo con carga nuclear Z y l del ECP
      Qnl=0.d0
      Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
      IF (Fulltimer_ECP) CALL cpu_time ( t1q )

      CALL Qtype1N(Kmod,Ccoef,Lmaxbase+1,necp(Z,l,w)+Lmaxbase+3,necp(Z,l,w)) !calcula integral radial
      CALL Anal_radial_int(1)

      IF (Fulltimer_ECP) THEN
         CALL cpu_time ( t2q )
         tQ1=tQ1 +t2q-t1q
      END IF


      dABC_LOCAL=0.d0
!%%%%%%%%%%%%
!FOCK
      dABC_LOCAL(1)=dABC_LOCAL(1)+ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)
!d/dxi 
      dABC_LOCAL(2)=dABC_LOCAL(2)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi+1,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(i,ii)
      if(lxi.ge.1) dABC_LOCAL(2)=dABC_LOCAL(2)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi-1,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lxi)
!d/dyi
      dABC_LOCAL(3)=dABC_LOCAL(3)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi+1,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(i,ii)
      if(lyi.ge.1) dABC_LOCAL(3)=dABC_LOCAL(3)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi-1,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lyi)
!d/dzi
      dABC_LOCAL(4)=dABC_LOCAL(4)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi+1,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(i,ii)
      if(lzi.ge.1) dABC_LOCAL(4)=dABC_LOCAL(4)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi-1,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lzi)
!d/dxj
      dABC_LOCAL(5)=dABC_LOCAL(5)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj+1,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(j,ji)
      if(lxj.ge.1) dABC_LOCAL(5)=dABC_LOCAL(5)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj-1,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lxj)
!d/dyj
      dABC_LOCAL(6)=dABC_LOCAL(6)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj+1,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(j,ji)
      if(lyj.ge.1) dABC_LOCAL(6)=dABC_LOCAL(6)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj-1,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lyj)
!d/dzj
      dABC_LOCAL(7)=dABC_LOCAL(7)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj+1,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(j,ji)
      if(lzj.ge.1) dABC_LOCAL(7)=dABC_LOCAL(7)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj-1,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lzj)
   END DO

   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tlocal=tlocal+t2-t1
   END IF
   RETURN
END FUNCTION dABC_LOCAL


DOUBLE PRECISION FUNCTION ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : Qnl,IzECP,angularint,pi,Fulltimer_ECP,Lmax,necp,aECP
   use subm_intECP   , only: OMEGA1, Q0, comb, Qtype1N

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji,w !terminos de la base
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dx1,dy1,dz1,dx2,dy2,dz2 !distancias del centro con ecp a cada nucleo
   INTEGER, INTENT(IN) :: k !numero de atomo
   INTEGER :: L, Z !L maximo del ecp, Z carga nuclear sin modificar
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   DOUBLE PRECISION :: Kmod,Ccoef, integral,auxcomb,auxdist,acum, auxdista, auxdistb, auxdistc, auxdistd, auxdiste, auxdistf
   INTEGER :: lmaxbase
   INTEGER :: ac,bc,cc,dc,ec,fc,lambda ! variables auxiliares
   DOUBLE PRECISION :: t1 ! auxiliares para timers
   IF (Fulltimer_ECP) CALL cpu_time ( t1 )


   ABC_LOCAL_loops=0.d0
   z=IzECP(k)
   L=Lmax(Z)
   lmaxbase=lxi+lyi+lzi+lxj+lyj+lzj
   Kvector=(/a(i,ii)*dx1+a(j,ji)*dx2,a(i,ii)*dy1+a(j,ji)*dy2,a(i,ii)*dz1+a(j,ji)*dz2/)
   Kvector=-2.d0*Kvector
   Kmod=sqrt(Kvector(1)**2.d0+Kvector(2)**2.d0+Kvector(3)**2.d0)
   integral=0.d0
   acum=0.d0

   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
      auxdista=dx1**(lxi-ac)
      IF (auxdista .NE. 0.d0) THEN !si el factor de distancia es 0 deja de calcular
         DO bc=0,lyi
            auxdistb=dy1**(lyi-bc)
            IF (auxdistb .NE. 0.d0) THEN
               DO cc=0,lzi
                  auxdistc=dz1**(lzi-cc)
                  IF (auxdistc .NE. 0.d0) THEN
                     DO dc=0,lxj
                        auxdistd=dx2**(lxj-dc)
                        IF (auxdistd .NE. 0.d0) THEN
                           DO ec=0,lyj
                              auxdiste=dy2**(lyj-ec)
                              IF (auxdiste .NE. 0.d0) THEN
                                 DO fc=0,lzj
                                    auxdistf=dz2**(lzj-fc)
                                    IF (auxdistf .NE. 0.d0) THEN

   auxdist=auxdista*auxdistb*auxdistc*auxdistd*auxdiste*auxdistf
   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
      IF ( Kmod .GT. 0.d0 ) THEN
         integral=integral + OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc) * Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)
         IF (Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda) .NE. Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)) &
         STOP " qnl1l2 = 0 in ABC_SEMILOCAL" !corta el calculo si integral radial es 0
      ELSE
         integral=integral + angularint(ac+dc,bc+ec,cc+fc) * Q0(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),Ccoef) *0.25d0/pi
!parche para el caso accidental en que K fuera = (0,0,0) por compensacion de a(i,ii)*dx1+a(j,ji)*dx2 (idem y,z)
!0.25d0/pi compensa el factor 4pi por el que se multiplica luego a la suma de las integrales
      END IF
   END DO
   auxcomb=comb(lxi,ac)*comb(lyi,bc)*comb(lzi,cc)*comb(lxj,dc)*comb(lyj,ec)*comb(lzj,fc)
   acum=acum + auxcomb*auxdist*integral

   integral=0.d0
   auxcomb=0.d0
   auxdist=0.d0

                                    END IF
                                 END DO
                              END IF
                           END DO
                        END IF
                     END DO
                  END IF
               END DO
            END IF
         END DO
      END IF
   END DO
   ABC_LOCAL_loops=aECP(z,L,w)*acum
   RETURN
END FUNCTION ABC_LOCAL_loops




FUNCTION dABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : Qnl1l2,necp,bECP,IzECP,Fulltimer_ECP,tsemilocal,tQ2,Lmax,expnumbersECP
   use subm_intECP   , only: Qtype2N, Anal_radial_int
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji,k
!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
!k atomo con ecp
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dxi,dyi,dzi,dxj,dyj,dzj !distancias de las nucleos con bases al nucleo con ecp
   INTEGER :: Z,l1max,l2max !Z= carga del nucleo, y momento angular de la base i y j
   DOUBLE PRECISION,DIMENSION (3) :: Kivector,Kjvector
   DOUBLE PRECISION :: Kimod, Kjmod,Ccoef
   INTEGER :: l,term !auxiliares ciclos
   DOUBLE PRECISION :: acumang,acumang1,acumang2,integral,auxcomb,auxdist,acum
   DOUBLE PRECISION, DIMENSION(7) :: dABC_SEMILOCAL !(<A|B|C>,<dA/dxA|B|C>,<dA/dyA|B|C>,<dA/dzA|B|C>,<A|B|dC/dxC>,<A|B|dC/dyC>,<A|B|dC/dzC>
!auxiliares
   DOUBLE PRECISION :: t1,t2, t1q,t2q,t1aux,t2aux !auxiliares para timers
   IF (Fulltimer_ECP) CALL cpu_time ( t1 )
   t1aux=0.d0
   t2aux=0.d0

   dABC_SEMILOCAL=0.d0
   integral=0.d0
   acum=0.d0
   Qnl1l2=0.d0
   acumang=0.d0
   acumang1=0.d0
   acumang2=0.d0
   auxcomb=0.d0
   auxdist=0.d0

   Z=IzECP(k)

   l1max=lxi+lyi+lzi
   l2max=lxj+lyj+lzj

   Kivector=-2.d0*a(i,ii)*(/dxi,dyi,dzi/)
   Kjvector=-2.d0*a(j,ji)*(/dxj,dyj,dzj/)

   Kimod=sqrt(Kivector(1)**2.d0+Kivector(2)**2.d0+Kivector(3)**2.d0)
   Kjmod=sqrt(Kjvector(1)**2.d0+Kjvector(2)**2.d0+Kjvector(3)**2.d0)

!aca deberia calcular integrales angulares

   DO l = 0 , Lmax(z)-1 !barre todos los l de la parte no local
      DO term=1, expnumbersECP(Z,l) !barre contracciones del ECP para el atomo con carga Z y momento angular l del ecp
         Qnl1l2=0.d0
         Ccoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
         IF (Fulltimer_ECP) CALL cpu_time ( t1q )
         call Qtype2N(Kimod,Kjmod,Ccoef,l1max+l+1,l2max+l+1,necp(Z,l,term)+l1max+l2max+4,necp(Z,l,term)) !agrega a la matriz Qnl1l2 los terminos correspondientes a un termino radiales.
         call Anal_radial_int(2)

         IF (Fulltimer_ECP) THEN
            CALL cpu_time ( t2q )
            tQ2=tQ2 +t2q-t1q
         END IF

!Fock
   dABC_SEMILOCAL(1)=dABC_SEMILOCAL(1)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)
!d/dxi
   dABC_SEMILOCAL(2)=dABC_SEMILOCAL(2)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi+1,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(i,ii)
   if(lxi.ge.1) dABC_SEMILOCAL(2)=dABC_SEMILOCAL(2)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi-1,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lxi)
!d/dyi
   dABC_SEMILOCAL(3)=dABC_SEMILOCAL(3)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi+1,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(i,ii)
   if(lyi.ge.1) dABC_SEMILOCAL(3)=dABC_SEMILOCAL(3)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi-1,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lyi)
!d/dzi
   dABC_SEMILOCAL(4)=dABC_SEMILOCAL(4)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi+1,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(i,ii)
   if(lzi.ge.1) dABC_SEMILOCAL(4)=dABC_SEMILOCAL(4)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi-1,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lzi)
!d/dxj
   dABC_SEMILOCAL(5)=dABC_SEMILOCAL(5)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj+1,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(j,ji)
   if(lxj.ge.1) dABC_SEMILOCAL(5)=dABC_SEMILOCAL(5)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj-1,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lxj)
!d/dyj
   dABC_SEMILOCAL(6)=dABC_SEMILOCAL(6)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj+1,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(j,ji)
   if(lyj.ge.1) dABC_SEMILOCAL(6)=dABC_SEMILOCAL(6)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj-1,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lyj)
!d/dzj
   dABC_SEMILOCAL(7)=dABC_SEMILOCAL(7)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj+1,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(j,ji)
   if(lzj.ge.1) dABC_SEMILOCAL(7)=dABC_SEMILOCAL(7)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj-1,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lzj)
      END DO
   END DO

   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tsemilocal=tsemilocal+t2-t1
   END IF

   RETURN
END FUNCTION dABC_SEMILOCAL





DOUBLE PRECISION FUNCTION ABC_SEMILOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : Qnl1l2,necp,Fulltimer_ECP,tsemilocal,Taux,aECP
   use subm_intECP   , only: comb,OMEGA2
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji,k
!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
!k atomo con ecp
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dxi,dyi,dzi,dxj,dyj,dzj !distancias de las nucleos con bases al nucleo con ecp
   INTEGER, INTENT(IN) :: l,term,Z
   INTEGER :: l1max,l2max !Z= carga del nucleo, y momento angular de la base i y j
   DOUBLE PRECISION,DIMENSION (3) :: Kivector,Kjvector
   DOUBLE PRECISION :: Kimod, Kjmod
   INTEGER :: ac,bc,cc,dc,ec,fc,lambdai, lambdaj,m,lambimin, lambjmin !auxiliares ciclos
   DOUBLE PRECISION :: acumang,acumang1,acumang2,integral,auxcomb,auxdist,acum,auxdista,auxdistb,auxdistc,auxdistd,auxdiste,auxdistf!auxiliares
   DOUBLE PRECISION :: t1,t2, t1aux,t2aux !auxiliares para timers
   IF (Fulltimer_ECP) CALL cpu_time ( t1 )
   t1aux=0.d0
   t2aux=0.d0

   ABC_SEMILOCAL_loops=0.d0
   integral=0.d0
   acum=0.d0
   acumang=0.d0
   acumang1=0.d0
   acumang2=0.d0
   auxcomb=0.d0
   auxdist=0.d0

   l1max=lxi+lyi+lzi
   l2max=lxj+lyj+lzj

   Kivector=-2.d0*a(i,ii)*(/dxi,dyi,dzi/)
   Kjvector=-2.d0*a(j,ji)*(/dxj,dyj,dzj/)

   Kimod=sqrt(Kivector(1)**2.d0+Kivector(2)**2.d0+Kivector(3)**2.d0)
   Kjmod=sqrt(Kjvector(1)**2.d0+Kjvector(2)**2.d0+Kjvector(3)**2.d0)

   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dxi)^lxi
      auxdista=dxi**(lxi-ac)
      IF (auxdista .NE. 0.d0) THEN !si el factor de distancia es 0 deja de calcular
         DO bc=0,lyi
            auxdistb=dyi**(lyi-bc)
            IF (auxdistb .NE. 0.d0) THEN
               DO cc=0,lzi
                  auxdistc=dzi**(lzi-cc)
                  IF (auxdistc .NE. 0.d0) THEN
                     DO dc=0,lxj
                        auxdistd=dxj**(lxj-dc)
                        IF (auxdistd .NE. 0.d0) THEN
                           DO ec=0,lyj
                              auxdiste=dyj**(lyj-ec)
                              IF (auxdiste .NE. 0.d0) THEN
                                 DO fc=0,lzj
                                    auxdistf=dzj**(lzj-fc)
                                    IF (auxdistf .NE. 0.d0) THEN


   IF (Fulltimer_ECP) CALL cpu_time ( t1aux )

   auxdist=auxdista*auxdistb*auxdistc*auxdistd*auxdiste*auxdistf
   auxcomb=comb(lxi,ac)*comb(lyi,bc)*comb(lzi,cc)*comb(lxj,dc)*comb(lyj,ec)*comb(lzj,fc)
   lambimin=0
   lambjmin=0

   IF (l-ac-bc-cc .GT. 0) lambimin=l-ac-bc-cc !lambda minimo que no anula la integral angular
   IF (l-dc-ec-fc .GT. 0) lambjmin=l-dc-ec-fc !lambda minimo que no anula la integral angular

   DO lambdai=ac+bc+cc+l,lambimin,-2
      DO lambdaj=dc+ec+fc+l,lambjmin,-2
         acumang=0.d0
         DO m=-l,l
            acumang=acumang+OMEGA2(Kivector,lambdai,l,m,ac,bc,cc)*OMEGA2(Kjvector,lambdaj,l,m,dc,ec,fc)
         END DO
         integral=integral+acumang*Qnl1l2(ac+bc+cc+dc+ec+fc+necp(Z,l,term),lambdai,lambdaj)
!         IF (Qnl1l2(ac+bc+cc+dc+ec+fc+necp(Z,l,term),lambdai,lambdaj) .EQ. 0.d0)  STOP " q = 0 in abc semiloc2"
         acumang=0.d0
      END DO
   END DO
   acum=acum + auxcomb*auxdist*integral
   integral=0.d0
   auxcomb=0.d0
   auxdist=0.d0

   IF (Fulltimer_ECP) CALL cpu_time ( t2aux )
   Taux=Taux+t2aux-t1aux

                                    END IF
                                 END DO
                              END IF
                           END DO
                        END IF
                     END DO
                  END IF
               END DO
            END IF
         END DO
      END IF
   END DO

   ABC_SEMILOCAL_loops=acum*aECP(Z,L,term)

   IF (Fulltimer_ECP) THEN
      CALL cpu_time ( t2 )
      tsemilocal=tsemilocal+t2-t1
   END IF

   RETURN
END FUNCTION ABC_SEMILOCAL_loops

end module subm_intECPG

