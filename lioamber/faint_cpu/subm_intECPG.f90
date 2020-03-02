!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INTECPG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This routines calculate Fock matrix elements for effective core              !
! potential (ECP) and it derivatives at same time.                             !
!                                                                              !
! EXTERNAL INPUT:                                                              !
!                                                                              !
! INTERNAL INPUT: system, basis set and ECP information.                       !
!   · system information:                                                      !
!        · natom: number of QM atoms.                                          !
!   · basis set information:                                                   !
!        · ncont(i): number of contractions of function i-th.                  !
!        · a(i,ii): basis function exponents.                                  !
!        · Cnorm(i,ii): basis function coefficients, renormalized.             !
!         · M: total number of basis functions                                 !
!        · Nuc(i): atomic index corresponding to function i.                   !
!        · Lxyz(i,1:3) exponents of the angular part for i function            !
!             · 1 is (x/r), 2 is (y/r) and 3 is (z/r)                          !
!   · ECP information:                                                         !
!        · ecptypes: number of atoms with an ECP                               !
!        · IzECP(i): original nuclear charge of the i-th atom                  ! 
!        · ZlistECP(k): nuclear charge of k-th atom defined in input with ECP  !
!        · Lmax(Z): max angular moment of ECP for the atom with nuclear        !
!          charge Z                                                            !
!        · expnumbersECP(Z,l): contraciones of ECP for nuclear charge Z,       !
!          angular moment l                                                    !
!        · nECP(Z,l): ECP polinomial part exponent                             !
!        · bECP(Z,l): ECP Gaussian part exponents.                             !
!        · aECP(Z,l): ECP coefficients.                                        !
!                                                                              !
!   · ECP calculation variables:                                               !
!        · cut2_0: cutoff value for 2 center integrals                         !
!        · cut3_0: cutoff value for 3 center integrals                         !
!                                                                              !
! OUTPUTS:                                                                     !
!   · VAAB: 2 center terms integrals of ECP (<A|A|B>, <B|A|A>)                 !
!   · VBAC: 3 center terms integrals of ECP (<B|A|C> & <B|A|B>)                !
!   · dHcore_AAB 2 center gradients                                            !
!   · dHcore_ABC 3 center gradients                                            !
!                                                                              !
! V 2.00 January 2020 Dbug and optimized                                       !
! V 1.00 November 2019 First version                                           !
!                                                                              !
! Nicolas Foglia                                                               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! References:                                                                  !
! J. Chem. Phys. 65, 3826 (1976); http://dx.doi.org/10.1063/1.432900           !
! J. Chem. Phys. 111, 8778 (1999); http://dx.doi.org/10.1063/1.480225          !
! Foglia N. O. <<Simulacion computacional de dinamica electronica y            !
! reactividad quımica>>. PhD Tesis. University of Buenos Aires, 2019           ! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

#include "../datatypes/datatypes.fh"

module subm_intECPG
contains
subroutine intECPG()
   use basis_data   , only: Nuc,M, ncont, a
   use garcha_mod, only : natom
   use ECP_mod, ONLY : Lxyz, ecptypes, IzECP, Cnorm, pi, ZlistECP, distx, disty, distz,ECPatoms_order, &
#ifdef FULL_CHECKS
   ECPatoms, &
#endif
   dHcore_AAB, dHcore_ABC,VAAB,VBAC, cut3_0, cut2_0 
   use subm_intECP   , only: AAB_LOCAL, AAB_SEMILOCAL, ABC_LOCAL, ABC_SEMILOCAL
   implicit none
   integer :: i, j, k !number of basis set function
   integer :: ii, ji !number of contraction
   integer :: kecp !atoms with ECP
   integer :: lxi,lxj,lyi,lyj,lzi,lzj !l?$  potencia de la parte angular de la base
   LIODBLE :: ABC, dABCpl, dABCpr
   LIODBLE :: acum
   LIODBLE :: exp_cut
   LIODBLE :: acuml, acumr, Distcoef, dxi, dyi, dzi, dxj, dyj, dzj, dx, dy, dz
   LIODBLE, dimension(4) :: dHcore_AAB_temp 
   LIODBLE, dimension(7) :: dHcore_ABC_temp, dHcore_ABC_temp_aux
   integer :: pos

   LIODBLE :: exp_Distcoef
#ifdef FULL_CHECKS
   integer :: l
#endif


   call g2g_timer_start('ECP_full')
   call g2g_timer_sum_start('ECP_full')

   VAAB=0.d0
   VBAC=0.d0
   dHcore_AAB=0.d0
   dHcore_ABC=0.d0
   dHcore_ABC_temp=0.d0


! Computting 2 center contributions !<Xa|Va|Xb>,<Xb|Va|Xa>. <Xb|Va|Xb> NOT included here.
   call g2g_timer_sum_start('ECP_2Centers')
   exp_cut = exp(-cut2_0)

   do i = 1, M
   do j = 1, i
      pos=i+(1-j)*(j-2*M)/2 !posicion en el vector que guarda la matriz de FOCK triangular
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
 
              IF (IzECP(nuc(i)) .EQ. ZlistECP(kecp)) THEN !calculo para ECP en i d<A|A|B>

                  DO ji=1, ncont(j) !barre contracciones de la base j
                     exp_Distcoef=exp(-Distcoef*a(j,ji))
                     if (exp_Distcoef.gt.exp_cut) then !cutoff for integrals
                        dHcore_AAB_temp=0.d0


                        DO ii=1, ncont(i) !ii barre contracciones de las funcion de base i
   dHcore_AAB_temp(1:4)=dHcore_AAB_temp(1:4)+dAAB_LOCAL(i,j,kecp,ii,ji,lxj,lyj,lzj,lxi,lyi,lzi,dx,dy,dz)*Cnorm(i,ii)
   dHcore_AAB_temp(1:4)=dHcore_AAB_temp(1:4)+dAAB_SEMILOCAL(i,j,ii,ji,kecp,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)*Cnorm(i,ii)
                        END DO

   VAAB(pos) = VAAB(pos)+ dHcore_AAB_temp(1)*Cnorm(j,ji)*4.d0*pi*exp_Distcoef

   dHcore_AAB(pos,1,1:3)=dHcore_AAB(pos,1,1:3)-dHcore_AAB_temp(2:4)*Cnorm(j,ji)*4.d0*pi*exp_Distcoef/0.529177D0 !di/dx,y,z
   dHcore_AAB(pos,2,1:3)=dHcore_AAB(pos,2,1:3)+dHcore_AAB_temp(2:4)*Cnorm(j,ji)*4.d0*pi*exp_Distcoef/0.529177D0 !dj/dx,y,z

   dHcore_AAB_temp=0.d0


                     end if
                  END DO
               END IF


               IF (IzECP(nuc(j)) .EQ. ZlistECP(kecp)) THEN !calculo para ECP en j d<B|A|A>
                  DO ii=1, ncont(i) ! barre contracciones de las funcion de base i
                     exp_Distcoef=exp(-Distcoef*a(i,ii))
                     if (exp_Distcoef.gt.exp_cut) then !cutoff for integrals
                        dHcore_AAB_temp=0.d0


                        DO ji=1, ncont(j) !barre contracciones de las funcion de base j
   dHcore_AAB_temp(1:4)=dHcore_AAB_temp(1:4)+ dAAB_LOCAL(j,i,kecp,ji,ii,lxi,lyi,lzi,lxj,lyj,lzj,-dx,-dy,-dz)*Cnorm(j,ji)
   dHcore_AAB_temp(1:4)=dHcore_AAB_temp(1:4)+ dAAB_SEMILOCAL(j,i,ji,ii,kecp,lxj,lyj,lzj,lxi,lyi,lzi,-dx,-dy,-dz)*Cnorm(j,ji)
                        END DO

   VAAB(pos) = VAAB(pos)+ dHcore_AAB_temp(1)*Cnorm(i,ii)*4.d0*pi*exp_Distcoef

   dHcore_AAB(pos,1,1:3)=dHcore_AAB(pos,1,1:3)+dHcore_AAB_temp(2:4)*Cnorm(i,ii)*4.d0*pi*exp_Distcoef/0.529177D0
   dHcore_AAB(pos,2,1:3)=dHcore_AAB(pos,2,1:3)-dHcore_AAB_temp(2:4)*Cnorm(i,ii)*4.d0*pi*exp_Distcoef/0.529177D0

   dHcore_AAB_temp=0.d0


                     end if
                  END DO
               END IF

            END IF
         END DO

   end if
   end do
   end do
   call g2g_timer_sum_pause('ECP_2Centers')


! Computing 3 center terms <Xb|Va|Xc> and 2 center terms <Xb|Va|Xb>
   call g2g_timer_sum_start('ECP_3Centers')
    exp_cut = exp(-cut3_0)
   DO i=1,M !barre funciones de la base
   DO j=1,i !barre funciones de la base
      pos=i+(1-j)*(j-2*M)/2 !posicion en el vector que guarda la matriz de FOCK triangular
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
                        exp_Distcoef=exp(-Distcoef)


   IF (exp_Distcoef.gt.exp_cut) then

   dHcore_ABC_temp_aux=dABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)*Cnorm(j,ji)*exp_Distcoef
   dHcore_ABC_temp=dHcore_ABC_temp+dHcore_ABC_temp_aux

   dHcore_ABC_temp_aux=4.d0*pi*dABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)*Cnorm(j,ji)*exp_Distcoef
   dHcore_ABC_temp=dHcore_ABC_temp+dHcore_ABC_temp_aux

   end if


                     END DO


   VBAC(pos) = VBAC(pos) + dHcore_ABC_temp(1)*Cnorm(i,ii)*4.d0*pi

   dHcore_ABC(pos,1,1:3)=dHcore_ABC(pos,1,1:3)+dHcore_ABC_temp(2:4)*Cnorm(i,ii)*4.d0*pi/0.529177D0
   dHcore_ABC(pos,2,1:3)=dHcore_ABC(pos,2,1:3)+dHcore_ABC_temp(5:7)*Cnorm(i,ii)*4.d0*pi/0.529177D0
   dHcore_ABC(pos,2+ECPatoms_order(k),1:3)=dHcore_ABC(pos,2+ECPatoms_order(k),1:3)- &
   dHcore_ABC_temp(2:4)*Cnorm(i,ii)*4.d0*pi/0.529177D0-dHcore_ABC_temp(5:7)*Cnorm(i,ii)*4.d0*pi/0.529177D0
   dHcore_ABC_temp=0.d0

                  END DO
               END IF
            END DO
         END IF
      END DO
   end do
   end do
   call g2g_timer_sum_pause('ECP_3Centers')


#ifdef FULL_CHECKS
!NAN in 2 center terms
   do i=1,M
      do j=1,M
         do l=1,3
             if (dHcore_AAB(pos,1,l).ne.dHcore_AAB(pos,1,l)) then
                write(*,*) "NAN en: dHcore_AAB", i,j,"1",l
                stop
             elseif (dHcore_AAB(pos,2,l).ne.dHcore_AAB(pos,2,l)) then
                write(*,*) "NAN en: dHcore_AAB", i,j,"2",l
                stop
             end if
          end do
       end do
    end do

!NAN in 3 center terms
   do i=1,M
      do j=1,M
         do k=1,2+ECPatoms
            do l=1,3
               if (dHcore_ABC(pos,k,l).ne.dHcore_ABC(pos,k,l)) then
                write(*,*) "NAN en: dHcore_ABC", i,j,k,l
                stop
               end if
            end do
          end do
       end do
    end do
#endif

   call g2g_timer_stop('ECP_full')
   call g2g_timer_sum_stop('ECP_full')
   return;
END SUBROUTINE intECPG

SUBROUTINE ECP_gradients(ff,rho,natom)
   use basis_data   , only: M, Nuc
   use ECP_mod, ONLY : dHcore_AAB, dHcore_ABC, ecptypes, IzECP, ZlistECP, ECPatoms_order
   IMPLICIT NONE
   integer         , intent(in)  :: natom
   LIODBLE, intent(out) :: ff(natom,3)
   LIODBLE, intent(in)  :: rho(:) !why this havent got dimention ??? Nick
   integer :: i,j,l,k, kecp, pos
   LIODBLE, dimension(natom) :: F_i
   call g2g_timer_start('ECP_grad')
   call g2g_timer_sum_start('ECP_grad')
   ff=0.d0
   do i=1,M
      do j=1,i
         pos=i+(1-j)*(j-2*M)/2
         do l=1,3
            F_i=0.d0
            do k=1,natom
               if (k.eq.nuc(i)) F_i(k)=F_i(k)+(dHcore_AAB(pos,1,l)+dHcore_ABC(pos,1,l))*0.529177D0
               if (k.eq.nuc(j)) F_i(k)=F_i(k)+(dHcore_AAB(pos,2,l)+dHcore_ABC(pos,2,l))*0.529177D0
               do kecp=1, ecptypes !barre atomos con ecp
                  if (IzECP(k) .EQ. ZlistECP(kecp)) THEN !solo calcula si el nucleo tiene ecp
                     F_i(k)=F_i(k)+dHcore_ABC(pos,2+ECPatoms_order(k),l)*0.529177D0
                  end if
               end do
               ff(k,l)=ff(k,l)+F_i(k)*rho(pos)
            enddo
         end do
      end do
   enddo
   call g2g_timer_stop('ECP_grad')
   call g2g_timer_sum_stop('ECP_grad')
   return
END SUBROUTINE ECP_gradients



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    2 Center terms <Xa|Va|Xb>    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

FUNCTION dAAB_LOCAL(i,j,k,ii,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz)
!Calcula el termino local del pseudopotencial centrado en i [<xi|Vi(LM)|xj>] y sus derivadas 
!los coef de la base se multiplican en la rutina que llama a esta
!i,j funcion de base
!ii,ji numero de contraccion de la funcion de base
!k atomo con ECP
USE basis_data, ONLY : a
USE ECP_mod, ONLY :nECP,bECP, aECP, ZlistECP, Lmax, expnumbersECP, Qnl
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
   INTEGER :: w !Auxiliares
   DOUBLE PRECISION :: Ccoef, Kmod, integral
   DOUBLE PRECISION, DIMENSION (4) :: acum
   DOUBLE PRECISION,DIMENSION (3) :: Kvector

   call g2g_timer_sum_start('ECP_2C_local')

   Z=ZlistECP(k)
   L=Lmax(Z)
   Lmaxbase=lx+ly+lz+kxi+kyi+kzi+1
   dAAB_LOCAL=0.d0
   acum=0.d0
   integral=0.d0
   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod= 2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)

   CALL AAB_LOCAL_angular(j,k,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz) ! angular integrals calculated here and stack in ECP_Ang_stack
   DO w =1, expnumbersECP(z,l) !barre todos los terminos del Lmaximo
      Qnl=0.d0
      Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)
      CALL Qtype1N(Kmod,Ccoef,Lmaxbase,necp(Z,l,w)+Lmaxbase,necp(Z,l,w)+kxi+kyi+kzi) ! angular integrals calculated here and stack in Qnl

!Fock terms
      acum(1)= acum(1) + AAB_LOCAL_loops(j,k,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz,w)
!d/dx terms
      acum(2)= acum(2) + AAB_LOCAL_loops(j,k,ji,lx+1,ly,lz,kxi,kyi,kzi,dx,dy,dz,w)*2.d0*a(j,ji)
      if (lx .gt.0) acum(2)= acum(2) - AAB_LOCAL_loops(j,k,ji,lx-1,ly,lz,kxi,kyi,kzi,dx,dy,dz,w)*dble(lx)
!d/dy terms
      acum(3)= acum(3) + AAB_LOCAL_loops(j,k,ji,lx,ly+1,lz,kxi,kyi,kzi,dx,dy,dz,w)*2.d0*a(j,ji)
      if (ly .gt.0) acum(3)= acum(3) - AAB_LOCAL_loops(j,k,ji,lx,ly-1,lz,kxi,kyi,kzi,dx,dy,dz,w)*dble(ly)
!d/dz terms
      acum(4)= acum(4) + AAB_LOCAL_loops(j,k,ji,lx,ly,lz+1,kxi,kyi,kzi,dx,dy,dz,w)*2.d0*a(j,ji)
      if (lz .gt.0) acum(4)= acum(4) - AAB_LOCAL_loops(j,k,ji,lx,ly,lz-1,kxi,kyi,kzi,dx,dy,dz,w)*dble(lz)

      dAAB_LOCAL=dAAB_LOCAL+aECP(z,L,w)*acum 
      acum=0.d0
   END DO

   call g2g_timer_sum_pause('ECP_2C_local')
   RETURN
END FUNCTION dAAB_LOCAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine AAB_LOCAL_angular(j,k,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : ZlistECP,Lmax, ECP_Ang_stack
   use subm_intECP   , only: comb, OMEGA1
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: j,k,ji,lx,ly,lz,kxi,kyi,kzi
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz
   INTEGER :: z,l
! Z carga nuclear
! l maximo valor del momento angular del pseudopotencial
   INTEGER :: lxi,lyi,lzi,lambda !Auxiliares
   DOUBLE PRECISION :: Kmod, integral
   DOUBLE PRECISION :: acum
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   INTEGER :: pos1

   Z=ZlistECP(k)
   L=Lmax(Z)
   acum=0.d0
   integral=0.d0
   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod= 2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)

!base
   DO lxi=0,lx !barre potencias por expansion del binomio de Newton (x - dx)^lx
      DO lyi=0,ly
         DO lzi=0,lz
            pos1=lxi+4*lyi+16*lzi+1
            DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
               ECP_Ang_stack(0,lambda,0,pos1,1)=OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi)
            END DO
         END DO
      END DO
   END DO

!x+1
   lxi=lx+1 !barre potencias por expansion del binomio de Newton (x - dx)^lx
   DO lyi=0,ly
      DO lzi=0,lz
         pos1=lxi+4*lyi+16*lzi+1
         DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
            ECP_Ang_stack(0,lambda,0,pos1,1)=OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi) 
         END DO
      END DO
   END DO

!y+1
   lyi=ly+1
   DO lxi=0,lx !barre potencias por expansion del binomio de Newton (x - dx)^lx
      DO lzi=0,lz
         pos1=lxi+4*lyi+16*lzi+1
         DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
            ECP_Ang_stack(0,lambda,0,pos1,1)=OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi)                            
         END DO
      END DO
   END DO

!z+1
   lzi=lz+1
   DO lxi=0,lx !barre potencias por expansion del binomio de Newton (x - dx)^lx
      DO lyi=0,ly
         pos1=lxi+4*lyi+16*lzi+1
         DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2
            ECP_Ang_stack(0,lambda,0,pos1,1)=OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi)                            
         END DO
      END DO
   END DO

END SUBROUTINE AAB_LOCAL_angular

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

DOUBLE PRECISION FUNCTION AAB_LOCAL_loops(j,k,ji,lx,ly,lz,kxi,kyi,kzi,dx,dy,dz,w)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : Qnl,necp, ZlistECP,Lmax, ECP_Ang_stack
   use subm_intECP   , only: comb, OMEGA1
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: j,k,ji,lx,ly,lz,kxi,kyi,kzi,w
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz
   INTEGER :: z,l
! Z carga nuclear
! l maximo valor del momento angular del pseudopotencial
   INTEGER :: lxi,lyi,lzi,lambda !Auxiliares
   DOUBLE PRECISION :: Kmod,distcoefx, distcoefy,distcoefz, integral
   DOUBLE PRECISION :: acum
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   INTEGER :: pos1

   Z=ZlistECP(k)
   L=Lmax(Z)
   AAB_LOCAL_loops=0.d0
   acum=0.d0
   integral=0.d0
   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod= 2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)

   DO lxi=0,lx !barre potencias por expansion del binomio de Newton (x - dx)^lx
      distcoefx=dx**(lx-lxi)
      DO lyi=0,ly
         distcoefy=dy**(ly-lyi)
         DO lzi=0,lz
            distcoefz=dz**(lz-lzi)
            pos1=lxi+4*lyi+16*lzi+1
            DO lambda=lxi+lyi+lzi+kxi+kyi+kzi,0,-2

#ifdef FULL_CHECKS
   if ((ECP_Ang_stack(0,lambda,0,pos1,1)).ne.OMEGA1(Kvector,lambda,lxi+kxi,lyi+kyi,lzi+kzi)) STOP "bad angular int 2C_LOC"
#endif
 
               integral=integral + ECP_Ang_stack(0,lambda,0,pos1,1) * &
               Qnl(lxi+lyi+lzi+kxi+kyi+kzi+nECP(Z,l,w),lambda)

            END DO
            acum= acum + integral*distcoefx * distcoefy * distcoefz *comb(lx,lxi) *comb(ly,lyi) * comb(lz,lzi)
            integral=0.d0
         END DO
      END DO
   END DO
   AAB_LOCAL_loops=acum
END FUNCTION AAB_LOCAL_loops

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

FUNCTION dAAB_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz)
! calcula el termino semi-local del pseudopotencial centrado en i

!  l                                   
!  Σ[<xi|lm> Vi(l-LM) <lm|xj>] 
! m=-l                               

! los coef de la base se multiplican en la rutina que llama a esta

   USE basis_data, ONLY : a !a(i,ni) exponente de la funcion de base i, contrccion ni
   USE ECP_mod, ONLY :ZlistECP,Lmax,nECP,bECP, expnumbersECP,Qnl
   use subm_intECP   , only: comb, OMEGA2, Aintegral, Qtype1N

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

   INTEGER :: l, term, lmaxbase !auxiliades para ciclos
   INTEGER :: Z,n !Z= carga del nucleo
   DOUBLE PRECISION :: acumang, acumint, AABx, AABy, AABz, Kmod,Ccoef

   call g2g_timer_sum_start('ECP_2C_S-local')

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

   CALL AAB_SEMILOCAL_angular(j,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz,Lmax(z)-1) ! angular integrals calculated here and stack in ECP_Ang_stack

   DO l = 0 , Lmax(z)-1 !barre todos los l de la parte no local
      DO term=1, expnumbersECP(z,l) !barre contracciones del ECP para el atomo con carga z y l del ecp
         Ccoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
         Qnl=0.d0
         CALL Qtype1N(Kmod,Ccoef,lmaxbase+l+1,necp(Z,l,term)+n+1,necp(Z,l,term)) !calcula integrales radiales
!              ͚ 
! Qnl(n,l) = ʃ Ml(k*r)*r^n * exp(-cr^2) dr
!  

!FOCK term
         dAAB_SEMILOCAL(1)=dAAB_SEMILOCAL(1)+&
         AAB_SEMILOCAL_loops(j,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz,l,term)
!d/dx term
         dAAB_SEMILOCAL(2)=dAAB_SEMILOCAL(2)+&
         AAB_SEMILOCAL_loops(j,ji,k,lxi,lyi,lzi,lxj+1,lyj,lzj,dx,dy,dz,l,term)*2.d0*a(j,ji)
         IF(lxj.gt.0) dAAB_SEMILOCAL(2)=dAAB_SEMILOCAL(2) &
         -AAB_SEMILOCAL_loops(j,ji,k,lxi,lyi,lzi,lxj-1,lyj,lzj,dx,dy,dz,l,term)*dble(lxj)
!d/dy term
         dAAB_SEMILOCAL(3)=dAAB_SEMILOCAL(3)+&
         AAB_SEMILOCAL_loops(j,ji,k,lxi,lyi,lzi,lxj,lyj+1,lzj,dx,dy,dz,l,term)*2.d0*a(j,ji)
         IF(lyj.gt.0) dAAB_SEMILOCAL(3)=dAAB_SEMILOCAL(3) &
         -AAB_SEMILOCAL_loops(j,ji,k,lxi,lyi,lzi,lxj,lyj-1,lzj,dx,dy,dz,l,term)*dble(lyj)
!d/dz term
         dAAB_SEMILOCAL(4)=dAAB_SEMILOCAL(4)& 
         +AAB_SEMILOCAL_loops(j,ji,k,lxi,lyi,lzi,lxj,lyj,lzj+1,dx,dy,dz,l,term)*2.d0*a(j,ji)
         IF(lzj.gt.0) dAAB_SEMILOCAL(4)=dAAB_SEMILOCAL(4)&
         -AAB_SEMILOCAL_loops(j,ji,k,lxi,lyi,lzi,lxj,lyj,lzj-1,dx,dy,dz,l,term)*dble(lzj)
      END DO
   END DO

   call g2g_timer_sum_pause('ECP_2C_S-local')
   RETURN
END FUNCTION dAAB_SEMILOCAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE AAB_SEMILOCAL_angular(j,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz,LMAX)
   USE basis_data, ONLY : a !a(i,ni) exponente de la funcion de base i, contrccion ni
   USE ECP_mod, ONLY :ZlistECP, ECP_Ang_stack
   use subm_intECP   , only: OMEGA2, Aintegral

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
   INTEGER, INTENT(IN) :: j,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,LMAX
! i,j funciones de la base
! ii,ji numero de contraccion de la funcion
! k atomo con ecp
! lx, ly, lz; i,j exponente de la parte angular de la base x^lx y^ly z^lz
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz ! dx, dy, dz distancia entre nucleos
   DOUBLE PRECISION, DIMENSION(3) :: Kvector

   INTEGER :: m, lx,ly,lz, lambda,l !auxiliades para ciclos
   INTEGER :: Z !Z= carga del nucleo
   DOUBLE PRECISION :: acumang, Kmod !auxiliares
   INTEGER :: lambmin !minimo valor de lambda para integral angular no nula

   INTEGER :: pos1,pos2
   Z=ZlistECP(k)

   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod=2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)

!base
   pos1=lxi+4*lyi+16*lzi+1
   DO lx=0,lxj !barre potencias por expansion del binomio de Newton (x - dx)^lxj
      DO ly=0,lyj
         DO lz=0,lzj
            pos2=lx+4*ly+16*lz+1
            lambmin=0
            DO l=0,LMAX

   IF (l-lx-ly-lz .GT. 0) lambmin=l-lx-ly-lz !minimo valor de lambda para integral angular no nula
   DO lambda=lx+ly+lz+l,lambmin,-2
      acumang=0.d0
      DO m=-l,l
         acumang=acumang+Aintegral(l,m,lxi,lyi,lzi)*OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
      END DO
      ECP_Ang_stack(l,lambda,1,pos1,pos2)=acumang
   END DO

            END DO          
         END DO
      END DO
   END DO

!x+1
   pos1=lxi+4*lyi+16*lzi+1
   lx=lxj+1 !barre potencias por expansion del binomio de Newton (x - dx)^lxj
   DO ly=0,lyj
      DO lz=0,lzj
         pos2=lx+4*ly+16*lz+1
         lambmin=0
         DO l=0,LMAX


   IF (l-lx-ly-lz .GT. 0) lambmin=l-lx-ly-lz !minimo valor de lambda para integral angular no nula
   DO lambda=lx+ly+lz+l,lambmin,-2
      acumang=0.d0
      DO m=-l,l
         acumang=acumang+Aintegral(l,m,lxi,lyi,lzi)*OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
      END DO
      ECP_Ang_stack(l,lambda,1,pos1,pos2)=acumang
   END DO
 

        END DO
      END DO
   END DO

!y+1
   pos1=lxi+4*lyi+16*lzi+1
   ly=lyj+1
   DO lx=0,lxj !barre potencias por expansion del binomio de Newton (x - dx)^lxj
      DO lz=0,lzj
         pos2=lx+4*ly+16*lz+1
         lambmin=0
         DO l=0,LMAX


   IF (l-lx-ly-lz .GT. 0) lambmin=l-lx-ly-lz !minimo valor de lambda para integral angular no nula
   DO lambda=lx+ly+lz+l,lambmin,-2
      acumang=0.d0
      DO m=-l,l
         acumang=acumang+Aintegral(l,m,lxi,lyi,lzi)*OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
      END DO
      ECP_Ang_stack(l,lambda,1,pos1,pos2)=acumang
   END DO


         END DO
      END DO
   END DO

!z+1
!base
   pos1=lxi+4*lyi+16*lzi+1
   DO lx=0,lxj !barre potencias por expansion del binomio de Newton (x - dx)^lxj
      DO ly=0,lyj
         lz=lzj+1
         pos2=lx+4*ly+16*lz+1
         lambmin=0
         DO l=0,LMAX


   IF (l-lx-ly-lz .GT. 0) lambmin=l-lx-ly-lz !minimo valor de lambda para integral angular no nula
   DO lambda=lx+ly+lz+l,lambmin,-2
      acumang=0.d0
      DO m=-l,l
         acumang=acumang+Aintegral(l,m,lxi,lyi,lzi)*OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
      END DO
      ECP_Ang_stack(l,lambda,1,pos1,pos2)=acumang
   END DO


         END DO
      END DO
   END DO


END SUBROUTINE AAB_SEMILOCAL_angular

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

DOUBLE PRECISION FUNCTION AAB_SEMILOCAL_loops(j,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx,dy,dz,l,term)
!parece q tenia un bug en la definicion de lambda, aca deberia estar corregido. TESTEAR!!!!!
   USE basis_data, ONLY : a !a(i,ni) exponente de la funcion de base i, contrccion ni
   USE ECP_mod, ONLY :ZlistECP,aECP,nECP, Qnl, ECP_Ang_stack
   use subm_intECP   , only: comb, OMEGA2, Aintegral, Qtype1N

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
   INTEGER, INTENT(IN) :: j,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,l,term
! i,j funciones de la base
! ii,ji numero de contraccion de la funcion
! k atomo con ecp
! lx, ly, lz; i,j exponente de la parte angular de la base x^lx y^ly z^lz
   DOUBLE PRECISION, INTENT(IN) :: dx,dy,dz ! dx, dy, dz distancia entre nucleos
   DOUBLE PRECISION, DIMENSION(3) :: Kvector

   INTEGER :: lx,ly,lz, lambda !auxiliades para ciclos
   INTEGER :: Z !Z= carga del nucleo
   DOUBLE PRECISION :: acumang, acumint, AABx, AABy, AABz, Kmod, auxdistx,auxdisty,auxdistz !auxiliares
   INTEGER :: lambmin !minimo valor de lambda para integral angular no nula

   INTEGER :: pos1, pos2
#ifdef FULL_CHECKS
   INTEGER :: m
#endif

   AAB_SEMILOCAL_loops=0.d0
   Z=ZlistECP(k)

   Kvector=(/-2.d0*dx,-2.d0*dy,-2.d0*dz/)*a(j,ji)
   Kmod=2.d0 * sqrt(dx**2.d0 + dy**2.d0 + dz**2.d0) *a(j,ji)
   AABx=0.d0
   AABy=0.d0
   AABz=0.d0
   acumint=0.d0
   acumang=0.d0

   pos1=lxi+4*lyi+16*lzi+1
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
                     pos2=lx+4*ly+16*lz+1


   IF (l-lx-ly-lz .GT. 0) lambmin=l-lx-ly-lz !minimo valor de lambda para integral angular no nula

   DO lambda=lx+ly+lz+l,lambmin,-2

#ifdef FULL_CHECKS
      acumang=0.d0
      DO m=-l,l
         acumang=acumang+Aintegral(l,m,lxi,lyi,lzi)*OMEGA2(Kvector,lambda,l,m,lx,ly,lz)
      END DO

      if (ECP_Ang_stack(l,lambda,1,pos1,pos2).ne.acumang) STOP "bad angular int 2C_SLOC"
#endif

      acumint=acumint+ECP_Ang_stack(l,lambda,1,pos1,pos2)*Qnl(necp(Z,l,term)+lx+ly+lz+lxi+lyi+lzi,lambda)*aECP(z,L,term)
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    3 Center terms <Xb|Va|Xc>    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


FUNCTION dABC_LOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY :expnumbersECP, Qnl,bECP,IzECP, Lmax,necp
   use subm_intECP   , only: OMEGA1, Q0, comb, Qtype1N, Qtype0N

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji !terminos de la base
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dx1,dy1,dz1,dx2,dy2,dz2 !distancias del centro con ecp a cada nucleo
   INTEGER, INTENT(IN) :: k !numero de atomo
   INTEGER :: L, Z !L maximo del ecp, Z carga nuclear sin modificar
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   DOUBLE PRECISION :: Kmod,Ccoef, integral,acum
   INTEGER :: lmaxbase,w
   DOUBLE PRECISION, DIMENSION(7) :: dABC_LOCAL !(<A|B|C>,<dA/dxA|B|C>,<dA/dyA|B|C>,<dA/dzA|B|C>,<A|B|dC/dxC>,<A|B|dC/dyC>,<A|B|dC/dzC>

   call g2g_timer_sum_start('ECP_3C_local')

   dABC_LOCAL=0.d0
   z=IzECP(k)
   L=Lmax(Z)
   lmaxbase=lxi+lyi+lzi+lxj+lyj+lzj
   Kvector=(/a(i,ii)*dx1+a(j,ji)*dx2,a(i,ii)*dy1+a(j,ji)*dy2,a(i,ii)*dz1+a(j,ji)*dz2/)
   Kvector=-2.d0*Kvector
   Kmod=sqrt(Kvector(1)**2.d0+Kvector(2)**2.d0+Kvector(3)**2.d0)
   integral=0.d0
   acum=0.d0

   call ABC_LOCAL_angular(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2) ! angular integrals calculated here and stack in ECP_Ang_stack


   DO w =1, expnumbersECP(z,l) !barre terminos del ECP para el atomo con carga nuclear Z y l del ECP
      Qnl=0.d0
      Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)

!calcula integral radial
      IF ( Kmod .GT. 0.d0 ) THEN
         CALL Qtype1N(Kmod,Ccoef,Lmaxbase+1,necp(Z,l,w)+Lmaxbase+3,necp(Z,l,w)) 
      ELSE
         CALL Qtype0N(lmaxbase+1+nECP(Z,l,w),lmaxbase+1,Ccoef)
!parche para el caso accidental en que K fuera = (0,0,0) por compensacion de a(i,ii)*dx1+a(j,ji)*dx2 (idem y,z)
!0.25d0/pi compensa el factor 4pi por el que se multiplica luego a la suma de las integrales
      END IF


!FOCK term
      dABC_LOCAL(1)=dABC_LOCAL(1)+ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)
!d/dxi term
      dABC_LOCAL(2)=dABC_LOCAL(2)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi+1,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(i,ii)
      if(lxi.ge.1) dABC_LOCAL(2)=dABC_LOCAL(2)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi-1,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lxi)
!d/dyi term
      dABC_LOCAL(3)=dABC_LOCAL(3)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi+1,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(i,ii)
      if(lyi.ge.1) dABC_LOCAL(3)=dABC_LOCAL(3)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi-1,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lyi)
!d/dzi term
      dABC_LOCAL(4)=dABC_LOCAL(4)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi+1,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(i,ii)
      if(lzi.ge.1) dABC_LOCAL(4)=dABC_LOCAL(4)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi-1,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lzi)
!d/dxj term
      dABC_LOCAL(5)=dABC_LOCAL(5)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj+1,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(j,ji)
      if(lxj.ge.1) dABC_LOCAL(5)=dABC_LOCAL(5)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj-1,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lxj)
!d/dyj term
      dABC_LOCAL(6)=dABC_LOCAL(6)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj+1,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(j,ji)
      if(lyj.ge.1) dABC_LOCAL(6)=dABC_LOCAL(6)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj-1,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lyj)
!d/dzj term
      dABC_LOCAL(7)=dABC_LOCAL(7)+&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj+1,dx1,dy1,dz1,dx2,dy2,dz2,w)*2.d0*a(j,ji)
      if(lzj.ge.1) dABC_LOCAL(7)=dABC_LOCAL(7)-&
      ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj-1,dx1,dy1,dz1,dx2,dy2,dz2,w)*dble(lzj)
   END DO

   call g2g_timer_sum_pause('ECP_3C_local')
   RETURN
END FUNCTION dABC_LOCAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE ABC_LOCAL_angular(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : IzECP,angularint, Fulltimer_ECP,Lmax, ECP_Ang_stack
   use subm_intECP   , only: OMEGA1, Q0, comb, Qtype1N

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji !terminos de la base
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dx1,dy1,dz1,dx2,dy2,dz2 !distancias del centro con ecp a cada nucleo
   INTEGER, INTENT(IN) :: k !numero de atomo
   INTEGER :: L, Z !L maximo del ecp, Z carga nuclear sin modificar
   DOUBLE PRECISION,DIMENSION (3) :: Kvector
   DOUBLE PRECISION :: Kmod, integral,acum
   INTEGER :: lmaxbase
   INTEGER :: ac,bc,cc,dc,ec,fc,lambda ! variables auxiliares
   DOUBLE PRECISION :: t1 ! auxiliares para timers
   INTEGER :: pos1, pos2
   IF (Fulltimer_ECP) CALL cpu_time ( t1 )


   z=IzECP(k)
   L=Lmax(Z)
   lmaxbase=lxi+lyi+lzi+lxj+lyj+lzj
   Kvector=(/a(i,ii)*dx1+a(j,ji)*dx2,a(i,ii)*dy1+a(j,ji)*dy2,a(i,ii)*dz1+a(j,ji)*dz2/)
   Kvector=-2.d0*Kvector
   Kmod=sqrt(Kvector(1)**2.d0+Kvector(2)**2.d0+Kvector(3)**2.d0)
   integral=0.d0
   acum=0.d0

   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
      DO bc=0,lyi
         DO cc=0,lzi
            DO dc=0,lxj
               DO ec=0,lyj
                  DO fc=0,lzj
   pos1=ac+4*bc+16*cc+1
   pos2=dc+4*ec+16*fc+1

   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
      IF ( Kmod .GT. 0.d0 ) THEN
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc)
      ELSE
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= angularint(ac+dc,bc+ec,cc+fc)
      END IF
   END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

!1
  ac=lxi+1 !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
  DO bc=0,lyi
     DO cc=0,lzi
        pos1=ac+4*bc+16*cc+1
        DO dc=0,lxj
           DO ec=0,lyj
              DO fc=0,lzj
   pos2=dc+4*ec+16*fc+1

   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
      IF ( Kmod .GT. 0.d0 ) THEN
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc)
      ELSE
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= angularint(ac+dc,bc+ec,cc+fc)
      END IF
   END DO
               END DO
            END DO
         END DO
      END DO
   END DO

!2
   bc=lyi+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
      DO cc=0,lzi
         pos1=ac+4*bc+16*cc+1
         DO dc=0,lxj
            DO ec=0,lyj
               DO fc=0,lzj
   pos2=dc+4*ec+16*fc+1

   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
      IF ( Kmod .GT. 0.d0 ) THEN
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc)
      ELSE
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= angularint(ac+dc,bc+ec,cc+fc)
      END IF
   END DO
               END DO
            END DO
         END DO
      END DO
   END DO
!3
   cc=lzi+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
      DO bc=0,lyi
         pos1=ac+4*bc+16*cc+1
         DO dc=0,lxj
            DO ec=0,lyj
               DO fc=0,lzj
   pos2=dc+4*ec+16*fc+1

   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
      IF ( Kmod .GT. 0.d0 ) THEN
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc)
      ELSE
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= angularint(ac+dc,bc+ec,cc+fc)
      END IF
   END DO
               END DO
            END DO
         END DO
      END DO
   END DO
!4
   dc=lxj+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
      DO bc=0,lyi
         DO cc=0,lzi
            pos1=ac+4*bc+16*cc+1
            DO ec=0,lyj
               DO fc=0,lzj
                  pos2=dc+4*ec+16*fc+1

   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
      IF ( Kmod .GT. 0.d0 ) THEN
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc)
      ELSE
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= angularint(ac+dc,bc+ec,cc+fc)
      END IF
   END DO
               END DO
            END DO
         END DO
      END DO
   END DO
!5
   ec=lyj+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
      DO bc=0,lyi
         DO cc=0,lzi
            pos1=ac+4*bc+16*cc+1
            DO dc=0,lxj
               DO fc=0,lzj
                  pos2=dc+4*ec+16*fc+1

   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
      IF ( Kmod .GT. 0.d0 ) THEN
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc)
      ELSE
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= angularint(ac+dc,bc+ec,cc+fc)
      END IF
   END DO
               END DO
            END DO
         END DO
      END DO
   END DO
!6
   fc=lzj+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
      DO bc=0,lyi
         DO cc=0,lzi
            pos1=ac+4*bc+16*cc+1
            DO dc=0,lxj
               DO ec=0,lyj
                  pos2=dc+4*ec+16*fc+1

   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
      IF ( Kmod .GT. 0.d0 ) THEN
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc)
      ELSE
         ECP_Ang_stack(0,lambda,0,pos1,pos2)= angularint(ac+dc,bc+ec,cc+fc)
      END IF
   END DO
               END DO
            END DO
         END DO
      END DO
   END DO

   RETURN
END SUBROUTINE ABC_LOCAL_angular

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

DOUBLE PRECISION FUNCTION ABC_LOCAL_loops(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dx1,dy1,dz1,dx2,dy2,dz2,w)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : Qnl,IzECP, Fulltimer_ECP,Lmax,necp,aECP, ECP_Ang_stack,  &
#ifdef FULL_CHECKS
   angularint, pi,  &
#endif
   bECP

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
   INTEGER :: pos1, pos2
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
   Ccoef=bECP(z,L,w)+a(i,ii)+a(j,ji)

   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dx1)^lxi
      auxdista=dx1**(lxi-ac)
      DO bc=0,lyi
         auxdistb=dy1**(lyi-bc)
         DO cc=0,lzi
            auxdistc=dz1**(lzi-cc)
            pos1=ac+4*bc+16*cc+1
            DO dc=0,lxj
               auxdistd=dx2**(lxj-dc)
               DO ec=0,lyj
                  auxdiste=dy2**(lyj-ec)
                  DO fc=0,lzj
                     auxdistf=dz2**(lzj-fc)
                     pos2=dc+4*ec+16*fc+1

   auxdist=auxdista*auxdistb*auxdistc*auxdistd*auxdiste*auxdistf
   DO lambda=ac+bc+cc+dc+ec+fc,0,-2
#ifdef FULL_CHECKS
      IF ( Kmod .GT. 0.d0 ) THEN
          IF (OMEGA1(Kvector,lambda,ac+dc,bc+ec,cc+fc).ne.ECP_Ang_stack(0,lambda,0,pos1,pos2)) &
          STOP "bad angular int 3C_LOC_I"
      ELSE
          IF (angularint(ac+dc,bc+ec,cc+fc).ne.ECP_Ang_stack(0,lambda,0,pos1,pos2)) STOP "bad angular int 3C_LOC_E"
          IF (Q0(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),Ccoef) *0.25d0/pi .ne. Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)) &
          STOP "bad radial int 3C_LOC_E"
      END IF
#endif
      integral=integral +ECP_Ang_stack(0,lambda,0,pos1,pos2)*Qnl(ac+bc+cc+dc+ec+fc+nECP(Z,l,w),lambda)
   END DO

   auxcomb=comb(lxi,ac)*comb(lyi,bc)*comb(lzi,cc)*comb(lxj,dc)*comb(lyj,ec)*comb(lzj,fc)
   acum=acum + auxcomb*auxdist*integral

   integral=0.d0
   auxcomb=0.d0
   auxdist=0.d0

                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
   ABC_LOCAL_loops=aECP(z,L,w)*acum
   RETURN
END FUNCTION ABC_LOCAL_loops

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

FUNCTION dABC_SEMILOCAL(i,j,ii,ji,k,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : Qnl1l2,necp,bECP,IzECP,Fulltimer_ECP,Lmax,expnumbersECP, usedQnl1l2
   use subm_intECP   , only: Qtype2N
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
   DOUBLE PRECISION :: t1q,t1aux,t2aux !auxiliares para timers

   call g2g_timer_sum_start('ECP_3C_S-local')
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

   call ABC_SEMILOCAL_angular(i,j,ii,ji,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,Lmax(z)-1) ! angular integrals calculated here and stack in ECP_Ang_stack

   DO l = 0 , Lmax(z)-1 !barre todos los l de la parte no local
      DO term=1, expnumbersECP(Z,l) !barre contracciones del ECP para el atomo con carga Z y momento angular l del ecp
         Qnl1l2=0.d0
         Ccoef=bECP(z,L,term)+a(i,ii)+a(j,ji)
         IF (Fulltimer_ECP) CALL cpu_time ( t1q )
         call Qtype2N(Kimod,Kjmod,Ccoef,l1max+l+1,l2max+l+1,necp(Z,l,term)+l1max+l2max+1,necp(Z,l,term)) !agrega a la matriz Qnl1l2 los terminos correspondientes a un termino radiales.

   usedQnl1l2=.false.

!Fock
   dABC_SEMILOCAL(1)=dABC_SEMILOCAL(1)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)
!d/dxi
   dABC_SEMILOCAL(2)=dABC_SEMILOCAL(2)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi+1,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(i,ii)
   if(lxi.ge.1) dABC_SEMILOCAL(2)=dABC_SEMILOCAL(2)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi-1,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lxi)
!d/dyi
   dABC_SEMILOCAL(3)=dABC_SEMILOCAL(3)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi+1,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(i,ii)
   if(lyi.ge.1) dABC_SEMILOCAL(3)=dABC_SEMILOCAL(3)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi-1,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lyi)
!d/dzi
   dABC_SEMILOCAL(4)=dABC_SEMILOCAL(4)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi+1,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(i,ii)
   if(lzi.ge.1) dABC_SEMILOCAL(4)=dABC_SEMILOCAL(4)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi-1,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lzi)
!d/dxj
   dABC_SEMILOCAL(5)=dABC_SEMILOCAL(5)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi,lxj+1,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(j,ji)
   if(lxj.ge.1) dABC_SEMILOCAL(5)=dABC_SEMILOCAL(5)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi,lxj-1,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lxj)
!d/dyj
   dABC_SEMILOCAL(6)=dABC_SEMILOCAL(6)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi,lxj,lyj+1,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(j,ji)
   if(lyj.ge.1) dABC_SEMILOCAL(6)=dABC_SEMILOCAL(6)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi,lxj,lyj-1,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lyj)
!d/dzj
   dABC_SEMILOCAL(7)=dABC_SEMILOCAL(7)+ &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi,lxj,lyj,lzj+1,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*2.d0*a(j,ji)
   if(lzj.ge.1) dABC_SEMILOCAL(7)=dABC_SEMILOCAL(7)- &
   ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi,lxj,lyj,lzj-1,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)*dble(lzj)


      END DO
   END DO

   call g2g_timer_sum_pause('ECP_3C_S-local')

   RETURN
END FUNCTION dABC_SEMILOCAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

SUBROUTINE ABC_SEMILOCAL_angular(i,j,ii,ji,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,lMAX)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : ECP_Ang_stack
   use subm_intECP   , only: OMEGA2
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji
!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dxi,dyi,dzi,dxj,dyj,dzj !distancias de las nucleos con bases al nucleo con ecp
   INTEGER, INTENT(IN) :: lMAX
   INTEGER :: l1max,l2max,l !Z= carga del nucleo, y momento angular de la base i y j
   DOUBLE PRECISION,DIMENSION (3) :: Kivector,Kjvector
   INTEGER :: ac,bc,cc,dc,ec,fc,lambdai, lambdaj,m,lambimin, lambjmin !auxiliares ciclos
   INTEGER :: pos1, pos2
   DOUBLE PRECISION :: acumang
   l1max=lxi+lyi+lzi
   l2max=lxj+lyj+lzj

   Kivector=-2.d0*a(i,ii)*(/dxi,dyi,dzi/)
   Kjvector=-2.d0*a(j,ji)*(/dxj,dyj,dzj/)
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dxi)^lxi
      DO bc=0,lyi
         DO cc=0,lzi
            pos1=ac+4*bc+16*cc+1
            DO dc=0,lxj
               DO ec=0,lyj
                  DO fc=0,lzj
                     pos2=dc+4*ec+16*fc+1
                     DO l=0,LMAX
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
         ECP_Ang_stack(l,lambdai,lambdaj,pos1,pos2)=acumang
         acumang=0.d0
      END DO
   END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

!1+1
   ac=lxi+1 !barre potencias por expansion del binomio de Newton (x - dxi)^lxi
   DO bc=0,lyi
      DO cc=0,lzi
         pos1=ac+4*bc+16*cc+1
         DO dc=0,lxj
            DO ec=0,lyj
               DO fc=0,lzj
                  pos2=dc+4*ec+16*fc+1
                  DO l=0,LMAX
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
         ECP_Ang_stack(l,lambdai,lambdaj,pos1,pos2)=acumang
         acumang=0.d0
      END DO
   END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

!2+
   bc=lyi+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dxi)^lxi
      DO cc=0,lzi
         pos1=ac+4*bc+16*cc+1
         DO dc=0,lxj
            DO ec=0,lyj
               DO fc=0,lzj
                 pos2=dc+4*ec+16*fc+1
                 DO l=0,LMAX
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
         ECP_Ang_stack(l,lambdai,lambdaj,pos1,pos2)=acumang
         acumang=0.d0
      END DO
   END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

!3+1
   cc=lzi+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dxi)^lxi
      DO bc=0,lyi
         pos1=ac+4*bc+16*cc+1
         DO dc=0,lxj
            DO ec=0,lyj
               DO fc=0,lzj
                  pos2=dc+4*ec+16*fc+1
                  DO l=0,LMAX
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
         ECP_Ang_stack(l,lambdai,lambdaj,pos1,pos2)=acumang
         acumang=0.d0
      END DO
   END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

!4+1
   dc=lxj+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dxi)^lxi
      DO bc=0,lyi
         DO cc=0,lzi
            pos1=ac+4*bc+16*cc+1
            DO ec=0,lyj
               DO fc=0,lzj
                  pos2=dc+4*ec+16*fc+1
                  DO l=0,LMAX
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
         ECP_Ang_stack(l,lambdai,lambdaj,pos1,pos2)=acumang
         acumang=0.d0
      END DO
   END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

!5+1
   ec=lyj+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dxi)^lxi
      DO bc=0,lyi
         DO cc=0,lzi
            pos1=ac+4*bc+16*cc+1
            DO dc=0,lxj
               DO fc=0,lzj
                 pos2=dc+4*ec+16*fc+1
                 DO l=0,LMAX
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
         ECP_Ang_stack(l,lambdai,lambdaj,pos1,pos2)=acumang
         acumang=0.d0
      END DO
   END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

!6+1
   fc=lzj+1
   DO ac=0,lxi !barre potencias por expansion del binomio de Newton (x - dxi)^lxi
      DO bc=0,lyi
         DO cc=0,lzi
            pos1=ac+4*bc+16*cc+1
            DO dc=0,lxj
               DO ec=0,lyj
                  pos2=dc+4*ec+16*fc+1
                  DO l=0,LMAX
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
         ECP_Ang_stack(l,lambdai,lambdaj,pos1,pos2)=acumang
         acumang=0.d0
      END DO
   END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
   RETURN
END SUBROUTINE ABC_SEMILOCAL_angular

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

DOUBLE PRECISION FUNCTION ABC_SEMILOCAL_loops(i,j,ii,ji,lxi,lyi,lzi,lxj,lyj,lzj,dxi,dyi,dzi,dxj,dyj,dzj,l,term,Z)
   USE basis_data, ONLY : a
   USE ECP_mod, ONLY : Qnl1l2,necp,Fulltimer_ECP, aECP, ECP_Ang_stack 
   use subm_intECP   , only: comb,OMEGA2
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,ii,ji
!i,j funciones de la base
!ii,ji numero de contraccion de la funcion
!k atomo con ecp
   INTEGER, INTENT(IN) :: lxi,lyi,lzi,lxj,lyj,lzj !potencias de la parte angular
   DOUBLE PRECISION, INTENT(IN) :: dxi,dyi,dzi,dxj,dyj,dzj !distancias de las nucleos con bases al nucleo con ecp
   INTEGER, INTENT(IN) :: l,term,Z
   INTEGER :: l1max,l2max !Z= carga del nucleo, y momento angular de la base i y j
   DOUBLE PRECISION,DIMENSION (3) :: Kivector,Kjvector
   DOUBLE PRECISION :: Kimod, Kjmod
   INTEGER :: ac,bc,cc,dc,ec,fc,lambdai, lambdaj,lambimin, lambjmin !auxiliares ciclos
   DOUBLE PRECISION :: acumang,acumang1,acumang2,integral,auxcomb,auxdist,acum,auxdista,auxdistb,auxdistc,auxdistd,auxdiste,auxdistf!auxiliares
   DOUBLE PRECISION :: t1aux,t2aux !auxiliares para timers
   INTEGER :: pos1, pos2

#ifdef FULL_CHECKS
   INTEGER :: m
#endif

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
      DO bc=0,lyi
         auxdistb=dyi**(lyi-bc)
         DO cc=0,lzi
            auxdistc=dzi**(lzi-cc)
            pos1=ac+4*bc+16*cc+1
            DO dc=0,lxj
               auxdistd=dxj**(lxj-dc)
               DO ec=0,lyj
                  auxdiste=dyj**(lyj-ec)
                  DO fc=0,lzj
                     auxdistf=dzj**(lzj-fc)
                     pos2=dc+4*ec+16*fc+1
   IF (Fulltimer_ECP) CALL cpu_time ( t1aux )

   auxdist=auxdista*auxdistb*auxdistc*auxdistd*auxdiste*auxdistf
   auxcomb=comb(lxi,ac)*comb(lyi,bc)*comb(lzi,cc)*comb(lxj,dc)*comb(lyj,ec)*comb(lzj,fc)
   lambimin=0
   lambjmin=0

   IF (l-ac-bc-cc .GT. 0) lambimin=l-ac-bc-cc !lambda minimo que no anula la integral angular
   IF (l-dc-ec-fc .GT. 0) lambjmin=l-dc-ec-fc !lambda minimo que no anula la integral angular

   DO lambdai=ac+bc+cc+l,lambimin,-2
      DO lambdaj=dc+ec+fc+l,lambjmin,-2

#ifdef FULL_CHECKS
         acumang=0.d0
         DO m=-l,l
            acumang=acumang+OMEGA2(Kivector,lambdai,l,m,ac,bc,cc)*OMEGA2(Kjvector,lambdaj,l,m,dc,ec,fc)
         END DO
         if ( ECP_Ang_stack(l,lambdai,lambdaj,pos1,pos2)-acumang .ne. 0.d0) STOP "bad angular int 3C_SL"
#endif

         integral=integral+ECP_Ang_stack(l,lambdai,lambdaj,pos1,pos2) &
         *Qnl1l2(ac+bc+cc+dc+ec+fc+necp(Z,l,term),lambdai,lambdaj)

         acumang=0.d0
      END DO
   END DO
   acum=acum + auxcomb*auxdist*integral
   integral=0.d0
   auxcomb=0.d0
   auxdist=0.d0
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

   ABC_SEMILOCAL_loops=acum*aECP(Z,L,term)

   RETURN
END FUNCTION ABC_SEMILOCAL_loops

end module subm_intECPG

