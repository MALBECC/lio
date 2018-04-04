!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dftb_subs
   implicit none

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dftb_init(M, OPEN)

   use dftb_data, only: MTB, MDFTB, end_bTB, Iend_TB, rho_aDFTB, rho_bDFTB

   implicit none

   integer, intent(in) :: M
   logical, intent(in) :: OPEN

   MDFTB=2*MTB+M
   allocate(Iend_TB(2,2*end_bTB), rho_aDFTB(MDFTB,MDFTB))
   if (OPEN) allocate (rho_bDFTB(MDFTB,MDFTB))

end subroutine dftb_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dftb_td_init (M,rho, rho_0, overlap, RMM5, dim3)

   use dftb_data, only: MTB, MDFTB, rho_aDFTB, rho_bDFTB, rhold_AOTB,          &
                        rhonew_AOTB
   implicit none
   integer, intent(in)       :: M
!carlos: dim3 is to declare the 3th dimension of matrices
   integer, intent(in)       :: dim3
   real*8, allocatable, intent(inout) :: overlap(:,:)
   real*8 , intent(in)  :: RMM5(M*(M+1)/2)
#ifdef TD_SIMPLE
   complex*8, intent(in)     :: rho_0(M,M,dim3)
   complex*8, intent(out)  :: rho(MDFTB,MDFTB,dim3)
#else
   complex*16, intent(in)    :: rho_0(M,M,dim3)
   complex*16, intent(out) :: rho(MDFTB,MDFTB,dim3)
#endif
   integer  :: ii,jj

   if (allocated(overlap)) deallocate(overlap)
   allocate(overlap(M,M), rhold_AOTB(MDFTB,MDFTB,dim3),                       &
            rhonew_AOTB(MDFTB,MDFTB,dim3))

   rhold_AOTB=0.0d0
   rhonew_AOTB=0.0d0

   call spunpack('L', M, RMM5, overlap)

      do jj=1, MDFTB
      do ii=1, MDFTB
         if (ii==jj) then
            rho(ii,jj,1)=cmplx(rho_aDFTB(ii,jj), 0.0D0)
         else
            rho(ii,jj,1)=cmplx(rho_aDFTB(ii,jj)/2.0d0, 0.0D0)
         end if
      end do
      end do

!carlos:Open shell option
   if(dim3==2) then
      do jj=1, MDFTB
      do ii=1, MDFTB
         if (ii==jj) then
            rho(ii,jj,2)=cmplx(rho_bDFTB(ii,jj), 0.0D0)
         else
            rho(ii,jj,2)=cmplx(rho_bDFTB(ii,jj)/2.0d0, 0.0D0)
         end if
      end do
      end do
   end if

   do jj= 1, M
   do ii= 1, M
      rho(ii+MTB,jj+MTB,:)=rho_0(ii,jj,:)
   end do
   end do

end subroutine dftb_td_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%i%%%%%%%%%%%%%%%%%%%%%%%%

subroutine getXY_DFTB(M_in,x_in,y_in,xmat,ymat)

   use dftb_data, only: MTB, MDFTB

   implicit none
   integer, intent(in)  :: M_in
   real*8 , intent(in)  :: x_in(M_in,M_in)
   real*8 , intent(in)  :: y_in(M_in,M_in)
   real*8 , intent(out) :: xmat(MDFTB,MDFTB)
   real*8 , intent(out) :: ymat(MDFTB,MDFTB)
   integer              :: ii,jj


   xmat=0.0d0
   ymat=0.0d0

   do ii=1, MTB
      xmat(ii,ii)=1.0d0
      xmat(MTB+M_in+ii,MTB+M_in+ii)=1.0d0

      ymat(ii,ii)=1.0d0
      ymat(MTB+M_in+ii,MTB+M_in+ii)=1.0d0
   end do

   do jj=1, M_in
   do ii=1, M_in
      xmat(MTB+ii, MTB+jj)=x_in(ii,jj)
      ymat(MTB+ii, MTB+jj)=y_in(ii,jj)
   end do
   end do
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine read_rhoDFTB(M, MM, RMM1, rhoalpha, rhobeta, OPEN)
   use dftb_data, only: MTB, MDFTB, rho_aDFTB, rho_bDFTB

   implicit none
   integer, intent(in)    :: M, MM
   real*8, intent(inout)  :: RMM1(MM)
   real*8, intent(inout)  :: rhoalpha(MM), rhobeta(MM)
   logical, intent(in)    :: OPEN
   integer                :: ii, jj, kk


  open(unit=1070,file='rhoTB.in')

   if (.not.OPEN) then
      do ii=1,MDFTB
      do jj=1,MDFTB
         read(1070,*) rho_aDFTB(ii,jj)
      end do
      end do

      do jj=1,M
      do kk=jj,M
               RMM1(kk+(2*M-jj)*(jj-1)/2)=rho_aDFTB(jj+MTB,kk+MTB)
      enddo
      enddo

   else if(OPEN) then

      do ii=1,MDFTB
      do jj=1,MDFTB
         read(1070,*) rho_aDFTB(ii,jj), rho_bDFTB(ii,jj)
      end do
      end do

      do jj=1,M
      do kk=jj,M
               rhoalpha(kk+(2*M-jj)*(jj-1)/2)=rho_aDFTB(jj+MTB,kk+MTB)
               rhobeta(kk+(2*M-jj)*(jj-1)/2)=rho_aDFTB(jj+MTB,kk+MTB)
      enddo
      enddo

      RMM1 = rhoalpha+rhobeta

      write(*,*) 'RHOTB readed'

   end if

    close(1070)
end subroutine read_rhoDFTB
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine construct_rhoDFTB(M, rho, rho_0 ,rho_DFTB, TBload, niter)

   use dftb_data , only: MTB, MDFTB

   implicit none
   integer, intent(in) :: M
   integer, intent(in) :: niter
   real*8, intent(in)  :: rho_0(M,M)
   real*8, intent(in)  :: rho_DFTB(MDFTB,MDFTB)
   logical, intent(in) :: TBload
   real*8, intent(out) :: rho(MDFTB,MDFTB)
   integer :: ii, jj

   if ((TBload.and.niter==1).or.(niter/=1)) then

      do ii=1,MDFTB
      do jj=ii+1,MDFTB
         rho(ii,jj)=rho_DFTB(ii,jj)/2
         rho(jj,ii)=rho(ii,jj)
      end do
      end do

   else if(.not.TBload.and.(niter==1)) then

      rho=0.0D0
      do ii=1, MTB
         rho(ii,ii)=1.0D0
         rho(MTB+M+ii,MTB+M+ii)=1.0D0
      end do
         rho(MTB+1:MTB+M, MTB+1:MTB+M)=rho_0
   end if

end subroutine construct_rhoDFTB
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine find_TB_neighbors(M_in,Nuc,natom)

   use dftb_data, only:Iend_TB

   implicit none

   integer              ::  ii, jj
   integer, intent(in)  ::  M_in, natom
   integer, intent(in)  ::  Nuc(M_in)


   jj=0
   do ii=1, M_in
      if (Nuc(ii)==1.or.Nuc(ii)==natom) then
         jj=jj+1
         Iend_TB(1,jj) = Nuc(ii)
         Iend_TB(2,jj) = ii
      end if
   end do

end subroutine find_TB_neighbors

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine build_chimera_DFTB (M_in,fock_in, fock_DFTB, natom, nshell, ncont)

   use dftb_data, only:MDFTB, MTB, Iend_TB, end_bTB, alfaTB, betaTB, gammaTB,Vbias_TB

   integer, intent(in)  :: M_in
   integer, intent(in)  :: natom
   integer, intent(in)  :: nshell (0:4)
   integer, intent(in)  :: ncont(M_in)
   real*8, intent(in)   :: fock_in (M_in, M_in)
   real*8, intent(out)  :: fock_DFTB (MDFTB, MDFTB)
   integer              :: ii, jj, kk, ns, np, l1, l2

   l1=0
   l2=0
   jj=0
   kk=0
   ns=nshell(0)
   np=nshell(1)

   fock_DFTB(:,:) = 0.0D0

   do ii=1, 2*end_bTB

      if (Iend_TB(1, ii) == 1) then

!         if ((ncont(Iend_TB(2,ii))==1).or.(ncont(Iend_TB(2,ii))==2)) then

            if (Iend_TB(2, ii)>ns.and.Iend_TB(2, ii)<=ns+np) then
               jj=jj+1
               if (jj==1) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =0.0d0 !-gammaTB
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =0.0d0 !-gammaTB
               else if(jj==2) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  = 0.0D0 !gammaTB !0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) = 0.0D0 !gammaTB !0.0d0
               else if(jj==3) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  = 0.0D0 !gammaTB !0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) = 0.0D0 !gammaTB !0.0d0
                  jj=0
               end if

            else if (Iend_TB(2, ii)>ns+np) then
               kk=kk+1
               if (kk==1) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =gammaTB
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =gammaTB
               else if (kk==2) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =gammaTB !0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =gammaTB !0.0d0
               else if (kk==3) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =gammaTB
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =gammaTB
               else if (kk==4) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =gammaTB !0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =gammaTB !0.0d0
               else if (kk==5) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =gammaTB !0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =gammaTB !0.0d0
               else if (kk==6) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =gammaTB
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =gammaTB
                  kk=0
               endif
            else
               fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  = gammaTB
               fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) = gammaTB

            end if
!            l1=l1+1

!         else

!            fock_DFTB(Iend_TB(2,ii)+MTB,MTB) = 0.0D0
!            fock_DFTB(MTB, Iend_TB(2,ii)+MTB) = 0.0D0

!         end if

      else if (Iend_TB(1, ii) == natom) then
            if (Iend_TB(2, ii)>ns.and.Iend_TB(2, ii)<=ns+np) then
               jj=jj+1
               if (jj==1) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  =0.0d0 !gammaTB
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) =0.0d0 !gammaTB
               else if(jj==2) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  =0.0D0 !gammaTB ! 0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) =0.0D0 !gammaTB !0.0d0
               else if(jj==3) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  =0.0D0 !gammaTB !0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) =0.0D0 !gammaTB !0.0d0
                  jj=0
               end if
            else if (Iend_TB(2, ii)>ns+np) then
               kk=kk+1
               if (kk==1) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB
               else if (kk==2) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB !0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB !0.0d0
               else if (kk==3) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB
               else if (kk==4) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB !0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB !0.0d0
               else if (kk==5) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB !0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB !0.0d0
               else if (kk==6) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB
                  kk=0
               endif
            else
               fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
               fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB

            end if
!            l2=l2+1
!         if ((ncont(Iend_TB(2,ii))==1) .or. (ncont(Iend_TB(2,ii))==2)) then

!            fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1)  = gammaTB
!            fock_DFTB(MTB+M_in+1, Iend_TB(2,ii)+MTB) = gammaTB
!         else
!            fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1)  = gammaTB
!            fock_DFTB(MTB+M_in+1, Iend_TB(2,ii)+MTB) = gammaTB
!         end if

      end if

   end do

   do ii=1,MTB
      fock_DFTB(ii,ii) = alfaTB-Vbias_TB/2.0d0
      fock_DFTB(MTB+M_in+ii, MTB+M_in+ii)= alfaTB+Vbias_TB/2.0d0

      if (ii<MTB) then !-end_bTB) then

         fock_DFTB(ii,ii+1) = betaTB
         fock_DFTB(ii+1,ii) = betaTB
         fock_DFTB(2*MTB+M_in-ii, 2*MTB+M_in-ii+1)= betaTB
         fock_DFTB(2*MTB+M_in-ii+1, 2*MTB+M_in-ii)= betaTB

      end if
   end do

!   do ii=1, end_bTB

!      fock_DFTB(MTB-end_bTB, MTB-end_bTB+ii)=betaTB
!      fock_DFTB(MTB-end_bTB+ii, MTB-end_bTB)=betaTB
!      fock_DFTB(MTB+M_in+end_bTB+1, MTB+M_in+end_bTB+1-ii)=betaTB
!      fock_DFTB(MTB+M_in+end_bTB+1-ii,MTB+M_in+end_bTB+1)=betaTB

!   end do

   fock_DFTB(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)


end subroutine build_chimera_DFTB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine extract_rhoDFT (M_in, rho_in, rho_out)

   use dftb_data, only: MDFTB, MTB

   implicit none
   integer, intent(in)  :: M_in
   real*8, intent(in)   :: rho_in(MDFTB,MDFTB)
   real*8, intent(out)  :: rho_out(M_in,M_in)
   integer              :: ii, jj

   rho_out=0.0D0

   do jj=1, M_in
   do ii=1, M_in
      rho_out(ii,jj)=rho_in(MTB+ii,MTB+jj)
   end do
   end do

end subroutine extract_rhoDFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine chimeraDFTB_evol(M_in,fock_in, fock_DFTB, natom, nshell,ncont, istep)

   use dftb_data, only:MDFTB, MTB, Iend_TB, end_bTB, alfaTB, betaTB, gammaTB,   &
                       Vbias_TB, start_tdtb, end_tdtb

   integer, intent(in)  :: M_in
   integer, intent(in)  :: natom
   integer, intent(in)  :: ncont(M_in)
   integer, intent(in)  :: istep
   integer, intent(in)  :: nshell (0:4)
   real*8, intent(in)   :: fock_in (M_in, M_in)
   real*8, intent(out)  :: fock_DFTB (MDFTB, MDFTB) !temporal dimensions
   real*8               :: pi=4.0D0*atan(1.0D0)
   real*8               :: lambda, t_step, f_t
   integer              :: ii, jj, kk, ns, np, l1, l2


   l1=0
   l2=0
   jj=0
   kk=0
   ns=nshell(0)
   np=nshell(1)

   print*, "istep, start, end", istep, start_tdtb, end_tdtb

   lambda=1.0d0/real(end_tdtb-start_tdtb)

   if (istep < start_tdtb) then
      f_t=1.0D0
   else if(istep >= start_tdtb .and. istep < end_tdtb) then
      t_step=real(istep-start_tdtb)
      f_t=(Cos(pi*lambda*t_step)+1.0D0)/2.0D0
   else if(istep >= end_tdtb) then
      f_t=0.0D0
   end if

   print*, "factor V", f_t

   fock_DFTB(:,:) = 0.0D0

   do ii=1, 2*end_bTB
      if (Iend_TB(1, ii) == 1) then
!         if ((ncont(Iend_TB(2,ii))==1) .or. (ncont(Iend_TB(2,ii))==2)) then

            if (Iend_TB(2, ii)>ns.and.Iend_TB(2, ii)<=ns+np) then
               jj=jj+1
               if (jj==1) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =0.0d0  !-gammaTB
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =0.0d0 !-gammaTB
               else if(jj==2) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =0.0D0 !gammaTB !0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =0.0D0 !gammaTB ! 0.0d0
               else if(jj==3) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =0.0D0 !gammaTB ! 0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =0.0D0 !gammaTB ! 0.0d0
                  jj=0
               end if
            else if (Iend_TB(2, ii)>ns+np) then
               kk=kk+1
               if (kk==1) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  = gammaTB
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) = gammaTB
               else if (kk==2) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =gammaTB ! 0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==3) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  = gammaTB
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) = gammaTB
               else if (kk==4) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =gammaTB ! 0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==5) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  =gammaTB ! 0.0d0
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==6) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  = gammaTB
                  fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) = gammaTB
                  kk=0
               endif
            else
               fock_DFTB(Iend_TB(2,ii)+MTB,MTB-l1)  = gammaTB
               fock_DFTB(MTB-l1, Iend_TB(2,ii)+MTB) = gammaTB

            end if
!            l1=l1+1
!         else

!            fock_DFTB(Iend_TB(2,ii)+MTB,MTB) = 0.0D0
!            fock_DFTB(MTB, Iend_TB(2,ii)+MTB) = 0.0D0
!         end if


      else if (Iend_TB(1, ii) == natom) then
            if (Iend_TB(2, ii)>ns.and.Iend_TB(2, ii)<=ns+np) then
               jj=jj+1
               if (jj==1) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  =0.0d0 !gammaTB
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) =0.0d0 !gammaTB
               else if(jj==2) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  =0.0D0 !gammaTB ! 0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) =0.0D0 !gammaTB ! 0.0d0
               else if(jj==3) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  =0.0D0 !gammaTB ! 0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) =0.0D0 !gammaTB ! 0.0d0
                  jj=0
               end if

            else if (Iend_TB(2, ii)>ns+np) then
               kk=kk+1
               if (kk==1) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB
               else if (kk==2) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  =gammaTB ! 0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==3) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB
               else if (kk==4) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  =gammaTB ! 0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==5) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  =gammaTB ! 0.0d0
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==6) then
                  fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB
                  kk=0
               endif
            else
               fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
               fock_DFTB(MTB+M_in+1+l2, Iend_TB(2,ii)+MTB) = gammaTB

            end if
!            l2=l2+1
!         if ((ncont(Iend_TB(2,ii))==1).or.(ncont(Iend_TB(2,ii))==2)) then
!            fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1)  = gammaTB
!            fock_DFTB(MTB+M_in+1,Iend_TB(2,ii)+MTB) = gammaTB
!         else
!            fock_DFTB(Iend_TB(2,ii)+MTB,MTB+M_in+1)  = gammaTB
!            fock_DFTB(MTB+M_in+1, Iend_TB(2,ii)+MTB) = gammaTB
!         end if

      end if

   end do


   do ii=1,MTB
      fock_DFTB(ii,ii) = alfaTB-(Vbias_TB/2.0d0)*f_t
      fock_DFTB(MTB+M_in+ii, MTB+M_in+ii)= alfaTB+(Vbias_TB/2.0d0)*f_t

      if (ii<MTB) then ! -end_bTB) then
         fock_DFTB(ii,ii+1) = betaTB
         fock_DFTB(ii+1,ii) = betaTB
         fock_DFTB(2*MTB+M_in-ii, 2*MTB+M_in-ii+1)= betaTB
         fock_DFTB(2*MTB+M_in-ii+1, 2*MTB+M_in-ii)= betaTB

      end if
   end do

!   do ii=1, end_bTB

!      fock_DFTB(MTB-end_bTB, MTB-end_bTB+ii)=betaTB
!      fock_DFTB(MTB-end_bTB+ii, MTB-end_bTB)=betaTB
!      fock_DFTB(MTB+M_in+end_bTB+1, MTB+M_in+end_bTB+1-ii)=betaTB
!      fock_DFTB(MTB+M_in+end_bTB+1-ii,MTB+M_in+end_bTB+1)=betaTB

!   end do

   fock_DFTB(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)


end subroutine chimeraDFTB_evol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TB_current (M_in,rhold,rhonew, overlap, TB_A, TB_B, TB_M)
   use dftb_data, only:MDFTB, MTB

   implicit none
   integer, intent(in)  :: M_in
   real*8, intent(in)   :: overlap(M_in, M_in)
#ifdef TD_SIMPLE
   complex*8, intent(in)   :: rhold(MDFTB,MDFTB)
   complex*8, intent(in)   :: rhonew(MDFTB,MDFTB)
#else
   complex*16, intent(in)   :: rhold(MDFTB,MDFTB)
   complex*16, intent(in)   :: rhonew(MDFTB,MDFTB)
#endif
   real*8, intent(out)  :: TB_A, TB_B, TB_M
   integer              :: ii, jj
   real*8, allocatable  :: delta_rho(:,:)
   real*8               :: qe

   allocate(delta_rho(MDFTB,MDFTB))

   delta_rho=real(rhonew)-real(rhold)


   TB_A=0.0D0
   TB_B=0.0D0
   TB_M=0.0D0

   do ii=1, MTB

      TB_A=delta_rho(ii,ii) + TB_A
      TB_B=delta_rho(MTB+M_in+ii, MTB+M_in+ii) + TB_B

   end do

   do ii=1,M_in
   do jj=1,M_in

      qe = delta_rho(MTB+ii,MTB+jj) * overlap(ii, jj)
      TB_M = qe + TB_M

   enddo
   enddo

end subroutine TB_current

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine write_rhoDFTB(M_in, OPEN)
   use dftb_data , only: rho_aDFTB, rho_bDFTB

   implicit none
   integer, intent(in) :: M_in
   logical, intent(in) :: OPEN
   integer :: ii, jj

   open(unit=1070,file='rhoTB.out')

   if (OPEN) then
      do ii=1,M_in
      do jj=1,M_in
         write(1070,*) rho_aDFTB(ii,jj), rho_bDFTB(ii,jj)
      end do
      end do
   else
      do ii=1,M_in
      do jj=1,M_in
         write(1070,*) rho_aDFTB(ii,jj)
      end do
      end do
   end if

   close(1070)

   write(*,*) 'RHOTB wrtted'

end subroutine write_rhoDFTB


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dftb_output(M, dim3, rho_aux, overlap, istep, Iz, natom, Nuc, OPEN)
   use dftb_data, only: rhold_AOTB, rhonew_AOTB, MTB, MDFTB
   implicit none

   logical, intent(in) :: OPEN
   integer, intent(in) :: M, istep, dim3
   integer, intent(in) :: natom
   integer, intent(in) :: Nuc(M)
   integer, intent(in) :: Iz(natom)
   real*8, intent(in)  :: overlap(M,M)
#ifdef TD_SIMPLE
   complex*8, intent(in)  :: rho_aux(MDFTB, MDFTB,dim3)
#else
   complex*16, intent(in) :: rho_aux(MDFTB, MDFTB,dim3)
#endif
   integer  :: n, ii, jj
   real*8   :: I_TB_A(dim3), I_TB_B(dim3), I_TB_M(dim3)
   real*8   :: chargeA_TB, chargeB_TB, chargeM_TB
   real*8   :: orb_charge, tot_orb_charge
   real*8   :: qe(natom)
   real*8   :: rhoscratch(M,M)

   if (istep==1) then

      open(unit=10101,file='currentTB')
      open(unit=20202,file='mullikenTB')

   else

      if (OPEN) then
         call TB_current(M,rhold_AOTB(:,:,1),rhonew_AOTB(:,:,1), overlap,      &
                         I_TB_A(1), I_TB_B(1), I_TB_M(1))
         call TB_current(M,rhold_AOTB(:,:,2),rhonew_AOTB(:,:,2), overlap,      &
                         I_TB_A(2), I_TB_B(2), I_TB_M(2))

         write(10101,*) "A", I_TB_A(1) + I_TB_A(2)
         write(10101,*) "B", I_TB_B(1) + I_TB_B(2)
         write(10101,*) "M", I_TB_M(1) + I_TB_M(2)
      else
         call TB_current(M,rhold_AOTB(:,:,1),rhonew_AOTB(:,:,1), overlap,      &
                         I_TB_A(1), I_TB_B(1), I_TB_M(1))

         write(10101,*) "A", I_TB_A(1)
         write(10101,*) "B", I_TB_B(1)
         write(10101,*) "M", I_TB_M(1)
      end if


      chargeA_TB=MTB
      chargeB_TB=MTB
      do ii=1, MTB
         chargeA_TB=chargeA_TB-rho_aux(ii,ii,1)
         chargeB_TB=chargeB_TB-rho_aux(MTB+M+ii,MTB+M+ii,1)
      end do

      if (OPEN) then
         do ii=1, MTB
            chargeA_TB=chargeA_TB-rho_aux(ii,ii,2)
            chargeB_TB=chargeB_TB-rho_aux(MTB+M+ii,MTB+M+ii,2)
         end do
      end if

      chargeM_TB=0.0D0
      do n=1,natom
         qe(n)=Iz(n)
      enddo

      rhoscratch=real(rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,1))

      call mulliken_calc(natom, M, rhoscratch, overlap, Nuc, qe)

      do n=1,natom
         chargeM_TB= chargeM_TB - qe(n)
      enddo

      if(OPEN) then
         rhoscratch=real(rho_aux(MTB+1:MTB+M,MTB+1:MTB+M,2))

         call mulliken_calc(natom, M, rhoscratch, overlap, Nuc, qe)

         do n=1,natom
            chargeM_TB= chargeM_TB - qe(n)
         enddo
      end if

      write(20202,*) "A", chargeA_TB
      write(20202,*) "B", chargeB_TB
      write(20202,*) "M", chargeM_TB

   endif


end subroutine dftb_output

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module dftb_subs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
