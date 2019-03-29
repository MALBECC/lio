!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module tbdft_subs
   implicit none

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine tbdft_init(M_in, Nuc, natom, OPEN)

   use tbdft_data, only: MTB, MTBDFT, end_bTB, Iend_TB, rhoa_TBDFT, rhob_TBDFT,    &
                        gammaW

   implicit none

   logical, intent(in)  :: OPEN
   integer, intent(in)  ::  M_in, natom
   integer, intent(in)  ::  Nuc(M_in)
   integer   ::  ii, jj

   MTBDFT=2*MTB+M_in
   allocate(Iend_TB(2,2*end_bTB), rhoa_TBDFT(MTBDFT,MTBDFT),gammaW(2*end_bTB))
   if (OPEN) allocate (rhob_TBDFT(MTBDFT,MTBDFT))

   open(unit=1001,file='gamma.in')

   do ii=1, 2*end_bTB
      read(1001,*) gammaW(ii)
   enddo

   close(1001)

   jj=0
   do ii=1, M_in
      if (Nuc(ii)==1.or.Nuc(ii)==natom) then
         jj=jj+1
         Iend_TB(1,jj) = Nuc(ii)
         Iend_TB(2,jj) = ii
      end if
   end do

end subroutine tbdft_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine tbdft_td_init (M_in,rho, rho_0, dim3)

   use tbdft_data, only: MTB, MTBDFT, rhoa_TBDFT, rhob_TBDFT, rhold_AOTB,          &
                        rhonew_AOTB
   implicit none
   integer, intent(in)     :: M_in
   integer, intent(in)     :: dim3
#ifdef TD_SIMPLE
   complex*8, intent(in)   :: rho_0(M_in,M_in,dim3)
   complex*8, intent(out)  :: rho(MTBDFT,MTBDFT,dim3)
#else
   complex*16, intent(in)  :: rho_0(M_in,M_in,dim3)
   complex*16, intent(out) :: rho(MTBDFT,MTBDFT,dim3)
#endif
   integer  :: ii,jj

   allocate(rhold_AOTB(MTBDFT,MTBDFT,dim3),                       &
            rhonew_AOTB(MTBDFT,MTBDFT,dim3))

   rhold_AOTB=0.0d0
   rhonew_AOTB=0.0d0

      do jj=1, MTBDFT
      do ii=1, MTBDFT
         if (ii==jj) then
            rho(ii,jj,1)=cmplx(rhoa_TBDFT(ii,jj), 0.0D0)
         else
            rho(ii,jj,1)=cmplx(rhoa_TBDFT(ii,jj)/2.0d0, 0.0D0)
         end if
      end do
      end do

!carlos:Open shell option
   if(dim3==2) then
      do jj=1, MTBDFT
      do ii=1, MTBDFT
         if (ii==jj) then
            rho(ii,jj,2)=cmplx(rhob_TBDFT(ii,jj), 0.0D0)
         else
            rho(ii,jj,2)=cmplx(rhob_TBDFT(ii,jj)/2.0d0, 0.0D0)
         end if
      end do
      end do
   end if

   do jj= 1, M_in
   do ii= 1, M_in
      rho(ii+MTB,jj+MTB,:)=rho_0(ii,jj,:)
   end do
   end do

end subroutine tbdft_td_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%i%%%%%%%%%%%%%%%%%%%%%%%%

subroutine getXY_TBDFT(M_in,x_in,y_in,xmat,ymat)

   use tbdft_data, only: MTB, MTBDFT

   implicit none
   integer, intent(in)  :: M_in
   real*8 , intent(in)  :: x_in(M_in,M_in)
   real*8 , intent(in)  :: y_in(M_in,M_in)
   real*8 , intent(out) :: xmat(MTBDFT,MTBDFT)
   real*8 , intent(out) :: ymat(MTBDFT,MTBDFT)
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
subroutine construct_rhoTBDFT(M, rho, rho_0 ,rho_TBDFT, niter, OPEN)

   use tbdft_data , only: MTB, MTBDFT

   implicit none
   logical, intent(in) :: OPEN
   integer, intent(in) :: M
   integer, intent(in) :: niter
   real*8, intent(in)  :: rho_0(M,M)
   real*8, intent(in)  :: rho_TBDFT(MTBDFT,MTBDFT)
   real*8, intent(out) :: rho(MTBDFT,MTBDFT)
   integer :: ii, jj
   real*8  :: ocup

   ocup=1.0d0
   if (OPEN) ocup=0.5d0

   if (niter/=1) then

      do ii=1,MTBDFT
      do jj=ii+1,MTBDFT
         rho(ii,jj)=rho_TBDFT(ii,jj)/2
         rho(jj,ii)=rho(ii,jj)
      end do
      end do

   else if(niter==1) then

      rho=0.0D0
      do ii=1, MTB
         rho(ii,ii)=ocup
         rho(MTB+M+ii,MTB+M+ii)=ocup
      end do
         rho(MTB+1:MTB+M, MTB+1:MTB+M)=rho_0
   end if

end subroutine construct_rhoTBDFT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine build_chimera_TBDFT (M_in,fock_in, fock_TBDFT, natom)

   use tbdft_data, only:MTBDFT, MTB, Iend_TB, end_bTB, alfaTB, betaTB, gammaTB,  &
                       Vbias_TB, gammaW

   integer, intent(in)  :: M_in
   integer, intent(in)  :: natom
   real*8, intent(in)   :: fock_in (M_in, M_in)
   real*8, intent(out)  :: fock_TBDFT (MTBDFT, MTBDFT)
   integer              :: ii, link

   fock_TBDFT(:,:) = 0.0D0

   do ii=1, 2*end_bTB
      if(Iend_TB(1,ii)==1) link = MTB
      if(Iend_TB(1,ii)==natom) link = MTB+M_in+1
      fock_TBDFT(Iend_TB(2,ii)+MTB,link)=gammaW(ii)*gammaTB
      fock_TBDFT(link,Iend_TB(2,ii)+MTB)=gammaW(ii)*gammaTB
   end do

   do ii=1,MTB
      fock_TBDFT(ii,ii) = alfaTB
      fock_TBDFT(MTB+M_in+ii, MTB+M_in+ii)= alfaTB

      if (ii<MTB) then

         fock_TBDFT(ii,ii+1) = betaTB
         fock_TBDFT(ii+1,ii) = betaTB
         fock_TBDFT(2*MTB+M_in-ii, 2*MTB+M_in-ii+1)= betaTB
         fock_TBDFT(2*MTB+M_in-ii+1, 2*MTB+M_in-ii)= betaTB

      end if
   end do

   fock_TBDFT(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)

end subroutine build_chimera_TBDFT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine extract_rhoDFT (M_in, rho_in, rho_out)

   use tbdft_data, only: MTBDFT, MTB

   implicit none
   integer, intent(in)  :: M_in
   real*8, intent(in)   :: rho_in(MTBDFT,MTBDFT)
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
subroutine chimeraTBDFT_evol(M_in,fock_in, fock_TBDFT, natom, istep)

   use tbdft_data, only:MTBDFT, MTB, Iend_TB, end_bTB, alfaTB, betaTB, gammaTB,   &
                       Vbias_TB, start_tdtb, end_tdtb, gammaW

   integer, intent(in)  :: M_in
   integer, intent(in)  :: natom
   integer, intent(in)  :: istep
   real*8, intent(in)   :: fock_in (M_in, M_in)
   real*8, intent(out)  :: fock_TBDFT (MTBDFT, MTBDFT) !temporal dimensions
   real*8               :: pi=4.0D0*atan(1.0D0)
   real*8               :: lambda, t_step, f_t
   integer              :: ii, link


   lambda=1.0d0/real(end_tdtb-start_tdtb)

   if (istep < start_tdtb) then
      f_t=0.0D0
   else if(istep >= start_tdtb .and. istep < end_tdtb) then
      t_step=real(istep-start_tdtb)
      f_t=(-Cos(pi*lambda*t_step)+1.0D0)/2.0D0
   else if(istep >= end_tdtb) then
      f_t=1.0D0
   end if

   fock_TBDFT(:,:) = 0.0D0

   do ii=1, 2*end_bTB
      if(Iend_TB(1,ii)==1) link = MTB
      if(Iend_TB(1,ii)==natom) link = MTB+M_in+1
      fock_TBDFT(Iend_TB(2,ii)+MTB,link)=gammaW(ii)*gammaTB
      fock_TBDFT(link,Iend_TB(2,ii)+MTB)=gammaW(ii)*gammaTB
   end do

   do ii=1,MTB
      fock_TBDFT(ii,ii) = alfaTB-(Vbias_TB/2.0d0)*f_t
      fock_TBDFT(MTB+M_in+ii, MTB+M_in+ii)= alfaTB+(Vbias_TB/2.0d0)*f_t

      if (ii<MTB) then
         fock_TBDFT(ii,ii+1) = betaTB
         fock_TBDFT(ii+1,ii) = betaTB
         fock_TBDFT(2*MTB+M_in-ii, 2*MTB+M_in-ii+1)= betaTB
         fock_TBDFT(2*MTB+M_in-ii+1, 2*MTB+M_in-ii)= betaTB

      end if
   end do


   fock_TBDFT(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)


end subroutine chimeraTBDFT_evol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TB_current (M_in,rhold,rhonew, overlap, TB_A, TB_B, TB_M)
   use tbdft_data, only:MTBDFT, MTB

   implicit none
   integer, intent(in)  :: M_in
   real*8, intent(in)   :: overlap(M_in, M_in)
#ifdef TD_SIMPLE
   complex*8, intent(in)   :: rhold(MTBDFT,MTBDFT)
   complex*8, intent(in)   :: rhonew(MTBDFT,MTBDFT)
#else
   complex*16, intent(in)   :: rhold(MTBDFT,MTBDFT)
   complex*16, intent(in)   :: rhonew(MTBDFT,MTBDFT)
#endif
   real*8, intent(out)  :: TB_A, TB_B, TB_M
   integer              :: ii, jj
   real*8, allocatable  :: delta_rho(:,:)
   real*8               :: qe

   allocate(delta_rho(MTBDFT,MTBDFT))

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
subroutine tbdft_scf_output(M_in, OPEN)
   use tbdft_data, only: rhoa_TBDFT, rhob_TBDFT, MTBDFT, MTB

   implicit none
   logical         , intent(in) :: OPEN
   integer         , intent(in) :: M_in
   double precision    :: rho_aux(MTBDFT, MTBDFT)
   double precision    :: chargeA_TB, chargeB_TB
   integer             :: ii


   if (OPEN) then
      rho_aux=rhoa_TBDFT+rhob_TBDFT
   else
      rho_aux=rhoa_TBDFT
   end if

   chargeA_TB=MTB
   chargeB_TB=MTB
   do ii=1, MTB
      chargeA_TB=chargeA_TB-rho_aux(ii,ii)
      chargeB_TB=chargeB_TB-rho_aux(MTB+M_in+ii,MTB+M_in+ii)
   end do

   open(unit=20202,file='mullikenTB')

   write(20202,*) "Mulliken TB  part A", chargeA_TB
   write(20202,*) "Mulliken TB  part B", chargeB_TB

   close(20202)
end subroutine tbdft_scf_output

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine tbdft_td_output(M_in, dim3, rho_aux, overlap, istep, Iz, natom, Nuc,&
                           OPEN)
   use tbdft_data, only: rhold_AOTB, rhonew_AOTB, MTB, MTBDFT
   implicit none

   logical, intent(in) :: OPEN
   integer, intent(in) :: M_in, istep, dim3
   integer, intent(in) :: natom
   integer, intent(in) :: Nuc(M_in)
   integer, intent(in) :: Iz(natom)
   real*8, intent(in)  :: overlap(M_in,M_in)
#ifdef TD_SIMPLE
   complex*8, intent(in)  :: rho_aux(MTBDFT, MTBDFT,dim3)
#else
   complex*16, intent(in) :: rho_aux(MTBDFT, MTBDFT,dim3)
#endif
   integer  :: n, ii, jj
   real*8   :: I_TB_A(dim3), I_TB_B(dim3), I_TB_M(dim3)
   real*8   :: chargeA_TB, chargeB_TB, chargeM_TB
   real*8   :: orb_charge, tot_orb_charge
   real*8   :: qe(natom)
   real*8   :: rhoscratch(M_in,M_in)

   if (istep==1) then

      open(unit=10101,file='currentTB')
      open(unit=20202,file='mullikenTB')

   else

      if (OPEN) then
         call TB_current(M_in,rhold_AOTB(:,:,1),rhonew_AOTB(:,:,1), overlap,   &
                         I_TB_A(1), I_TB_B(1), I_TB_M(1))
         call TB_current(M_in,rhold_AOTB(:,:,2),rhonew_AOTB(:,:,2), overlap,   &
                         I_TB_A(2), I_TB_B(2), I_TB_M(2))

         write(10101,*) "Current TB  part A", I_TB_A(1) + I_TB_A(2)
         write(10101,*) "Current TB  part B", I_TB_B(1) + I_TB_B(2)
         write(10101,*) "Current DFT part M", I_TB_M(1) + I_TB_M(2)
      else
         call TB_current(M_in,rhold_AOTB(:,:,1),rhonew_AOTB(:,:,1), overlap,   &
                         I_TB_A(1), I_TB_B(1), I_TB_M(1))

         write(10101,*) "Current TB  part A", I_TB_A(1)
         write(10101,*) "Current TB  part B", I_TB_B(1)
         write(10101,*) "Current DFT part M", I_TB_M(1)
      end if


      chargeA_TB=MTB
      chargeB_TB=MTB
      do ii=1, MTB
         chargeA_TB=chargeA_TB-rho_aux(ii,ii,1)
         chargeB_TB=chargeB_TB-rho_aux(MTB+M_in+ii,MTB+M_in+ii,1)
      end do

      if (OPEN) then
         do ii=1, MTB
            chargeA_TB=chargeA_TB-rho_aux(ii,ii,2)
            chargeB_TB=chargeB_TB-rho_aux(MTB+M_in+ii,MTB+M_in+ii,2)
         end do
      end if

      chargeM_TB=0.0D0
      do n=1,natom
         qe(n)=Iz(n)
      enddo

      rhoscratch=real(rho_aux(MTB+1:MTB+M_in,MTB+1:MTB+M_in,1))

      call mulliken_calc(natom, M_in, rhoscratch, overlap, Nuc, qe)

      do n=1,natom
         chargeM_TB= chargeM_TB + qe(n)
      enddo

      if(OPEN) then
         rhoscratch=real(rho_aux(MTB+1:MTB+M_in,MTB+1:MTB+M_in,2))

         call mulliken_calc(natom, M_in, rhoscratch, overlap, Nuc, qe)

         do n=1,natom
            chargeM_TB= chargeM_TB + qe(n)
         enddo
      end if

      write(20202,*) "Mulliken TB   part A", chargeA_TB
      write(20202,*) "Mulliken TB   part B", chargeB_TB
      write(20202,*) "Mulliken DFT  part M", chargeM_TB

   endif


end subroutine tbdft_td_output
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module tbdft_subs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
