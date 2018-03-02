!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dftb_subs
   implicit none

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dftb_init(M)

   use dftb_data, only: MTB, MDFTB, end_basis, Iend, rho_DFTB

   implicit none

   integer, intent(in) :: M

   MDFTB=2*MTB+M
   allocate(Iend(2,2*end_basis), rho_DFTB(MDFTB,MDFTB))
   
end subroutine dftb_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dftb_td_init (M,rho, rho_0, overlap, RMM5)

   use dftb_data, only: MTB, MDFTB, rho_DFTB, rhold_AO, rhonew_AO
   implicit none
   integer, intent(in)       :: M
   real*8, allocatable, intent(inout) :: overlap(:,:)
   real*8 , intent(in)  :: RMM5(M*(M+1)/2)
#ifdef TD_SIMPLE
   complex*8, intent(in)     :: rho_0(M,M)
   complex*8, intent(out)  :: rho(MDFTB,MDFTB)
#else
   complex*16, intent(in)    :: rho_0(M,M)
   complex*16, intent(out) :: rho(MDFTB,MDFTB)
#endif
   integer  :: ii,jj

   if (allocated(overlap)) deallocate(overlap)
   allocate(overlap(M,M), rhold_AO(MDFTB,MDFTB), rhonew_AO(MDFTB,MDFTB))
   
   rhold_AO=0.0d0
   rhonew_AO=0.0d0

   call spunpack('L', M, RMM5, overlap)

      do jj=1, MDFTB
      do ii=1, MDFTB
         if (ii==jj) then
            rho(ii,jj)=cmplx(rho_DFTB(ii,jj), 0.0D0)
         else
            rho(ii,jj)=cmplx(rho_DFTB(ii,jj)/2.0d0, 0.0D0)
         end if
      end do
      end do

      do jj= 1, M
      do ii= 1, M
         rho(ii+MTB,jj+MTB)=rho_0(ii,jj)
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
subroutine find_neighbors(M_in,Nuc,natom)

   use dftb_data, only:Iend
  
   implicit none
   
   integer              ::  ii, jj
   integer, intent(in)  ::  M_in, natom
   integer, intent(in)  ::  Nuc(M_in)

   
   jj=0
   do ii=1, M_in
      if (Nuc(ii)==1.or.Nuc(ii)==natom) then
         jj=jj+1
         Iend(1,jj) = Nuc(ii)
         Iend(2,jj) = ii
      end if
   end do   

end subroutine find_neighbors

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine build_chimera (M_in,fock_in, chimerafock, natom, nshell, ncont)

   use dftb_data, only:MDFTB, MTB, Iend, end_basis, alfaTB, betaTB, gammaTB,Vbias

   integer, intent(in)  :: M_in
   integer, intent(in)  :: natom
   integer, intent(in)  :: nshell (0:4)
   integer, intent(in)  :: ncont(M_in)
   real*8, intent(in)   :: fock_in (M_in, M_in)
   real*8, intent(out)  :: chimerafock (MDFTB, MDFTB)
   integer              :: ii, jj, kk, ns, np, l1, l2

   l1=0
   l2=0
   jj=0
   kk=0
   ns=nshell(0)
   np=nshell(1)

   chimerafock(:,:) = 0.0D0

   do ii=1, 2*end_basis

      if (Iend(1, ii) == 1) then

!         if ((ncont(Iend(2,ii))==1).or.(ncont(Iend(2,ii))==2)) then

            if (Iend(2, ii)>ns.and.Iend(2, ii)<=ns+np) then
               jj=jj+1
               if (jj==1) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =0.0d0 !-gammaTB
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =0.0d0 !-gammaTB
               else if(jj==2) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  = 0.0D0 !gammaTB !0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) = 0.0D0 !gammaTB !0.0d0
               else if(jj==3) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  = 0.0D0 !gammaTB !0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) = 0.0D0 !gammaTB !0.0d0
                  jj=0
               end if

            else if (Iend(2, ii)>ns+np) then
               kk=kk+1
               if (kk==1) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =gammaTB
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =gammaTB
               else if (kk==2) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =gammaTB !0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =gammaTB !0.0d0
               else if (kk==3) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =gammaTB
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =gammaTB
               else if (kk==4) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =gammaTB !0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =gammaTB !0.0d0
               else if (kk==5) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =gammaTB !0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =gammaTB !0.0d0
               else if (kk==6) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =gammaTB
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =gammaTB
                  kk=0
               endif
            else
               chimerafock(Iend(2,ii)+MTB,MTB-l1)  = gammaTB
               chimerafock(MTB-l1, Iend(2,ii)+MTB) = gammaTB

            end if
!            l1=l1+1

!         else

!            chimerafock(Iend(2,ii)+MTB,MTB) = 0.0D0
!            chimerafock(MTB, Iend(2,ii)+MTB) = 0.0D0

!         end if         

      else if (Iend(1, ii) == natom) then
            if (Iend(2, ii)>ns.and.Iend(2, ii)<=ns+np) then
               jj=jj+1
               if (jj==1) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  =0.0d0 !gammaTB
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) =0.0d0 !gammaTB
               else if(jj==2) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  =0.0D0 !gammaTB ! 0.0d0
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) =0.0D0 !gammaTB !0.0d0
               else if(jj==3) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  =0.0D0 !gammaTB !0.0d0
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) =0.0D0 !gammaTB !0.0d0
                  jj=0
               end if 
            else if (Iend(2, ii)>ns+np) then
               kk=kk+1
               if (kk==1) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB
               else if (kk==2) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB !0.0d0
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB !0.0d0
               else if (kk==3) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB
               else if (kk==4) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB !0.0d0
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB !0.0d0 
               else if (kk==5) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB !0.0d0 
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB !0.0d0 
               else if (kk==6) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB
                  kk=0
               endif
            else
               chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
               chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB

            end if
!            l2=l2+1
!         if ((ncont(Iend(2,ii))==1) .or. (ncont(Iend(2,ii))==2)) then

!            chimerafock(Iend(2,ii)+MTB,MTB+M_in+1)  = gammaTB
!            chimerafock(MTB+M_in+1, Iend(2,ii)+MTB) = gammaTB
!         else
!            chimerafock(Iend(2,ii)+MTB,MTB+M_in+1)  = gammaTB
!            chimerafock(MTB+M_in+1, Iend(2,ii)+MTB) = gammaTB
!         end if

      end if

   end do

   do ii=1,MTB
      chimerafock(ii,ii) = alfaTB-Vbias/2.0d0
      chimerafock(MTB+M_in+ii, MTB+M_in+ii)= alfaTB+Vbias/2.0d0

      if (ii<MTB) then !-end_basis) then

         chimerafock(ii,ii+1) = betaTB
         chimerafock(ii+1,ii) = betaTB
         chimerafock(2*MTB+M_in-ii, 2*MTB+M_in-ii+1)= betaTB
         chimerafock(2*MTB+M_in-ii+1, 2*MTB+M_in-ii)= betaTB

      end if
   end do

!   do ii=1, end_basis

!      chimerafock(MTB-end_basis, MTB-end_basis+ii)=betaTB
!      chimerafock(MTB-end_basis+ii, MTB-end_basis)=betaTB
!      chimerafock(MTB+M_in+end_basis+1, MTB+M_in+end_basis+1-ii)=betaTB
!      chimerafock(MTB+M_in+end_basis+1-ii,MTB+M_in+end_basis+1)=betaTB

!   end do

   chimerafock(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)


end subroutine build_chimera

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
subroutine chimera_evol(M_in,fock_in, chimerafock, natom, nshell,ncont, istep)

   use dftb_data, only:MDFTB, MTB, Iend, end_basis, alfaTB, betaTB, gammaTB,   &
                       Vbias, start_tdtb, end_tdtb

   integer, intent(in)  :: M_in
   integer, intent(in)  :: natom
   integer, intent(in)  :: ncont(M_in)
   integer, intent(in)  :: istep
   integer, intent(in)  :: nshell (0:4)
   real*8, intent(in)   :: fock_in (M_in, M_in)
   real*8, intent(out)  :: chimerafock (MDFTB, MDFTB) !temporal dimensions
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

   chimerafock(:,:) = 0.0D0

   do ii=1, 2*end_basis
      if (Iend(1, ii) == 1) then
!         if ((ncont(Iend(2,ii))==1) .or. (ncont(Iend(2,ii))==2)) then

            if (Iend(2, ii)>ns.and.Iend(2, ii)<=ns+np) then
               jj=jj+1
               if (jj==1) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =0.0d0  !-gammaTB
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =0.0d0 !-gammaTB
               else if(jj==2) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =0.0D0 !gammaTB !0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =0.0D0 !gammaTB ! 0.0d0
               else if(jj==3) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =0.0D0 !gammaTB ! 0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =0.0D0 !gammaTB ! 0.0d0
                  jj=0
               end if
            else if (Iend(2, ii)>ns+np) then
               kk=kk+1
               if (kk==1) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  = gammaTB
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) = gammaTB
               else if (kk==2) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =gammaTB ! 0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==3) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  = gammaTB
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) = gammaTB
               else if (kk==4) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =gammaTB ! 0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==5) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  =gammaTB ! 0.0d0
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==6) then
                  chimerafock(Iend(2,ii)+MTB,MTB-l1)  = gammaTB
                  chimerafock(MTB-l1, Iend(2,ii)+MTB) = gammaTB
                  kk=0
               endif
            else
               chimerafock(Iend(2,ii)+MTB,MTB-l1)  = gammaTB
               chimerafock(MTB-l1, Iend(2,ii)+MTB) = gammaTB

            end if
!            l1=l1+1
!         else

!            chimerafock(Iend(2,ii)+MTB,MTB) = 0.0D0
!            chimerafock(MTB, Iend(2,ii)+MTB) = 0.0D0
!         end if


      else if (Iend(1, ii) == natom) then
            if (Iend(2, ii)>ns.and.Iend(2, ii)<=ns+np) then
               jj=jj+1
               if (jj==1) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  =0.0d0 !gammaTB
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) =0.0d0 !gammaTB
               else if(jj==2) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  =0.0D0 !gammaTB ! 0.0d0
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) =0.0D0 !gammaTB ! 0.0d0
               else if(jj==3) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  =0.0D0 !gammaTB ! 0.0d0
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) =0.0D0 !gammaTB ! 0.0d0
                  jj=0
               end if

            else if (Iend(2, ii)>ns+np) then
               kk=kk+1
               if (kk==1) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB
               else if (kk==2) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  =gammaTB ! 0.0d0
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==3) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB
               else if (kk==4) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  =gammaTB ! 0.0d0
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==5) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  =gammaTB ! 0.0d0
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) =gammaTB ! 0.0d0
               else if (kk==6) then
                  chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
                  chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB
                  kk=0
               endif
            else
               chimerafock(Iend(2,ii)+MTB,MTB+M_in+1+l2)  = gammaTB
               chimerafock(MTB+M_in+1+l2, Iend(2,ii)+MTB) = gammaTB

            end if
!            l2=l2+1
!         if ((ncont(Iend(2,ii))==1).or.(ncont(Iend(2,ii))==2)) then
!            chimerafock(Iend(2,ii)+MTB,MTB+M_in+1)  = gammaTB
!            chimerafock(MTB+M_in+1,Iend(2,ii)+MTB) = gammaTB
!         else
!            chimerafock(Iend(2,ii)+MTB,MTB+M_in+1)  = gammaTB
!            chimerafock(MTB+M_in+1, Iend(2,ii)+MTB) = gammaTB
!         end if

      end if

   end do


   do ii=1,MTB
      chimerafock(ii,ii) = alfaTB-(Vbias/2.0d0)*f_t
      chimerafock(MTB+M_in+ii, MTB+M_in+ii)= alfaTB+(Vbias/2.0d0)*f_t

      if (ii<MTB) then ! -end_basis) then
         chimerafock(ii,ii+1) = betaTB
         chimerafock(ii+1,ii) = betaTB
         chimerafock(2*MTB+M_in-ii, 2*MTB+M_in-ii+1)= betaTB
         chimerafock(2*MTB+M_in-ii+1, 2*MTB+M_in-ii)= betaTB

      end if
   end do

!   do ii=1, end_basis

!      chimerafock(MTB-end_basis, MTB-end_basis+ii)=betaTB
!      chimerafock(MTB-end_basis+ii, MTB-end_basis)=betaTB
!      chimerafock(MTB+M_in+end_basis+1, MTB+M_in+end_basis+1-ii)=betaTB
!      chimerafock(MTB+M_in+end_basis+1-ii,MTB+M_in+end_basis+1)=betaTB

!   end do

   chimerafock(MTB+1:MTB+M_in, MTB+1:MTB+M_in) = fock_in(:,:)


end subroutine chimera_evol

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

subroutine dftb_output(M, rho_aux, overlap, istep, Iz, natom, Nuc)
   use dftb_data, only: rhold_AO, rhonew_AO, MTB, MDFTB
   implicit none
   integer, intent(in) :: M, istep
   integer, intent(in) :: natom
   integer, intent(in) :: Nuc(M)
   integer, intent(in) :: Iz(natom)
   real*8, intent(in)  :: overlap(M,M)
#ifdef TD_SIMPLE
   complex*8, intent(in)  :: rho_aux(MDFTB, MDFTB)
#else
   complex*16, intent(in) :: rho_aux(MDFTB, MDFTB)
#endif
   integer  :: n, ii, jj
   real*8   :: I_TB_A, I_TB_B, I_TB_M
   real*8   :: chargeA_TB, chargeB_TB, chargeM_TB
   real*8   :: orb_charge, tot_orb_charge
   real*8   :: qe(natom)
   real*8   :: rhoscratch(M,M) 

   if (istep==1) then

      open(unit=10101,file='currentTB')
      open(unit=20202,file='mullikenTB')

   else

      call TB_current(M,rhold_AO,rhonew_AO, overlap, I_TB_A, I_TB_B, I_TB_M)

      write(10101,*) "A", I_TB_A
      write(10101,*) "B", I_TB_B
      write(10101,*) "M", I_TB_M

      chargeA_TB=MTB
      chargeB_TB=MTB
      do ii=1, MTB
         chargeA_TB=chargeA_TB-rho_aux(ii,ii)
         chargeB_TB=chargeB_TB-rho_aux(MTB+M+ii,MTB+M+ii)
      end do

      chargeM_TB=0.0D0
      do n=1,natom
         qe(n)=Iz(n)
      enddo

      rhoscratch=real(rho_aux(MTB+1:MTB+M,MTB+1:MTB+M))

      call mulliken_calc(natom, M, rhoscratch, overlap, Nuc, qe)

      do n=1,natom
         chargeM_TB= chargeM_TB - qe(n)
      enddo

      write(20202,*) "A", chargeA_TB
      write(20202,*) "B", chargeB_TB
      write(20202,*) "M", chargeM_TB

   endif


end subroutine dftb_output

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module dftb_subs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
