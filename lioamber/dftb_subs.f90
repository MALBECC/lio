!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dftb_subs
   implicit none

contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dftb_init(M)

   use dftb_data, only: MTB, MDFTB, end_bTB, Iend_TB, rho_DFTB

   implicit none

   integer, intent(in) :: M

   MDFTB=2*MTB+M
   allocate(Iend_TB(2,2*end_bTB), rho_DFTB(MDFTB,MDFTB))

end subroutine dftb_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dftb_td_init (M,rho, rho_0, overlap, RMM5)

   use dftb_data, only: MTB, MDFTB, rho_DFTB, rhold_AOTB, rhonew_AOTB
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
   allocate(overlap(M,M), rhold_AOTB(MDFTB,MDFTB), rhonew_AOTB(MDFTB,MDFTB))

   rhold_AOTB=0.0d0
   rhonew_AOTB=0.0d0

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

subroutine dftb_output(M, rho_aux, overlap, istep, Iz, natom, Nuc)
   use dftb_data, only: rhold_AOTB, rhonew_AOTB, MTB, MDFTB
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

      call TB_current(M,rhold_AOTB,rhonew_AOTB, overlap, I_TB_A, I_TB_B, I_TB_M)

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
