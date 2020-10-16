subroutine ESAget(Xexc,Eexc,Cscf,Ndim,NCO,Nvirt,M,nstates)
use garcha_mod,   only: Pmat_vec, nunp, iz, pc, r, d
use properties,   only: dipole
use excited_data, only: state_LR, read_other
   implicit none
   
   integer, intent(in) :: Ndim, NCO, Nvirt, M, nstates
   LIODBLE, intent(in) :: Xexc(Ndim,nstates), Eexc(nstates)
   LIODBLE, intent(in) :: Cscf(M,M)

   integer :: MM, ii, ef_stat, istat, jstat, io, jo, av, bv
   integer :: imo, jmo, amo, bmo, NCOc, ix, jx
   LIODBLE :: coef, EE
   LIODBLE, allocatable :: DipV(:,:), DipAO(:,:,:), DipMO(:,:,:), Xp(:), OsSt(:)
   LIODBLE, allocatable :: scratch(:,:), DipEx(:,:), tdd(:), gsdip(:), Ene(:)

   write(*,"(1X,A)") "EXCITED STATE ABSORPTION"
   write(*,"(1X,A,2X,I3)") "Using state in ESAfosc:", state_LR

   allocate(Xp(Ndim)); 
   ! Select which state is used in the perturbation
   if ( read_other ) then
      print*, "reading other state"
      open(unit=456,file="vectors_root")
      read(456,*) EE
      do ii=1,Ndim
          read(456,*) Xp(ii)
      enddo
      close(456)
      istat = 0
   else
      Xp = Xexc(:,state_LR)
      istat = state_LR
      EE = Eexc(state_LR)
   endif

   MM = M*(M+1)/2
   allocate(DipV(3,MM),DipAO(3,M,M),DipMO(3,M,M))
   DipV = 0.0d0; DipAO = 0.0d0; DipMO = 0.0d0

   ! Dipole moment integral matrix AO: in vector form
   call get_DipMatrix(DipV)

   ! Dipole moment integral matrix AO: in matrix form
   do ii=1,3
      call spunpack('L',M,DipV(ii,:),DipAO(ii,:,:))
   enddo
   deallocate(DipV)

   ! Dipole moment integral matrix MO
   allocate(scratch(M,M))
   do ii=1,3
      call dgemm('N','N',M,M,M,1.0d0,DipAO(ii,:,:),M,Cscf,M,0.0d0,scratch,M)
      call dgemm('T','N',M,M,M,1.0d0,Cscf,M,scratch,M,0.0d0,DipMO(ii,:,:),M)
   enddo
   deallocate(scratch,DipAO)

   ef_stat = nstates - istat
   allocate(tdd(3),DipEx(ef_stat,3),gsdip(3),Ene(ef_stat),OsSt(ef_stat))

   ! Ground State Dipole Moment: electronic + nuclear
   call dipole(gsdip, Pmat_vec, 2*NCO+Nunp, r, d, Iz, pc)
   gsdip = gsdip / 2.54d0 ! Debye to a.u.

   ! Dipole Moment between excited states
   NCOc = NCO + 1
   do jstat=1, ef_stat
      tdd = 0.0d0
      do io=1,NCO
      do av=1,Nvirt
         ix  = (io - 1) * Nvirt + av
         imo = NCOc - io
         amo = NCO  + av
         do jo=1,NCO
         do bv=1,Nvirt
            jx = (jo - 1) * Nvirt + bv
            jmo = NCOc - jo
            bmo = NCO  + bv
            coef = Xp(ix) * Xexc(jx,jstat+istat)

            if (imo==jmo.and.amo/=bmo) then
               tdd=tdd+coef*DipMO(:,amo,bmo)
            else if (imo/=jmo.and.amo==bmo) then
               tdd=tdd-coef*DipMO(:,imo,jmo)
            else if (imo==jmo.and.amo==bmo) then
               tdd=tdd+coef*(gsdip-DipMO(:,imo,imo)+DipMO(:,amo,amo))
            end if
         enddo ! jo
         enddo ! bv
      enddo ! io
      enddo ! av
      DipEx(jstat,:) = tdd 
      Ene(jstat) = Eexc(jstat+istat) - EE
   enddo
   deallocate(DipMO,tdd,gsdip,Xp)

   ! Oscillator Strenght Calculation and Print Results
   call PrintESA(DipEx,Ene,OsSt,ef_stat,istat)
   deallocate(DipEx,Ene,OsSt)
end subroutine ESAget
