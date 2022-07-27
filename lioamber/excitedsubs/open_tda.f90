subroutine open_linear_response(MatCoefA,VecEneA,MatCoefB,VecEneB,M,Mlr, &
                                NvirtA,NCOA,NdimA,NvirtB,NCOB,NdimB)
! This routine perform a Davidson diagonalization in order to obtain
! Excitation Energies and Transition vectors of excited states in open shell

! Inputs
! MatCoef: Molecular Orbitals Coefficients
! VecEne : Molecular Orbitals Energies

! Ouputs
! No outputs in open shell 
use garcha_mod  , only: npas
use excited_data, only: nstates, fittExcited, use_last, guessLR, max_subs
   implicit none
   integer, intent(in)  :: M, Mlr, NvirtA, NCOA, NdimA, NvirtB, NCOB, NdimB
   LIODBLE, intent(in)  :: MatCoefA(M,Mlr), VecEneA(Mlr), MatCoefB(M,Mlr), VecEneB(Mlr)

   character(len=8) :: char_max
   integer :: maxIter, iter, vec_dim, first_vec, newvec
   integer :: Subdim, Ndim_tot, Ndim_max
   LIODBLE, dimension(:,:), allocatable :: H,Hb,eigvec
   LIODBLE, dimension(:,:), allocatable :: AX_a, AX_b, tvecMOa, tvecMOb
   LIODBLE, dimension(:,:), allocatable :: RitzVecA, ResMatA
   LIODBLE, dimension(:,:), allocatable :: RitzVecB, ResMatB
   LIODBLE, dimension(:), allocatable :: eigval, val_old, Osc
   logical :: conv
   conv = .false. ! Convergence criteria bool

   print*, ""
   print*,"======================================="
   print*,"   LINEAR RESPONSE OPEN SHELL - TDA"
   print*,"======================================="
   print*, ""

   ! SET DIMENSION OF ALPHA AND BETA MATRICES
   Ndim_tot = NdimA + NdimB
   Ndim_max = max(NdimA,NdimB)

   ! Check Input Variables
   if ( nstates > Ndim_tot ) then
      print*, "Number of Excited States That You Want is Bigger than &
               & Dimension Matrix "
      nstates = Ndim_tot
      print*, "We will calculate", nstates
   endif

   ! Davidson Initialization
   allocate(val_old(nstates),Osc(nstates))
   ! TODO: 4 is the initial vectors for each excitated state ( INPUT ).
   val_old = 1.0d0; vec_dim = 4 * nstates
   first_vec = 0; newvec = 0
   if ( vec_dim >= Ndim_tot ) then
      vec_dim = Ndim_tot; Subdim = Ndim_tot
      maxIter = 1; max_subs = Ndim_tot
   else
      if ( max_subs > Ndim_tot ) max_subs = Ndim_tot
      maxIter = 50; Subdim = vec_dim
   endif
   allocate(tvecMOa(Ndim_max,max_subs),tvecMOb(Ndim_max,max_subs))
   allocate(AX_a(Ndim_max,max_subs),AX_b(Ndim_max,max_subs))

   ! This line is to avoid warnings.
   allocate(eigval(1))

   ! Initial Guess
   call open_vec_init(VecEneA,VecEneB,tvecMOa,tvecMOb,Ndim_max,&
                      vec_dim,Mlr,NCOA,NCOB)

   allocate(RitzVecA(Ndim_max,nstates),ResMatA(Ndim_max,nstates))
   allocate(RitzVecB(Ndim_max,nstates),ResMatB(Ndim_max,nstates))
   RitzVecA = 0.0d0; RitzVecB = 0.0d0

   ! Print Information
   write (char_max, '(i8)') max_subs
   write(*,"(1X,A,5X,I3,2X,I3,2X,I3,2X,I3,2X,I3,2X,I3)") "NCOA, NCOB, NVIRTA, NVIRTB, M, Mlr",NCOA,NCOB,NvirtA,NvirtB,M,Mlr
   write(*,"(1X,A,5X,I5,2X,I5,2X,I5)") "DIMENSION OF ALPHA, BETA AND TOTAL", NdimA, NdimB, Ndim_tot
   write(*,"(1X,A,23X,A)") "MAX SUBSPACE",adjustl(char_max)
   write(*,"(1X,A,20X,I2)") "MAX ITERATIONES",maxIter
   write(*,"(1X,A,3X,I4)") "NUMBER OF INITIAL TRIALS VECTORS",vec_dim

   ! DAVIDSON START
   do iter=1,maxIter
      call g2g_timer_start("Iteration LR")
      write(*,"(A)") " "
      write(*,"(1X,A,6X,I2)") "ITERATION:",iter
      write(*,"(1X,A,7X,I4)") "SUBSPACE:",Subdim
      write(*,"(1X,A,1X,I4)") "VECTORS INSIDE:",vec_dim

      ! This routine calculate Fock and form The Excited Matrix.
      call open_solve_focks(MatCoefA,tvecMOa,MatCoefB,tvecMOb,AX_a,AX_b, &
                            M,Mlr,NCOA,NvirtA,NCOB,NvirtB,Ndim_max,max_subs, &
                            vec_dim,Subdim,first_vec)
      ! AX_a,b += (Ea-Ei)*Xia 
      call addInt(AX_a,VecEneA,tvecMOa,Ndim_max,Mlr,Subdim,NCOA,vec_dim,first_vec)
      call addInt(AX_b,VecEneB,tvecMOb,Ndim_max,Mlr,Subdim,NCOB,vec_dim,first_vec)

      ! We obtain subspace matrix
      allocate(H(Subdim,Subdim),Hb(Subdim,Subdim))
      call dgemm('T','N',Subdim,Subdim,Ndim_max,1.0d0,tvecMOa,Ndim_max,AX_a,&
                  Ndim_max,0.0d0,H,Subdim)
      call dgemm('T','N',Subdim,Subdim,Ndim_max,1.0d0,tvecMOb,Ndim_max,AX_b,&
                  Ndim_max,0.0d0,Hb,Subdim)
      H = H + Hb ! alpha + beta

      ! Diagonalization of Subspace Matrix
      if(allocated(eigvec)) deallocate(eigvec)
      if(allocated(eigval)) deallocate(eigval)
      allocate(eigvec(Subdim,Subdim),eigval(Subdim))
      call diagonH(H,Subdim,eigval,eigvec)
      deallocate(H,Hb)

      ! Change subspace to full space - Obtain Ritz Vector
      call dgemm('N','N',Ndim_max,nstates,Subdim,1.0d0,tvecMOa,Ndim_max,&
                 eigvec,Subdim,0.0d0,RitzVecA,Ndim_max)
      call dgemm('N','N',Ndim_max,nstates,Subdim,1.0d0,tvecMOb,Ndim_max,&
                 eigvec,Subdim,0.0d0,RitzVecB,Ndim_max)

      ! When do not perform davidson diagonalization
      if ( maxIter == 1 ) then
         call open_OscStr(RitzVecA,RitzVecB,eigval,MatCoefA,MatCoefB, &
                          Osc,M,Mlr,NCOA,NCOB,NvirtA,NvirtB,Ndim_max,nstates)

         call open_PrintResults(RitzVecA,RitzVecB,eigval,Osc,Ndim_max,nstates,Mlr,NCOa,NCOb)
         exit
      endif

      ! Obtain residual vectors
      call residual(eigval,eigvec,RitzVecA,AX_a,ResMatA,Ndim_max,Subdim,nstates)
      call residual(eigval,eigvec,RitzVecB,AX_b,ResMatB,Ndim_max,Subdim,nstates)

      ! Check Convergence and append new vectors
      conv = .false.; newvec = 0
      call open_new_vectors(ResMatA,ResMatB,eigval,VecEneA,VecEneB,tvecMOa,tvecMOb, &
                       val_old,Ndim_max,Subdim,nstates,Mlr,NCOA,NCOB,NdimA,NdimB,newvec,conv)
      if ( conv .eqv. .true. ) then
         write(*,*) ""
         write(*,"(1X,A,I2,1X,A)") "CONVERGED IN:",iter,"ITERATIONS"
         call open_OscStr(RitzVecA,RitzVecB,eigval,MatCoefA,MatCoefB, &
                          Osc,M,Mlr,NCOA,NCOB,NvirtA,NvirtB,Ndim_max,nstates)

         call open_PrintResults(RitzVecA,RitzVecB,eigval,Osc,Ndim_max,nstates,Mlr,NCOa,NCOb)
         exit
      else
         ! Actualization of vectors index
         first_vec = Subdim
         Subdim = Subdim + newvec
         vec_dim = newvec
      endif
      call g2g_timer_stop("Iteration LR")
   enddo ! END DAVIDSON ITERATION

   ! Free Memory
   deallocate(RitzVecA,RitzVecB,eigval,eigvec,ResMatA,ResMatB,AX_a,AX_b,tvecMOa,tvecMOb)
   deallocate(val_old,Osc)

end subroutine open_linear_response

subroutine open_vec_init(EneA,EneB,VecA,VecB,Ndim,vecnum,M,NCOA,NCOB)
implicit none

   integer, intent(in)   :: Ndim, vecnum, M, NCOA, NCOB
   LIODBLE, intent(in)   :: EneA(M), EneB(M)
   LIODBLE, intent(out)  :: VecA(Ndim,vecnum), VecB(Ndim,vecnum)

   integer :: ii, occ, virt, cont, start, NdimA, NdimB, NdimT
   integer :: NvirtA, NvirtB
   integer, dimension(:), allocatable :: ind
   LIODBLE, dimension(:), allocatable :: deltaE

   ! Set indexes again
   NvirtA = M-NCOA
   NvirtB = M-NCOB
   NdimA = NCOA * NvirtA
   NdimB = NCOB * NvirtB
   NdimT = NdimA + NdimB

   allocate(deltaE(NdimT),ind(NdimT))
   ! Calculate delta molecular orbital energies for alpha
   do ii=1,NdimA
      cont = ii - 1
      occ = NCOA - (cont/NvirtA)
      virt = mod(cont,NvirtA) + NCOA + 1
      deltaE(ii) = EneA(virt) - EneA(occ)
   enddo

   ! Calculate delta molecular orbital energies for beta
   do ii=1,NdimB
      cont = ii - 1
      occ = NCOB - (cont/NvirtB)
      virt = mod(cont,NvirtB) + NCOB + 1
      deltaE(NdimA+ii) = EneB(virt) - EneB(occ)
   enddo

   ! Sorting Energies
   ind = 0
   call eigsort(deltaE,ind,NdimT)

!  Initial trial vectors
   VecA = 0.0d0; VecB = 0.0d0
   do ii=1,vecnum
      if (ind(ii) <= NdimA) then
         VecA(ind(ii),ii) = 1.0d0
      else
         VecB(ind(ii)-NdimA,ii) = 1.0d0
      endif
   enddo
   deallocate(ind,deltaE)
end subroutine open_vec_init

