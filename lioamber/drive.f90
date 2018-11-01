!------------------------------------------------------------------------------!
! READS AND DRIVES SUBROUTINE                                                  !
! reads and prepare basis sets and geometry to perform electronic structure
! calculations using density functional theory and gaussian basis sets for
! the expansion of MO, charge density and exchange correlation potential
!
! Started 13 March 1992, Dario Estrin
!
!------------------------------------------------------------------------------!

subroutine drive(ng2, ngDyn, ngdDyn, iostat)
      USE garcha_mod, ONLY : a,c, basis, done, done_fit, natomc, nnps,         &
      nnpp, nnpd, nns, nnp, nnd, atmin, jatc, ncf, lt, at, ct, nnat, nshell,   &
      nuc, ncont, nlb, nshelld, cd, ad, Nucd, ncontd, nld, Nucx, indexii,      &
      ncontx, cx, ax, indexiid, X, RMM, rhoalpha,rhobeta, af, charge,          &
      basis_set, fitting_set, e_, e_2, e3, exists, NORM, fcoord,               &
      fmulliken, natom, frestart, M, FAC, Iexch, int_basis, max_func, &
      frestartin, Md, NCO, nng, npas, Nr, used, STR, omit_bas, Nr2,   &
      wang, wang2, wang3, VCINP, OPEN, whatis, Num, Iz, pi,             &
      Rm2, rqm, rmax, Nunp, nl, nt, ng, ngd, restart_freq,             &
      writexyz, number_restr, restr_pairs,restr_index,restr_k,restr_w,restr_r0,&
      mulliken, MO_coef_at, MO_coef_at_b, use_libxc, ex_functional_id, &
      ec_functional_id, gpu_level

      USE ECP_mod, ONLY : ecpmode, asignacion
      USE fileio , ONLY : read_coef_restart, read_rho_restart
      use td_data    , only: td_do_pop
      use fileio_data, only: verbose, rst_dens
      use ghost_atoms_subs, only: summon_ghosts

      IMPLICIT NONE
      LOGICAL :: basis_check
      CHARACTER*255 :: int_basis_file, fit_basis_file
      CHARACTER*255 :: liohome
      CHARACTER*255 :: inp_line
      CHARACTER :: inp_char
      CHARACTER*3 :: simb

      INTEGER, INTENT(IN) :: ng2, ngDyn, ngdDyn
      INTEGER, INTENT(INOUT) :: iostat
      REAL*8 :: atmint, iprob
      REAL*8 :: xnorm
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: restart_coef, restart_coef_b, &
                                             restart_dens, restart_adens,  &
                                             restart_bdens
      INTEGER :: NCOa, NCOb, ncon, nraw
      INTEGER :: is,ip,id, index
      INTEGER :: igpu, ios, NBAS, iatom
      INTEGER :: M3, M5, M7, M9, M11, M13, M15, M17, M18, M18b, MM, MMd !punteros
      INTEGER :: Nd, Ndim, nopt
      INTEGER :: i,ikk,j,k,k1,kk,kl,kkk,l,l1,l2,n,NN, No !auxiliars
!------------------------------------------------------------------------------!
! parameters for 2 basis sets, normal, density
! ng maximum # of contracted functions , nl maximum # of primitives in
! a contraction
! c, cd ce , linear coefficients
! a , ad, ae exponents
! Nuc , indicates center , ncont # of primitives in the correspondent
! contraction, and nx ny nz exponents of x , y and z in the cartesian
! gaussian
! ------------------
!
! Dimensions for Overlap and Fock matrices
!
!
! Angular momenta : up to f functions ( could be easily extended if
! necessary)
!
! --- Defaults ---------------------------------------------------------
!
! NORM true , expansion in normalized gaussians, so normalization factor
! included in coefficients

      nopt=0
!
! calls generator of table for incomplete gamma functions
!
       call GENERF
       call GENERFS
       call GRIDLIO
       npas=0
!-----------------------------------------------------------------------

! reads input file
   if (writexyz) open(unit=18,file=fcoord)
   if ((mulliken).or.(td_do_pop.gt.0)) open(unit=85,file=fmulliken)
   if (restart_freq.gt.0) open(unit=88,file=frestart)
!c---------------------------------------------------------
!c POINTERS -----------------------------------------------
!c
!c first P
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2

!c now F alpha
      M3=MM+1
!c now S, F beta also uses the same position after S was used
      M5=MM+M3
!c now G
      M7=MM+M5
!c now Gm
      M9=M7+MMd
!c now H
      M11=M9+MMd
!c W ( eigenvalues ), also this space is used in least squares
      M13=MM+M11
!c aux ( vector for ESSl)
      M15=M+M13
!c Least squares
      M17=MM+M15
!c vectors of MO alpha
      M18=MMd+M17
!c vectors of MO beta
!c
!c Density matrix  construction - For closed shell only <<<<=========
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (ecpmode) then !agregadas por Nick para lectura de ECP
           call lecturaECP()   !lee parametros
           CALL allocate_ECP() !allocatea la matriz de Fock de p-potenciales y el vector con los terminos de 1 electron sin corregir
           CALL ReasignZ() !reasigna las cargas de los nucleos removiendo la carga del core
        end if
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c DIMENSION TESTS -----------------------------------------------
!c
      Ndim=5*M*(M+1)/2+3*Md*(Md+1)/2+M+M*NCO!+M*Ngrid
      iprob=0
      if (Ndim.gt.ng2) then
        write(*,*) 'DIMENSION PROBLEMS WITH DYNAMICAL VECTOR NG2',Ndim,ng2
        iprob=1
      endif
!c
      if (iprob.eq.1) then
         write(*,*) 'PAUSE IS A DELETED FEATURE'
      endif

      ! Gets the number of occupied orbitals in a closed shell system (or
      ! Spin Up in an open shell system).
      call get_nco(Iz, natom, nco, NUNP, charge, OPEN, iostat)

      ! Allocates and initialises rhoalpha and rhobeta
      if(OPEN) then
        allocate(rhoalpha(M*(M+1)/2),rhobeta(M*(M+1)/2))
      else
        allocate(rhoalpha(1),rhobeta(1))
      endif
      rhoalpha(:) = 0.0d0
      rhobeta(:)  = 0.0d0


      ! Reads coefficient restart and builds density matrix. The MO
      ! coefficients are read in the same order as basis sets.
      ! Then vectors are put in dynamical allocation (to be used later)
      if (VCINP) then
         call g2g_timer_start('restart_read')

         open(unit=89, file=frestartin)
         if (rst_dens .gt. 0) then
            allocate(restart_dens(M, M))
            if (.not. OPEN) then
               call read_rho_restart(restart_dens, M, 89)
               call sprepack('L', M, RMM(1), restart_dens)
            else
               allocate(restart_adens(M,M), restart_bdens(M,M))
               call read_rho_restart(restart_adens, restart_bdens, M, 89)
               restart_dens = restart_adens + restart_bdens
               call sprepack('L', M, RMM(1)  , restart_dens)
               call sprepack('L', M, rhoalpha, restart_adens)
               call sprepack('L', M, rhobeta , restart_bdens)
               deallocate(restart_adens, restart_bdens)
            endif
            deallocate(restart_dens)
         else
            allocate(restart_dens(M, M), restart_coef(M, NCO))
            if (.not.OPEN) then
               call read_coef_restart(restart_coef, restart_dens, M, NCO, 89)

               kk = 0
               do k=1, NCO
               do i=1, M
                  kk = kk + 1
                  MO_coef_at(kk) = restart_coef(indexii(i), k)
               enddo
               enddo
            else
               NCOa = NCO
               NCOb = NCO + Nunp
               allocate(restart_coef_b(M, NCOb), restart_adens(M,M), &
                        restart_bdens(M,M))

               call read_coef_restart(restart_coef, restart_coef_b, &
                                      restart_dens, restart_adens,  &
                                      restart_bdens, M, NCOa, NCOb, 89)
               kk = 0
               do k=1, NCOa
               do i=1, M
                  kk = kk + 1
                  MO_coef_at(kk) = restart_coef(indexii(i), k)
               enddo
               enddo

               kk = 0
               do k=1, NCOb
               do i=1, M
                  kk = kk + 1
                  MO_coef_at_b(kk) = restart_coef_b(indexii(i), k)
               enddo
               enddo
               deallocate(restart_coef_b)
            endif

            ! Reorders by s, p, d.
            k = 0
            do j=1, M
            do i=j, M
               k = k + 1
               RMM(k)      = restart_dens(indexii(i), indexii(j))
               if (i.ne.j) then
                  RMM(k)      = RMM(k)*2.D0
               endif
            enddo
            enddo

            if (OPEN) then
               k = 0
               do j=1, M
               do i=j, M
                  k = k + 1
                  rhoalpha(k) = restart_adens(indexii(i), indexii(j))
                  rhobeta(k)  = restart_bdens(indexii(i), indexii(j))
                  if (i.ne.j) then
                     rhoalpha(k) = rhoalpha(k)*2.0D0
                     rhobeta(k)  = rhobeta(k)*2.0D0
                  endif
               enddo
               enddo
               deallocate(restart_adens, restart_bdens)
            endif
            deallocate(restart_dens, restart_coef)
         endif

         close(89)
         call g2g_timer_stop('restart_read')
      endif
      ! End of restart.

!c------- G2G Initialization ---------------------

      call g2g_parameter_init(NORM,natom,natom,ngDyn, &
                             rqm,Rm2,Iz,Nr,Nr2,Nuc, &
                             M,ncont,nshell,c,a, &
                             RMM,M5,M3,rhoalpha,rhobeta, &
                             NCO,OPEN,Nunp,nopt,Iexch, &
                             e_, e_2, e3, wang, wang2, wang3, &
			                    use_libxc, ex_functional_id, ec_functional_id)
      call summon_ghosts(Iz, natom, verbose)

      if (gpu_level .ne. 0) then
         call aint_parameter_init(Md, ncontd, nshelld, cd, ad, Nucd, af, RMM, &
                                  M9, M11, STR, FAC, rmax, Iz, gpu_level)
      endif
      allocate(X(M,4*M))

!--------------------------------------------------------------------------------------
      IF (number_restr.GT.0) THEN
! Distance Restrain parameters read
       ALLOCATE( restr_pairs(2,number_restr), restr_index(number_restr)&
       ,restr_k(number_restr), restr_w(number_restr),&
       restr_r0(number_restr))
        call read_restrain_params()
      END IF
!--------------------------------------------------------------------------------------

 100  format (A8)
 200  format ('  Basis set corresponding to Z = ', I3,' was not used.')
      return
end subroutine drive

subroutine get_nco(atom_Z, n_atoms, n_orbitals, n_unpaired, charge, open_shell,&
                   ext_status)
   use ghost_atoms_subs, only: adjust_ghost_charge
   implicit none
   integer, intent(in)  :: n_atoms, n_unpaired, charge, atom_Z(n_atoms)
   logical, intent(in)  :: open_shell
   integer, intent(out) :: n_orbitals, ext_status

   integer :: icount, nuc_charge, electrons

   nuc_charge = 0
   ext_status = 0
   do icount = 1, n_atoms
      nuc_charge = nuc_charge + atom_Z(icount)
   enddo
   call adjust_ghost_charge(atom_Z, n_atoms, nuc_charge)

   electrons = nuc_charge - charge
   if ((.not. open_shell) .and. (mod(electrons, 2) .ne. 0)) then
      write(*,'(A)') "  ERROR - DRIVE: Odd number of electrons in a &
                     &closed-shell calculation."
      write(*,'(A)') "  Please check system charge."
      ext_status = 1
      return
   endif

   if ((mod(electrons, 2) .eq. 0) .and. (mod(n_unpaired, 2) .ne. 0)) then
      write(*,'(A)') "  ERROR - DRIVE: Even number of electrons but odd &
                     &number of unpaired electrons."
      write(*,'(A)') "  Please check system charge or number of unpaired &
                     &electrons."
      ext_status = 2
      return
   endif

   if ((mod(electrons, 2) .ne. 0) .and. (mod(n_unpaired, 2) .eq. 0)) then
      write(*,'(A)') "  ERROR - DRIVE: Odd number of electrons but even &
                     &number of unpaired electrons."
      write(*,'(A)') "  Please check system charge or number of unpaired &
                     &electrons."
      ext_status = 2
      return
   endif

   n_orbitals = (electrons - n_unpaired) / 2
end subroutine get_nco
