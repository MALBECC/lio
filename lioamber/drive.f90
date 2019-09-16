!------------------------------------------------------------------------------!
! READS AND DRIVES SUBROUTINE                                                  !
! reads and prepare basis sets and geometry to perform electronic structure
! calculations using density functional theory and gaussian basis sets for
! the expansion of MO, charge density and exchange correlation potential
!
! Started 13 March 1992, Dario Estrin
!
!------------------------------------------------------------------------------!

subroutine drive(iostat)
   use garcha_mod, only: X, rhoalpha, rhobeta,  charge, e_, e_2, e3, Pmat_vec, &
                         fcoord, fmulliken, natom, frestart, Iexch, frestartin,&
                         NCO, npas, Nr, Nr2, wang, wang2, wang3, VCINP, OPEN,  &
                         Iz, Rm2, rqm, Nunp, restart_freq, writexyz, gpu_level,&
                         number_restr, restr_pairs, restr_index, restr_k,      &
                         restr_w, restr_r0, mulliken, MO_coef_at, MO_coef_at_b,&
                         use_libxc, ex_functional_id, ec_functional_id,        &
                         pi, Fmat_vec, Fmat_vec2, Ginv_vec, Hmat_vec, becke
   use basis_data, only: nshell, nshelld, ncont, ncontd, indexii, a, c, ad, cd,&
                         af, M, Md, rmax, norm, nuc, nucd
   use ECP_mod     , only: ecpmode
   use fileio      , only: read_coef_restart, read_rho_restart
   use td_data     , only: td_do_pop
   use fileio_data , only: verbose, rst_dens
   use math_data   , only: FAC, STR
   use liosubs_math, only: init_math
   use ghost_atoms_subs, only: summon_ghosts
   use tbdft_data  , only: MTB, tbdft_calc


   implicit none
   integer, intent(inout) :: iostat

   double precision, allocatable :: restart_coef(:,:), restart_coef_b(:,:), &
                                    restart_dens(:,:), restart_adens(:,:),  &
                                    restart_bdens(:,:)
   integer :: NCOa, NCOb
   integer :: M_f, NCO_f,i0
   integer :: icount, kcount, jcount

   ! Calls generator of table for incomplete gamma functions
   call init_math()
   call GRIDLIO()

   ! Opens files for IO
   if (writexyz) open(unit = 18,file = fcoord)
   if (restart_freq .gt. 0) open(unit = 88, file = frestart)

   if (ecpmode) then !agregadas por Nick para lectura de ECP
      call lecturaECP()   !lee parametros
      call allocate_ECP() !allocatea la matriz de Fock de p-potenciales y el vector con los terminos de 1 electron sin corregir
      call ReasignZ() !reasigna las cargas de los nucleos removiendo la carga del core
   end if

   ! Gets the number of occupied orbitals in a closed shell system (or
   ! Spin Up in an open shell system).
   call get_nco(Iz, natom, nco, NUNP, charge, OPEN, iostat)

!TBDFT: Updating M and NCO for TBDFT calculations
    if (tbdft_calc) then
       M_f = M+2*MTB
       NCO_f=NCO+MTB
       i0 = MTB
    else
       M_f = M
       NCO_f=NCO
       i0=0
    end if

   ! Allocates and initialises rhoalpha and rhobeta
   if (OPEN) then
      allocate(rhoalpha(M*(M+1)/2),rhobeta(M*(M+1)/2))
   else
      allocate(rhoalpha(1),rhobeta(1))
   endif
   rhoalpha(:) = 0.0D0
   rhobeta(:)  = 0.0D0

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
            call sprepack('L', M, Pmat_vec, restart_dens)
         else
            allocate(restart_adens(M,M), restart_bdens(M,M))
            call read_rho_restart(restart_adens, restart_bdens, M, 89)
            restart_dens = restart_adens + restart_bdens
            call sprepack('L', M, Pmat_vec, restart_dens)
            call sprepack('L', M, rhoalpha, restart_adens)
            call sprepack('L', M, rhobeta , restart_bdens)
            deallocate(restart_adens, restart_bdens)
         endif
         deallocate(restart_dens)
      else
         allocate(restart_dens(M, M), restart_coef(M, NCO_f))
         MO_coef_at = 0.0D0
         if (.not. OPEN) then
            call read_coef_restart(restart_coef, restart_dens, M, NCO_f, 89)
            do icount = 1, M
            do jcount = 1, NCO_f
               MO_coef_at(i0+icount, jcount) = &
                                restart_coef(indexii(icount), jcount)
            enddo
            enddo
         else
            MO_coef_at_b = 0.0D0
            NCOa = NCO_f
            NCOb = NCO_f + Nunp
            allocate(restart_coef_b(M, NCOb), restart_adens(M,M), &
                     restart_bdens(M,M))

            call read_coef_restart(restart_coef, restart_coef_b, &
                                   restart_dens, restart_adens,  &
                                   restart_bdens, M, NCOa, NCOb, 89)
            do icount = 1, M
            do jcount = 1, NCOa
               MO_coef_at(i0+icount, jcount) = &
                                restart_coef(indexii(icount), jcount)
            enddo
            enddo

            do icount = 1, M
            do jcount = 1, NCOb
               MO_coef_at_b(i0+icount, jcount) = &
                                restart_coef_b(indexii(icount), jcount)
            enddo
            enddo
            deallocate(restart_coef_b)
         endif

         ! Reorders by s, p, d.
         kcount = 0
         do jcount = 1     , M
         do icount = jcount, M
            kcount = kcount + 1
            Pmat_vec(kcount) = restart_dens(indexii(icount),indexii(jcount))
            if (icount .ne. jcount) then
               Pmat_vec(kcount) = Pmat_vec(kcount) * 2.0D0
            endif
         enddo
         enddo

         if (OPEN) then
            kcount = 0
            do jcount = 1     , M
            do icount = jcount, M
               kcount = kcount + 1
               rhoalpha(kcount) = restart_adens(indexii(icount),indexii(jcount))
               rhobeta(kcount)  = restart_bdens(indexii(icount),indexii(jcount))
               if (icount .ne. jcount) then
                  rhoalpha(kcount) = rhoalpha(kcount) * 2.0D0
                  rhobeta(kcount)  = rhobeta(kcount)  * 2.0D0
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

   ! G2G and AINT(GPU) Initializations
   call g2g_parameter_init(NORM, natom, natom, M, rqm, Rm2, Iz, Nr, Nr2,  &
                           Nuc, M, ncont, nshell, c, a, Pmat_vec, Fmat_vec,   &
                           Fmat_vec2, rhoalpha, rhobeta, NCO, OPEN, Nunp, 0,  &
                           Iexch, e_, e_2, e3, wang, wang2, wang3, use_libxc, &
                           ex_functional_id, ec_functional_id, becke)
   call summon_ghosts(Iz, natom, verbose)

   if (gpu_level .ne. 0) call aint_parameter_init(Md, ncontd, nshelld, cd, ad, &
                                                  Nucd, af, Ginv_vec, Hmat_vec,&
                                                  STR, FAC, rmax, Iz, gpu_level)
  ! TO-DO: Relocate this.
  allocate(X(M , 4*M))
  npas = 0

  ! Restraints initialization
  if (number_restr .gt. 0) then
     allocate(restr_pairs(2,number_restr), restr_index(number_restr), &
              restr_k(number_restr), restr_w(number_restr), restr_r0(number_restr))
     call read_restrain_params()
  endif
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
