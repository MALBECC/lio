!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% PYTHON_INTERFACE.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Funciones de utilidad para facilitar el acceso desde Python                  !
! Proporciona funciones auxiliares simples para interfaz Python               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

module python_interface
   use iso_c_binding
   implicit none
   
   public :: python_hello, python_get_version, python_lio_init, python_lio_calculate, &
            python_lio_scf, python_lio_gradients

contains

!==============================================================================!
! PYTHON_HELLO: Función de prueba simple
!==============================================================================!
subroutine python_hello() bind(C, name="python_hello")
   implicit none
   
   write(*,*) "👋 Hola desde el módulo Fortran de LIO!"
   write(*,*) "✅ La interfaz Python-Fortran está funcionando correctamente"

end subroutine python_hello

!==============================================================================!
! PYTHON_GET_VERSION: Obtener información de versión
!==============================================================================!
subroutine python_get_version(major, minor, patch) bind(C, name="python_get_version")
   implicit none
   
   ! Parámetros
   integer(c_int), intent(out) :: major, minor, patch
   
   major = 2020
   minor = 1
   patch = 0
   
   write(*,*) "📋 Versión LIO: ", major, ".", minor, ".", patch

end subroutine python_get_version

!==============================================================================!
! PYTHON_LIO_INIT: Inicialización completa del sistema LIO
!==============================================================================!
subroutine python_lio_init(natom, atomic_numbers, coordinates, charge, &
                           basis_name, iostat) bind(C, name="python_lio_init")
   use garcha_mod, only: natom_global => natom, ntatom, nsol, Iz, r, rqm, pc, charge_global => charge
   use basis_data, only: basis_set, fitting_set, int_basis
   implicit none
   
   ! Declaraciones externas
   external :: lio_defaults, init_lio_common
   
   ! Parámetros de entrada
   integer(c_int), intent(in), value :: natom
   integer(c_int), intent(in) :: atomic_numbers(natom)
   real(c_double), intent(in) :: coordinates(3*natom)  ! x1,y1,z1,x2,y2,z2,...
   integer(c_int), intent(in), value :: charge
   character(c_char), intent(in) :: basis_name(*)
   integer(c_int), intent(out) :: iostat
   
   ! Variables locales
   character(len=20) :: basis_str
   integer :: i, null_pos
   
   iostat = 0
   
   write(*,*) "🔧 python_lio_init: Inicializando sistema LIO completo"
   write(*,*) "   Átomos:", natom
   write(*,*) "   Carga:", charge
   
   ! Convertir string de basis desde C
   null_pos = 1
   do while (basis_name(null_pos) /= c_null_char .and. null_pos <= 20)
      basis_str(null_pos:null_pos) = basis_name(null_pos)
      null_pos = null_pos + 1
   end do
   basis_str(null_pos:) = ' '
   
   write(*,*) "   Basis set:", trim(basis_str)
   
   ! Configurar defaults de LIO
   call lio_defaults()
   
   ! Configurar variables globales
   natom_global = natom
   ntatom = natom  ! Sin solvente
   nsol = 0
   charge_global = charge
   
   ! Configurar basis sets usando el parámetro de entrada
   ! Si basis_str está vacío, usar valor por defecto
   if (len_trim(basis_str) > 0) then
      basis_set = trim(basis_str)
      write(*,*) "   Usando basis set especificado:", trim(basis_set)
   else
      ! Mantener valor por defecto de basis_data (DZVP)
      write(*,*) "   Usando basis set por defecto:", trim(basis_set)
   endif
   
   ! fitting_set mantiene su valor por defecto (DZVP Coulomb Fitting)
   write(*,*) "   Fitting set:", trim(fitting_set)
   int_basis = .true.  ! Usar basis set interno
   
   ! Alocar y configurar arrays de átomos
   if (allocated(Iz)) deallocate(Iz)
   if (allocated(r)) deallocate(r)
   if (allocated(rqm)) deallocate(rqm)
   if (allocated(pc)) deallocate(pc)
   
   allocate(Iz(natom))
   allocate(r(natom, 3))
   allocate(rqm(natom, 3))
   allocate(pc(natom))
   
   ! Asignar números atómicos
   do i = 1, natom
      Iz(i) = atomic_numbers(i)
   end do
   
   ! Asignar coordenadas (convertir de array 1D a 2D)
   do i = 1, natom
      r(i, 1) = coordinates(3*(i-1) + 1)
      r(i, 2) = coordinates(3*(i-1) + 2)
      r(i, 3) = coordinates(3*(i-1) + 3)
      
      ! Copiar coordenadas QM (mismo que r para cálculo completo)
      rqm(i, 1) = r(i, 1)
      rqm(i, 2) = r(i, 2)
      rqm(i, 3) = r(i, 3)
   end do
   
   ! Inicializar cargas puntuales (cero para átomos QM)
   pc = 0.0d0
   
   write(*,*) "   Coordenadas configuradas:"
   do i = 1, natom
      write(*,'(A,I3,A,I3,A,3F10.6)') "     Átomo", i, " (Z=", Iz(i), "):", r(i,:)
   end do
   
   ! Llamar inicialización común de LIO
   call init_lio_common(natom, Iz, nsol, 0)
   
   write(*,*) "✅ python_lio_init: Sistema LIO inicializado completamente"

end subroutine python_lio_init

!==============================================================================!
! PYTHON_LIO_CALCULATE: Calcular integrales con sistema inicializado
!==============================================================================!
subroutine python_lio_calculate(energy, iostat) bind(C, name="python_lio_calculate")
   use subm_int1, only: int1
   use garcha_mod, only: natom, Iz, r
   use basis_data, only: M  ! Usar el número real de funciones de base
   implicit none
   
   ! Parámetros
   real(c_double), intent(out) :: energy
   integer(c_int), intent(out) :: iostat
   
   ! Variables locales para int1
   real*8 :: En
   real*8, allocatable :: Fmat(:), Hmat(:), Smat(:,:)
   real*8, allocatable :: d(:,:), r_local(:,:)
   integer :: MM, i, j, MM2
   
   iostat = 0
   energy = 0.0d0
   
   write(*,*) "🔧 python_lio_calculate: Calculando integrales"
   
   ! Usar el número real de funciones de base desde basis_data
   MM = M  ! M se configura durante init_lio_common
   MM2 = MM * MM
   
   write(*,*) "   Número real de funciones de base (M):", MM
   
   ! Alocar matrices temporales
   allocate(Fmat(MM2))      ! Array 1D
   allocate(Hmat(MM2))      ! Array 1D
   allocate(Smat(MM, MM))   ! Array 2D
   allocate(d(natom, natom))
   allocate(r_local(natom, 3))
   
   ! Copiar coordenadas
   do i = 1, natom
      r_local(i, :) = r(i, :)
   end do
   
   ! Inicializar matrices
   Fmat = 0.0d0
   Hmat = 0.0d0
   Smat = 0.0d0
   
   ! Calcular matriz de distancias al cuadrado entre átomos
   do i = 1, natom
      do j = 1, natom
         if (i == j) then
            d(i, j) = 0.0d0
         else
            d(i, j) = (r_local(i,1) - r_local(j,1))**2 + &
                      (r_local(i,2) - r_local(j,2))**2 + &
                      (r_local(i,3) - r_local(j,3))**2
         endif
      end do
   end do
   
   write(*,*) "   Llamando int1 con:", natom, "átomos y", MM, "funciones base"
   
   ! Llamar int1 con todos los parámetros requeridos
   call int1(En, Fmat, Hmat, Smat, d, r_local, Iz, natom, natom)
   
   energy = En
   
   write(*,*) "✅ Energía nuclear calculada:", energy, "Hartree"
   
   ! Limpiar memoria
   deallocate(Fmat, Hmat, Smat, d, r_local)

end subroutine python_lio_calculate

!==============================================================================!
! PYTHON_LIO_SCF: Realizar cálculo SCF completo
!==============================================================================!
subroutine python_lio_scf(total_energy, iostat) bind(C, name="python_lio_scf")
   use garcha_mod, only: natom, Iz, r, NCO, OPEN, npas, converge, noconverge, &
                         Smat, RealRho, sqsm, Eorbs, Eorbs_b
   use basis_data, only: MM, M, nshell, Nuc
   use tbdft_data, only: MTB, tbdft_calc
   use typedef_operator, only: operator
   use initial_guess_subs, only: get_initial_guess
   use SCF_aux, only: standard_coefs
   use fileio_data, only: verbose
   implicit none
   
   ! Declaraciones externas
   external :: SCF
   
   ! Parámetros
   real(c_double), intent(out) :: total_energy
   integer(c_int), intent(out) :: iostat
   
   ! Variables locales
   real*8 :: E
   type(operator) :: rho_aop, fock_aop, rho_bop, fock_bop
   real*8, allocatable :: Xmat(:,:), Hvec(:), Rhovec(:)
   real*8, allocatable :: Rhoalpha(:), Rhobeta(:)
   integer :: i, j, k, M_f
   
   iostat = 0
   total_energy = 0.0d0
   E = 0.0d0
   
   ! Configurar verbose para obtener más información de debug
   verbose = 4
   
   ! Configurar variables del SCF que pueden estar sin inicializar
   npas = 1
   converge = 0
   noconverge = 0
   
   ! Calcular M_f para TBDFT compatibility
   M_f = M
   if (tbdft_calc /= 0) M_f = M + MTB
   
   write(*,*) "🔧 python_lio_scf: Iniciando cálculo SCF"
   write(*,*) "   Funciones de base (M):", M
   write(*,*) "   M_f (M + MTB):", M_f
   write(*,*) "   Matriz size (MM):", MM
   write(*,*) "   Orbitales ocupados (NCO):", NCO
   write(*,*) "   Open shell:", OPEN
   write(*,*) "   Verbose level:", verbose
   
   ! CRÍTICO: Alocar matrices globales que SCF necesita (como en liomain.f90)
   if (.not.allocated(Smat))      allocate(Smat(M,M))
   if (.not.allocated(RealRho))   allocate(RealRho(M,M))
   if (.not.allocated(sqsm))      allocate(sqsm(M,M))
   if (.not.allocated(Eorbs))     allocate(Eorbs(M_f))
   if (.not.allocated(Eorbs_b))   allocate(Eorbs_b(M_f))
   
   ! Alocar matrices
   allocate(Xmat(M, M))
   allocate(Hvec(MM))
   allocate(Rhovec(MM))
   allocate(Rhoalpha(MM))
   allocate(Rhobeta(MM))
   
   ! Inicializar operadores
   allocate(rho_aop%data_AO(M, M))
   allocate(fock_aop%data_AO(M, M))
   if (OPEN) then
      allocate(rho_bop%data_AO(M, M))
      allocate(fock_bop%data_AO(M, M))
   endif
   
   ! Inicializar matrices en cero
   Xmat = 0.0d0
   Hvec = 0.0d0
   Rhovec = 0.0d0
   rho_aop%data_AO = 0.0d0
   fock_aop%data_AO = 0.0d0
   
   ! Para shell cerrado, usar matriz identidad como guess inicial de transformación
   do i = 1, M
      Xmat(i, i) = 1.0d0
   end do
   
   write(*,*) "   Obteniendo guess inicial..."
   
   ! Obtener guess inicial para la matriz de densidad
   if (OPEN) then
      call get_initial_guess(M, MM, NCO, NCO, Xmat, Hvec, Rhovec, Rhoalpha, &
                            Rhobeta, .true., natom, Iz, nshell, Nuc)
      
      ! Convertir vectores a matrices de operadores
      k = 0
      do i = 1, M
         do j = 1, i
            k = k + 1
            rho_aop%data_AO(i, j) = Rhoalpha(k)
            rho_aop%data_AO(j, i) = Rhoalpha(k)
            rho_bop%data_AO(i, j) = Rhobeta(k)
            rho_bop%data_AO(j, i) = Rhobeta(k)
         end do
      end do
   else
      call get_initial_guess(M, MM, NCO, NCO, Xmat, Hvec, Rhovec, Rhoalpha, &
                            Rhobeta, .false., natom, Iz, nshell, Nuc)
      
      ! Convertir vector a matriz de operador
      k = 0
      do i = 1, M
         do j = 1, i
            k = k + 1
            rho_aop%data_AO(i, j) = Rhovec(k)
            rho_aop%data_AO(j, i) = Rhovec(k)
         end do
      end do
   endif
   
   write(*,*) "   Ejecutando ciclo SCF..."
   
   ! Llamar función SCF principal
   call SCF(E, fock_aop, rho_aop, fock_bop, rho_bop)
   
   total_energy = E
   
   write(*,*) "✅ Cálculo SCF completado"
   write(*,*) "   Energía total:", total_energy, "Hartree"
   
   ! Limpiar solo memoria que nosotros alocamos
   deallocate(Xmat, Hvec, Rhovec)
   deallocate(Rhoalpha, Rhobeta)
   
   ! NO limpiar los objetos operator aquí, pueden estar siendo manejados por LIO
   ! if (allocated(rho_aop%data_AO)) deallocate(rho_aop%data_AO)
   ! if (allocated(fock_aop%data_AO)) deallocate(fock_aop%data_AO)
   ! if (allocated(rho_bop%data_AO)) deallocate(rho_bop%data_AO)
   ! if (allocated(fock_bop%data_AO)) deallocate(fock_bop%data_AO)

end subroutine python_lio_scf

!==============================================================================!
! PYTHON_LIO_GRADIENTS: Calcular gradientes (fuerzas negativas)
!==============================================================================!
subroutine python_lio_gradients(gradients, iostat) bind(C, name="python_lio_gradients")
   use garcha_mod, only: natom
   implicit none
   
   ! Declaraciones externas
   external :: dft_get_qm_forces
   
   ! Parámetros
   real(c_double), intent(out) :: gradients(3*natom)  ! Array 1D: dx1,dy1,dz1,dx2,dy2,dz2,...
   integer(c_int), intent(out) :: iostat
   
   ! Variables locales
   real*8, allocatable :: forces(:,:)  ! Array 2D para dft_get_qm_forces
   integer :: i, idx, alloc_stat
   
   iostat = 0
   
   write(*,*) "🔧 python_lio_gradients: Calculando gradientes"
   write(*,*) "   Átomos:", natom
   
   ! Alocar array para fuerzas (3, natom) con verificación de errores
   allocate(forces(3, natom), stat=alloc_stat)
   if (alloc_stat /= 0) then
      write(*,*) "❌ Error alocando memoria para fuerzas"
      iostat = -1
      return
   end if
   
   ! Inicializar array de fuerzas
   forces = 0.0d0
   
   write(*,*) "   Llamando dft_get_qm_forces..."
   
   ! Calcular fuerzas usando la función de LIO
   call dft_get_qm_forces(forces)
   
   write(*,*) "   Fuerzas calculadas exitosamente"
   
   ! Convertir de matriz 2D a array 1D para Python
   ! Los gradientes son las fuerzas negativas: grad = -force
   idx = 0
   do i = 1, natom
      idx = idx + 1
      gradients(idx) = -forces(1, i)  ! -fx
      idx = idx + 1
      gradients(idx) = -forces(2, i)  ! -fy
      idx = idx + 1
      gradients(idx) = -forces(3, i)  ! -fz
   end do
   
   write(*,*) "✅ Gradientes calculados para", natom, "átomos"
   write(*,*) "   Norma del gradiente:", sqrt(sum(gradients**2))
   
   ! Limpiar memoria
   deallocate(forces)

end subroutine python_lio_gradients

!==============================================================================!
! PYTHON_LIO_FINALIZE: Limpieza completa del sistema LIO (versión segura)
!==============================================================================!
subroutine python_lio_finalize() bind(C, name="python_lio_finalize")
   use garcha_mod, only: Smat, RealRho, sqsm, Eorbs, Eorbs_b, &
                         MO_coef_at, MO_coef_at_b, r, v, rqm, Em, Rm, &
                         pc, Iz, d, Fmat_vec, Fmat_vec2, Pmat_vec,    &
                         Hmat_vec, Ginv_vec, Gmat_vec, Pmat_en_wgt
   use properties    , only: properties_finalise
   use ECP_mod       , only: ecpmode
   use basis_subs    , only: basis_deinit
   use converger_subs, only: converger_finalise
   use dftd3         , only: dftd3_finalise
   use lj_switch     , only: ljs_finalise
   use ceed_subs     , only: ceed_finalize
   implicit none
   
   ! Declaración externa
   external :: generalECP
   
   write(*,*) "🧹 python_lio_finalize: Limpieza segura del sistema LIO"
   
   ! Deallocar matrices opcionales (con verificación)
   if (allocated(Smat))         deallocate(Smat)
   if (allocated(RealRho))      deallocate(RealRho)
   if (allocated(sqsm))         deallocate(sqsm)
   if (allocated(Eorbs))        deallocate(Eorbs)
   if (allocated(Eorbs_b))      deallocate(Eorbs_b)
   if (allocated(MO_coef_at))   deallocate(MO_coef_at)
   if (allocated(MO_coef_at_b)) deallocate(MO_coef_at_b)
   
   ! Deallocar vectores con verificación segura
   if (allocated(Fmat_vec))    deallocate(Fmat_vec)
   if (allocated(Fmat_vec2))   deallocate(Fmat_vec2)
   if (allocated(Pmat_vec))    deallocate(Pmat_vec)
   if (allocated(Hmat_vec))    deallocate(Hmat_vec)
   if (allocated(Ginv_vec))    deallocate(Ginv_vec)
   if (allocated(Gmat_vec))    deallocate(Gmat_vec)
   if (allocated(Pmat_en_wgt)) deallocate(Pmat_en_wgt)
   
   ! Deallocar arrays básicos con verificación
   if (allocated(r))   deallocate(r)
   if (allocated(v))   deallocate(v)
   if (allocated(rqm)) deallocate(rqm)
   if (allocated(Em))  deallocate(Em)
   if (allocated(Rm))  deallocate(Rm)
   if (allocated(pc))  deallocate(pc)
   if (allocated(Iz))  deallocate(Iz)
   if (allocated(d))   deallocate(d)
   
   ! Finalizar submódulos de forma segura
   call basis_deinit()
   
   ! ECP solo si está habilitado
   if (ecpmode) call generalECP(4)
   
   ! Finalizar GPU/G2G de forma segura
   call aint_deinit()
   call g2g_timer_summary()
   call g2g_deinit()
   
   ! Otros módulos
   call converger_finalise()
   call dftd3_finalise()
   call ljs_finalise()
   call ceed_finalize()
   call properties_finalise()
   
   write(*,*) "✅ python_lio_finalize: Sistema LIO limpiado correctamente"

end subroutine python_lio_finalize

end module python_interface
