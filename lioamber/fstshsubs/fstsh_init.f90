subroutine fstsh_init(  )
use fstsh_data  , only: tsh_file, all_states, type_coupling, Sovl_old, Sovl_now, &
                        a_old, c_old, r_old, C_scf_old, tsh_nucStep, tsh_Enstep, &
                        current_state, phases_old
use excited_data, only: root, nstates
use basis_data  , only: M,  max_c_per_atom
use garcha_mod  , only: natom, ntatom
   implicit none

   ! Open file and print header
   open(unit=tsh_file,file="tsh_lio.log")
   write(tsh_file,*) "=============================="
   write(tsh_file,*) " TRAJECTORY SURFFACE HOPPING"
   write(tsh_file,*) "=============================="

   ! Set and Allocated Variables
   all_states    = nstates + 1
   current_state = root + 1
   allocate(Sovl_old(M,M),Sovl_now(M,M))
   allocate(a_old(M,max_c_per_atom),c_old(M,max_c_per_atom))
   allocate(r_old(ntatom,3)); r_old = 0.0d0
   allocate(C_scf_old(M,M)); C_scf_old = 0.0d0
   allocate(phases_old(all_states)); phases_old = 1.0d0
   tsh_nucStep = 0
   Sovl_old = 0.0d0; Sovl_now = 0.0d0; a_old = 0.0d0; c_old = 0.0d0;

   ! Print Information
   write(tsh_file,*) " "
   write(tsh_file,*) "Information:"
   write(tsh_file,"(1X,A,I1,A)") "All states= ", all_states, " (1=GS)"
   write(tsh_file,"(1X,A,I1,A)") "Current state= ", current_state
   write(tsh_file,"(1X,A,I1,A)") "Type of Coupling= ", type_coupling
   write(tsh_file,"(1X,A,I2,A)") "Electronic Steps= ", tsh_Enstep

end subroutine fstsh_init
