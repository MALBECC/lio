#include "../datatypes/datatypes.fh"
module LJ_switch
   implicit none
   private
   public :: doing_ljs
   public :: ljs_initialise
   public :: ljs_input_read
   public :: ljs_finalise

   public :: ljs_settle_mm
   public :: ljs_substract_mm
   public :: ljs_gradients_qmmm
   public :: ljs_get_energy

   public :: ljs_add_fock_terms
   public :: ljs_add_fock_terms_op
   
   
contains

function doing_ljs() result(is_doing)
   use LJ_switch_data, only: n_lj_atoms
   
   implicit none
   logical :: is_doing

   is_doing = .false.
   if (n_lj_atoms > 0) is_doing = .true.
   return
end function doing_ljs

#include "ljs_dEdQ.f90"
#include "ljs_fock_terms.f90"
#include "ljs_init_end.f90"
#include "ljs_input.f90"
#include "ljs_mm_interface.f90"

end module LJ_switch
