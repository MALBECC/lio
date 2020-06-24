!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%% LJ SWITCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This module contains the subroutines needed for a variable LJ parameter model!
! for QM atoms in a QM/MM scheme. See lj_switch_data.f90 for data definitions. !
!                                                                              !
! This implementation is self-consistent, meaning that the model used also     !
! contributes with Fock matrix elements; this is done in order to avoid        !
! Mulliken charge derivatives in gradient calculation, which can be a mess.    !
! See Kuechler's work (http://dx.doi.org/10.1063/1.4937166) for a similar idea.!
!                                                                              !
! In essence, LJ parameters for a given QM atom are interpolated between two   !
! states with a given reference Mulliken atomic charge. Special care must be   !
! taken when defining those states. The interpolation is performed with a Fermi!
! distribution function, as shown below:                                       !     
!                   sigma = sigma_1 + (sigma_2 - sigma_1) /                    !  
!                                     (1 + exp(-k * ( q - (q2 + q1)/2 ))       !
!                                                                              !
! The same calculation scheme applies to epsilon.                              !
! When using this module, the following input structure should be added to the !
! LIO input file:                                                              !
!                                                                              !
! {LJSWITCH}                                                                   !
!   index charge_1 charge_2 sigma_1 sigma_2 epsilon_1 epsilon_2                !
!    ...                                                                       !
!   index charge_1 charge_2 sigma_1 sigma_2 epsilon_1 epsilon_2                !
! {END}                                                                        !
!                                                                              !    
! Where index is the index of the QM atom, and charge, sigma and epsilon are   !
! those correspoding to references states 1 and 2. Each line belongs to a      !
! different atom, and no blank lines should be found in that block.            !
! Sigma is to be given in angstrom and epsilon in kJ/mol.                      !
!                                                                              !
! Through out the module, values are stored in atomic units. In particular,    !
! epsilon and its derivatives are always stored as *4, in order to ease LJ     !
! calculations. Subroutines found in ljs_mm_interface are actually OUTSIDE     !
! this module, so that they can be accessed from other programs.               !
!                                                                              !
! Within this module, the following subroutines can be found:                  !
!  * ljs_initialise:     Initialises QM-only data for this module, receiving   !
!                        the atomic numbers (iz) and atom-of-function (nuc)    !
!                        arrays as input and returning a boolean as output if  !
!                        this module is indeed used.                           !
!  * ljs_finalise:       Deallocates this modules data structures.             !
!  * ljs_input_read:     Reads the input file LJSWITCH segment, as specified   !
!                        above. It receives the input file UID as input.       !
!  * ljs_add_fock_terms: Adds the LJS terms to Fock matrix elements. Receives  !
!                        both density and overlap matrices as input, and       !
!                        outputs Fock matrix terms and energy. The open shell  !
!                        variant receives density matrices alpha and beta, and !
!                        outputs Fock alpha and beta.                          !
!  * ljs_gradients_qmmm: Calculates the cartesian gradient terms for LJS       !
!                        contributions. It receives the global (QM+MM)         !
!                        positions array (r) and the number of QM atoms (natom)!
!                        as inputs, returning the QM-region gradients and      !
!                        MM-region gradients as separate outputs.              !
!  * ljs_get_energy:     Calculates the LJS energy terms. This subroutine is   !
!                        not actually used anywhere, but given as an extra tool!
!                        if needed.                                            !
!  * ljs_get_dEdQ:       This subroutine is PRIVATE within the module, accessed!
!                        only by_ljs_add_fock_terms. It calculates the term    !
!                        dEdQ for a given QM atom (received as input), and     !
!                        returns the energy contribution and dEdQ term for said!
!                        atom.                                                 !
!                                                                              !
! In addition, the following subroutines are included in ljs_mm_interface.f90: !
!  * ljs_set_params:   This subroutine receives the number of different MM atom!
!                      types, and arrays containing the epsilon and sigma      !
!                      values for those types (in a.u.). Then copies these data!
!                      to the internal structures.                             !
!  * ljs_settle_mm:    This subroutine receives the number of QM atoms, the    !
!                      number of MM atoms, the position of QM atoms, the       !
!                      the positions of MM atoms, an array containing the MM   !
!                      atom types for the QM atoms, and an array containing the!
!                      MM atom types for the MM atoms. Then, it reorganises    !
!                      this information within the internal data structure.    !
!                      Atoms positions are expected to be in a.u.              !
!  * ljs_substract_mm: This calculates the classical energy correction for the !
!                      MM software, in case said software calculates the LJ    !
!                      contributions before entering QM/MM. As input, it       !
!                      receives the number of QM atoms, the number of MM atoms !
!                      the position of QM atoms and the position of MM atoms   !
!                      (both in a.u.). It returns the corrections to energy,   !
!                      MM-region gradients and QM-region gradients.            !
!                                                                              !
! For further sample usages, check the lio/tests/AMBER and /lioamber/utests.   !
! Also, take a look the the LIO/AMBER interface in lio/dat/ambermod.           !
!                                                                              !
! First written by: Federico Pedron, Jun/2020                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../datatypes/datatypes.fh"
module LJ_switch
   implicit none
   private
   public :: ljs_initialise
   public :: ljs_input_read
   public :: ljs_finalise

   public :: ljs_gradients_qmmm
   public :: ljs_get_energy

   public :: ljs_add_fock_terms
   public :: ljs_add_fock_terms_op
   
   
contains

#include "ljs_dEdQ.f90"
#include "ljs_fock_terms.f90"
#include "ljs_init_end.f90"
#include "ljs_input.f90"
#include "ljs_gradients_qmmm.f90"

end module LJ_switch
