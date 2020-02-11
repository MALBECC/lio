! Common    parameters for all of the quantum system.
! These are defaults to set some arrays.
!
! ntq: Number of atoms in part of the system treated quantum-mechanically.
! ntc: Number of 'solvent' classically treated molecules.
! ns : number of atoms in each solvent molecules.
! nt : Number total atoms.
! ng0: number of functions, ngd0 : number aux. funct.
! nl : number of primitives in a given contraction.

!   integer :: ntq, ntc, nss, nt, ng0, ng, nl
!   integer :: ngd0, ngd, ntqss, norbit, ngrid
!   parameter (ntq    = 200)
!   parameter (ntc    = 0)
!   parameter (nss    = 1)
!   parameter (nt     = ntq + nss * ntc)
!   parameter (ng0    = 100)
 !  parameter (ng     = ntq * ng0)
 !  parameter (nl     = 13)
!   parameter (ngd0   = 500)
!   parameter (ngd    = ntq * ngd0)
!   parameter (ntqss  = ntq + nss)
!   parameter (norbit = 800)
!   parameter (Ngrid  = 0)

! Things to update in gpu/cuda/excnum.cu in case this file is edited:
! FORTRAN_MAX_ATOMS = nt
! FORTRAN_NG = ng
! FORTRAN_NL = nl
