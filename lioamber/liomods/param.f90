! Common    parameters for all of the quantum system.
! These are defaults to set some arrays.
!
! ntq: Number of atoms in part of the system treated quantum-mechanically.
! ntc: Number of 'solvent' classically treated molecules.
! ns : number of atoms in each solvent molecules.
! nt : Number total atoms.
! ng0: number of functions, ngd0 : number aux. funct.
! nl : number of primitives in a given contraction.

   integer :: ntq, ntc, nss, nt, ng0, ng, nl
   integer :: ngd0, ngd, ntqss, norbit, ngrid
   parameter (ntq    = 200)
   parameter (ntc    = 0)
   parameter (nss    = 1)
   parameter (nt     = ntq + nss * ntc)
   parameter (ng0    = 100)
   parameter (ng     = ntq * ng0)
   parameter (nl     = 13)
   parameter (ngd0   = 500)
   parameter (ngd    = ntq * ngd0)
   parameter (ntqss  = ntq + nss)
   parameter (norbit = 800)
   parameter (Ngrid  = 0)

! Things to update in gpu/cuda/excnum.cu in case this file is edited:
! FORTRAN_MAX_ATOMS = nt
! FORTRAN_NG = ng
! FORTRAN_NL = nl


! Numerical    parameters:
   real(kind=8) :: PI, PI32, RPI, PI5, PI52
   real(kind=8) :: PISS, PIS32, RPIS, PIS5, PIS52
   parameter(PI32 = 5.56832799683170698D0, PI  = 3.14159265358979312D0, &
             RPI  = 1.77245385090551588D0, PI5 = 34.9868366552497108D0, &
             PI52 = 34.9868366552497108D0)
   parameter(PIS32 = 5.56832799683170698E0, PISS = 3.14159265358979312E0, &
             RPIS  = 1.77245385090551588E0, PIS5 = 34.9868366552497108E0, &
             PIS52 = 34.9868366552497108E0)

