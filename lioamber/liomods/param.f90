C
C     DATOS COMUNES A TODAS LAS RUTINAS DEL
C     SISTEMA QUANTICO (SOLUTO)
C

c PARAMETER FILE - FOR BASIS FUNCTIONS, ETC
c ntq : Number of atoms in part of the system treated
c quantum-mechanically
c
c ntc : Number of 'solvent' classically treated molecules
c ns : number of atoms in each solvent molecule
c
c nt : Number total atoms,
c ng0 : number of functions, ngd0 : number aux. funct.
c nl : number of primitives in a given contraction
c


      integer ntq, ntc, nss, nt, ng0, ng, nl
      integer ngd0, ngd, ntqss, norbit, ngrid
      parameter (ntq=200,ntc=0,nss=1)
      parameter (nt=ntq+nss*ntc)
      parameter (ng0=100,ng=ntq*ng0,nl=13)
cng0 era 100 y ngd0 500, Nick
      parameter (ngd0=500,ngd=ntq*ngd0)
      parameter (ntqss=ntq+nss)
      parameter (norbit=800,Ngrid=0)

c !!!!!!! Actualizar en gpu/cuda/excnum.cu !!!!!
c FORTRAN_MAX_ATOMS = nt
c FORTRAN_NG = ng
c FORTRAN_NL = nl
c !!!!!!!!
