c This subroutine calculates the solvent-solvent 
c INDUCED DIPOLES, evaluated using a modification of Sprik's method
c of 4 point charges to represent the induced dipole
c
c DOES NOT EVALUATE INTERACTION ENERGIES
c the output of this subroutine is used in order to evaluate
c energies in mmsol.f
c also the same induced point charges are used in intsol.f
c
c In order to evaluate the induced charges, it is necessary to
c calculate the electrical field due to the other solvent
c molecules (done here since it is trivial) and the electrical
c field fue to the solvent (done with a call to the
c subroutine efield.f)
c
c To be used when the option sol$ (in namelist &int0) is selected.
c
c
c Dario Estrin, Buenos Aires February 1994.
c----------------------------------------------------------------------------
c
      subroutine mmpol(NORM,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,
     >   natom,Nsol,natsol,Iz,alpha,pc,pci,r)
c
c Nsol : # solvent molecules
c natsol : # atoms in each solvent molecule
c pc : partial charges
c rs : nuclear coordinates 
c Es : solvent contributions to the energy
c
      INCLUDE 'param'
      parameter (pi=3.14159265358979312D0)
c     parameter (a2=0.2672774D0)
      parameter (a2=0.10D0)
      implicit real*8 (a-h,o-z)
      dimension pc(natsol),r(nt,3)
      dimension Iz(natom)
c
      dimension pci(nt,4),alpha(nss)
      dimension px(nt),py(nt),pz(nt),xi(3),Ef(3)
      dimension rt(4,3)
c
      common /tetra/ rt
c-----------------------------------------------------------------
c UNITS/CONSTANTS
c
c constant to go from charges (electron) and distances (angstroms)
c to kcal/mol in electrostatic terms
c solvent internal coordinates force constants
c in kcal/mol, angstroms and degrees
c water data
c in kcal/mol, angstroms and degrees
c     data RF,Ra,d0,theta0 /450.0D0,55.0D0,0.9572D0,104.52D0/
c in hartress, bohr and degrees
      data RF,Ra,d0,theta0 /0.71701721D0,0.087635437D0,1.8088466D0,
     >                       104.52D0/
c kcal/mol and angstroms
c     C1=336.26649
c hartrees and bohr
      C1 = 1.0D0
c The Lennard-Jones parameters should be given in consistent units
c (kcal/mol angstroms if using the first set, or hartrees, bohr
c which is recommended).
c---------------------------------------------------------------
c
c rt tetraedral coordinates of the 4 point charges
c displacements from center
       rt(1,1)=-a2
       rt(1,2)=-a2
       rt(1,3)=a2
       rt(2,1)=a2
       rt(2,2)=a2
       rt(2,3)=a2
       rt(3,1)=a2
       rt(3,2)=-a2
       rt(3,3)=-a2
       rt(4,1)=-a2
       rt(4,2)=a2
       rt(4,3)=-a2
c-----------------
c      rt(1,1)=0.0
c      rt(1,2)=0.0
c      rt(1,3)=0.0
c      rt(2,1)=a2
c      rt(2,2)=0.0
c      rt(2,3)=0.0
c      rt(3,1)=0.0
c      rt(3,2)=a2
c      rt(3,3)=0.0
c      rt(4,1)=0.0
c      rt(4,2)=0.0
c      rt(4,3)=a2
c
c----------------------------------------------------------------

c
c Calculation of the electrical field at the polarizable atoms
c
      j1=natom
      do 10 i1=1,Nsol
      do 10 k1=1,natsol
       j1=j1+1
c     
       xi(1)=r(j1,1)
       xi(2)=r(j1,2)
       xi(3)=r(j1,3)
c
c calls subroutine which calculates the field due to the quantum-part
c of the system.
       call efield(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
     >            c,a,RMM,xi,V,Ef)
c
       do 11 i2=1,i1-1
       do 11 k2=1,natsol
        j2=natom + natsol*(i2-1)+k2
c field due to the other solvent molecules charges
c this goes for the 4 point charges that represent the induced dipole
c and the permanent charge
c
       do 11 n1=1,4
c
       tx=r(j1,1)-(r(j2,1)+rt(n1,1))
       ty=r(j1,2)-(r(j2,2)+rt(n1,2))
       tz=r(j1,3)-(r(j2,3)+rt(n1,3))
       dd1=tx**2 + ty**2 + tz**2
       dd=sqrt(dd1)
       tv=pci(j2,n1)/dd
       V=V-tv
c
       te=tv/dd1
       Ef(1)=Ef(1)+ te*tx
       Ef(2)=Ef(2)+ te*ty
       Ef(3)=Ef(3)+ te*tz
c
 11   continue
c
       do 16 i2=i1+1,Nsol
       do 16 k2=1,natsol
        j2=natom + natsol*(i2-1)+k2
c field due to the other solvent molecules charges
c this goes for the 4 point charges that represent the induced dipole
c and the permanent charge
c
       do 16 n1=1,4
c
       tx=r(j1,1)-(r(j2,1)+rt(n1,1))
       ty=r(j1,2)-(r(j2,2)+rt(n1,2))
       tz=r(j1,3)-(r(j2,3)+rt(n1,3))
       dd1=tx**2 + ty**2 + tz**2
       dd=sqrt(dd1)
       tv=pci(j2,n1)/dd
       V=V-tv
c
       te=tv/dd1
       Ef(1)=Ef(1)+ te*tx
       Ef(2)=Ef(2)+ te*ty
       Ef(3)=Ef(3)+ te*tz
c
 16   continue
c
c-------------------------------------------------
c up to here it calculated the electrical field at site j1
       px(j1)=alpha(k1)*Ef(1)
       py(j1)=alpha(k1)*Ef(2)
       pz(j1)=alpha(k1)*Ef(3)
 10   continue
c
c CALCULATION OF NEW INDUCED CHARGES (following Sprik)
c
      Upol=0.0D0
c it's necessary to solve a system of 4 equations with four unknowns
      j1=natom
      do 19 i1=1,Nsol
      do 19 k1=1,natsol
       j1=j1+1
c
       qa=pc(k1)
       qx=px(j1)/a2
       qy=py(j1)/a2
       qz=pz(j1)/a2
c
       pci(j1,1)=(qa-qx-qy+qz)/4.0D0
       pci(j1,2)=(qa+qx+qy+qz)/4.0D0
       pci(j1,3)=(qa+qx-qy-qz)/4.0D0
       pci(j1,4)=(qa-qx+qy-qz)/4.0D0
c
c     px(j1)=pci(j1,1)*rt(1,1)+pci(j1,2)*rt(2,1)+
c    >  pci(j1,3)*rt(3,1)+pci(j1,4)*rt(4,1)
c
c     py(j1)=pci(j1,1)*rt(1,2)+pci(j1,2)*rt(2,2)+
c    >  pci(j1,3)*rt(3,2)+pci(j1,4)*rt(4,2)
c
c     pz(j1)=pci(j1,1)*rt(1,3)+pci(j1,2)*rt(2,3)+
c    >  pci(j1,3)*rt(3,3)+pci(j1,4)*rt(4,3)
c     write(*,*) 'INDUCED CALCULATED'
c     write(*,*) px(j1),py(j1),pz(j1)
c
c     pci(j1,1)=pc(k1)
c     pci(j1,2)=qx
c     pci(j1,3)=qy
c     pci(j1,4)=qz
c        do  i=1,4
c         write(*,*) pci(j1,i)
c        enddo
c
      if (alpha(k1).ne.0.0D0) then
       Upol=Upol+(px(j1)**2+py(j1)**2+pz(j1))**2/(2.0D0*alpha(k1))
      endif
      write(*,*) 'Upol=',Upol
 19   continue
c
      return
      end
