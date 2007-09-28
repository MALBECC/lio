c This subroutine calculates the solvent-solvent interaction energies
c To be used when the option sol$ (in namelist &int0) is selected.
c
c It needs coordinates for all atoms in all solvent molecules
c partial charges 
c Lennard-Jones parameters
c
c Dario Estrin, Buenos Aires February 1994.
c----------------------------------------------------------------------------
c
      subroutine mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
      implicit real*8 (a-h,o-z)
c
c Nsol : # solvent molecules
c natsol : # atoms in each solvent molecule
c pc : partial charges
c rs : nuclear coordinates 
c Es : solvent contributions to the energy
c
      INCLUDE 'param'
      parameter (pi=3.14159265358979312D0)
      dimension r(nt,3),Em(ntq+nss),Rm(ntq+nss)
      dimension Iz(natom),pc(nss),rr(3),rr2(3)
c
c
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
c     data RF,Ra,d0,theta0 /0.71701721D0,0.087635437D0,1.8088466D0,
c    >                       104.52D0/
c     data RF,Ra,d0,theta0 /0.20078511D0,2.669526D-05,1.8088466D0,
c    >                       104.52D0/
c datos 'reales'
c
      data RF,Ra,d0,theta0 /0.22D0,2.250000D-05,1.8088466D0,
     >                       104.52D0/
c datos para que sea rigida
c     data RF,Ra,d0,theta0 /100.0D0,100.0D0,1.8088466D0,
c    >                       104.52D0/
c kcal/mol and angstroms
c     C1=336.26649
c hartrees and bohr
      C1 = 1.0D0
c The Lennard-Jones parameters should be given in consistent units
c (kcal/mol angstroms if using the first set, or hartrees, bohr
c which is recommended).
c
c datos TIP4P ----------------
c corresponde a desplazar carga negativa 0.15A en direccion de los H
       alpha=0.7439762D0
       alpha2=0.1280119D0
c caso SPC
c       alpha=1.00D0
c       alpha2=0.00D0
c -----------------------------
c
c Lennard Jones parameters should be given in a.u (energy and distance)
c and are saved in a vector containing Em and Rm for solute atoms
c first and then for solvent atoms (Size = natom+natsol)
c
c
      Es1=0.0D0
      Es2=0.0D0
      Ei=0.0D0
      Elj=0.0D0
      Ese=0.0D0
c
      write(*,*) 'en mmsol'
c INTERMOLECULAR SOLVENT-SOLVENT TERMS ------------------------------------
c general
c
      j1=natom
      do 10 i1=1,Nsol
      do 10 k1=1,natsol
       j1=j1+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (k1.eq.1) then
        rr(1)=alpha*r(j1,1)+alpha2*(r(j1+1,1)+r(j1+2,1))
        rr(2)=alpha*r(j1,2)+alpha2*(r(j1+1,2)+r(j1+2,2))
        rr(3)=alpha*r(j1,3)+alpha2*(r(j1+1,3)+r(j1+2,3))
c
       else
        rr(1)=r(j1,1)
        rr(2)=r(j1,2)
        rr(3)=r(j1,3)
       endif
c       

       do 10 i2=i1+1,Nsol
       do 10 k2=1,natsol
        j2=natom + natsol*(i2-1)+k2
c
       if (k2.eq.1) then
        rr2(1)=alpha*r(j2,1)+alpha2*(r(j2+1,1)+r(j2+2,1))
        rr2(2)=alpha*r(j2,2)+alpha2*(r(j2+1,2)+r(j2+2,2))
        rr2(3)=alpha*r(j2,3)+alpha2*(r(j2+1,3)+r(j2+2,3))
c
       else
        rr2(1)=r(j2,1)
        rr2(2)=r(j2,2)
        rr2(3)=r(j2,3)
       endif
c       
c 
       dd=(r(j1,1)-r(j2,1))**2 + (r(j1,2)-r(j2,2))**2+ 
     >   (r(j1,3)-r(j2,3))**2 
       dd=sqrt(dd) 
c 
        dd2=(rr(1)-rr2(1))**2 + (rr(2)-rr2(2))**2+ 
     >   (rr(3)-rr2(3))**2 
        dd2=sqrt(dd2) 
     
c Lennard-Jones
        np1=k1+natom
        np2=k2+natom
        epsilon=sqrt(Em(np1)*Em(np2))
        sigma=0.50D0*(Rm(np1)+Rm(np2))
c
        t1=sigma**6
        B=4.0D0*epsilon*t1
        A=B*t1
c
c
c ELECTROSTATICAL TERMS
c considering all atoms with point charges
       Es1=Es1+pc(k1)*pc(k2)/dd2
       Es2=Es2+ A/dd**12 -B/dd**6
c
 10    continue
c 
c----------------------------------------------------------------------
c INTRAMOLECULAR SOLVENT TERMS --------------------------------
c for water
c
      if (natsol.eq.3) then
       do n=1,Nsol
        j=natom + natsol*(n-1)
c
        x1=-r(j+1,1)+r(j+2,1)
        y1=-r(j+1,2)+r(j+2,2)
        z1=-r(j+1,3)+r(j+2,3)
        d1=x1**2 + y1**2 + z1**2
        d1=sqrt(d1)
c
        x2=-r(j+1,1)+r(j+3,1)
        y2=-r(j+1,2)+r(j+3,2)
        z2=-r(j+1,3)+r(j+3,3)
        d2=x2**2 + y2**2 + z2**2
        d2=sqrt(d2)
c
        r1r2 = x1*x2 + y1*y2 + z1*z2
        t1=r1r2/(d1*d2)
        theta= acos(t1)/pi*180.0D0
       write(*,*) d1,d2,theta
c
c
        Ei=Ei+RF*(d1-d0)**2+RF*(d2-d0)**2+Ra*(theta-theta0)**2
       enddo

c      write(*,*) 'Ei =',Ei
       endif
c
c-----------------------------------------------------------------
c 
c SOLVENT-SOLUTE ELECTROSTATIC AND LJ  INTERACTIONS
c
       do 125 j1=1,natom
        j2=natom
       do 125 nc=1,Nsol
       do 125 ni=1,natsol
        j2=j2+1
c
c TIP4P case
c para el O desplazamiento del sitio
c
       if (ni.eq.1) then
        rr(1)=alpha*r(j2,1)+alpha2*(r(j2+1,1)+r(j2+2,1))
        rr(2)=alpha*r(j2,2)+alpha2*(r(j2+1,2)+r(j2+2,2))
        rr(3)=alpha*r(j2,3)+alpha2*(r(j2+1,3)+r(j2+2,3))
c
       else
        rr(1)=r(j2,1)
        rr(2)=r(j2,2)
        rr(3)=r(j2,3)
       endif
c       
c
       tx=r(j1,1)-r(j2,1)
       ty=r(j1,2)-r(j2,2)
       tz=r(j1,3)-r(j2,3)
       dd1=tx**2 + ty**2 + tz**2
       dd=sqrt(dd1)
c
       tx=r(j1,1)-rr(1)
       ty=r(j1,2)-rr(2)
       tz=r(j1,3)-rr(3)
       dd2=tx**2+ty**2+tz**2
       dd2=sqrt(dd2)
c
        Ese=Ese+Iz(j1)*pc(ni)/dd2
c
        npi=natom+ni
        epsilon=sqrt(Em(npi)*Em(j1))
        sigma=0.50D0*(Rm(npi)+Rm(j1))
c
        t1=sigma**6
        B=4.0D0*epsilon*t1
        A=B*t1
c
       Elj = Elj+ A/dd**12 -B/dd**6
 125   continue
c
       Ese=Ese*C1
c
c------------------------------------------------------------------
c
       Es=Elj + Ese + Es1 + Es2 + Ei 

       write(*,*) 'Elj,Ese',Elj,Ese
       write(*,*) 'Es1,Es2',Es1,Es2
       write(*,*) 'Ei=',Ei,'Es=',Es
       write(*,*) 'Es1+Es2=',Es1+Es2
       write(*,*) 'Es=',Es
       return
       end
