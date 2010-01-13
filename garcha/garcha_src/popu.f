c popu subroutine ----------------------------------
c calculates population in fragment orbitals in
c a larger system
c Dario Estrin, 1995
c in file frag, there should be Nfrag (number of atoms in
c fragment), and after that the ubication of the atoms in
c the input file, ordered by their appearence (given by the
c order of basis sets)
c---------------------------------------------------
       subroutine popu(NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >     Md,RMM,X,XX,ngdDyn,Nuc1,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
c
      implicit real*8 (a-h,o-z)
      logical NORM,ATRHO,VCINP,DIRECT,EXTR,dens,write
      logical OPEN,SVD,SHFT,GRAD,BSSE,integ,field,sol,free
      logical exter,contr
      integer nopt,iconst,igrid,igrid2
      INCLUDE 'param'
      dimension r(nt,3),nshell(0:3)
      dimension c(ng,nl),a(ng,nl),Nuc(M),Nuc1(M),ncont(M),Iz(nt)
c
c
      dimension indx(170),cof(170,170),rn(170),pop(170),kt(ng,ng)
      dimension cof0(170,170),N1(15),Nff(15),Nf1(14),Nord(14)
c
       dimension RMM(*),xi(3),q(nt)
c
c auxiliars
c X scratch space in matrices
      dimension XX(ngdDyn,ngdDyn)
c
c test
      dimension aux(ng),ds(nt)
c
      COMMON /TABLE/ STR(880,0:21)
c     common /HF/ nopt,OPEN,NMAX,NCO,ATRHO,VCINP,DIRECT,
c    >             IDAMP,EXTR,SHFT,SHI,GOLD,told,write,Nunp
c
      common /index/ index(ng)
c
c------------------------------------------------------------------
c
c Pointers
c

      sq2=sqrt(2.D0)
      MM=M*(M+1)/2 
      MM2=M**2
      MMd=Md*(Md+1)/2
      Md2=2*Md
      M2=2*M
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, F also uses the same position after S was used
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
      M13=M11+MM
c aux ( vector for ESSl)
      M15=M13+M
c Least squares
      M17=M15+MM
c vectors of MO
      M18=M17+MMd
c weights (in case of using option )
      M19=M18+M*NCO
c
      Nel=2*NCO+Nunp
c------------------------------------------------
c Initializations/Defaults
c
      write(*,*) ' FRAGMENT POPULATION CALCULATION  '
c
c
c H H core, 1 electron matrix elements
c
c
      call int1(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,c,a,RMM,En)
c
c calculates matrix of coefficients for new overlap kt(i,j)
cc
      do i=1,M
       Nuc1(index(i))=Nuc(i)
      enddo
c
       do i=1,M
        do j=1,M
         if (j.le.i) then
          k=i+(M2-j)*(j-1)/2
          kt(index(i),index(j))=M5+k-1
         else
          k=j+(M2-i)*(i-1)/2
          kt(index(i),index(j))=M5+k-1
         endif
        enddo
       enddo
c
       contr=.false.
       if (contr) then
       read(17,*)
c Calculation of contribution of atoms to each orbital
       do k=1,NCO+7
c
       read(17,*)
       read (17,*) (XX(l,k),l=1,M)
c
        do n=1,natom
         q(n)=0
        enddo
c
        do i=1,M
        do j=1,M
         t0=XX(i,k)*XX(j,k)*RMM(kt(i,j))
          q(Nuc1(i))=q(Nuc1(i))+t0
        enddo
        enddo
c
c
        write(*,*) 'CONTRIBUTIONS TO MO ',k
         write(*,770)
c
c
        do n=1,natom
         write(*,760) n,Iz(n),q(n)
        enddo
c
        enddo
c
c---------------------------------------------
       goto 999
       endif
c 
c----- Calculation of populations in MO of a fragment
c
c
      open(unit=16,file='frag')
c read coefficients of fragment MO
      read(16,*) Nfrag
      read(16,*) (Nord(n),n=1,Nfrag)
c
       Nf=0
c--------------------------
c 
      do i=1,M
       Nuc1(index(i))=Nuc(i)
      enddo
c
      do n=1,Nfrag
      Nff(n)=0
      do i=1,M
       if (Nuc1(i).eq.Nord(n)) then
        N1(n)=i
        Nff(n)=Nff(n)+1
       endif
      enddo
      N1(n)=N1(n)-Nff(n)+1
      enddo
c----------------------------
      do n=1,Nfrag
c      read(16,*) N1(n),Nff(n)
       Nf=Nf+Nff(n)
      enddo
c
      do n=1,Nfrag
       Nf1(n)=1
       do na=1,n-1
       Nf1(n)=Nf1(n)+Nff(na)
       enddo
      enddo
c
c
      do n=1,Nf
       read(16,*)
       pop(n)=0.0D0
       read (16,*) (cof(l,n),l=1,Nf)
      enddo
c
      do n=1,Nf
      do l=1,Nf
       cof0(l,n)=cof(l,n)
      enddo
      enddo
c
       call ludcmp(cof,Nf,170,indx,d)
c
c calculates matrix of coefficients for new overlap kt(i,j)
       do i=1,M
        do j=1,M
         if (j.le.i) then
          k=i+(M2-j)*(j-1)/2
          kt(index(i),index(j))=M5+k-1
         else
          k=j+(M2-i)*(i-1)/2
          kt(index(i),index(j))=M5+k-1
         endif
        enddo
       enddo
c  loop over occupied MO of total system
       do 155 kk=1,NCO
c
        k=0
       do n=1,Nfrag
        do i=1,Nff(n)
        k=k+1
        rn(k)=XX(i+N1(n)-1,kk)
        enddo
        enddo
c
c calls linear systems of equations
       call lubksb(cof,Nf,170,indx,rn)
c
c calculates population of iith MO of fragment
      N1(Nfrag+1)=M+1
      Nff(Nfrag+1)=0
      do 156 ii=1,Nf
c
c calculation with j not belonging to the fragment
       j=1
      do 101 na=1,Nfrag+1
      do 102 while (j.lt.N1(na))
c
      do nb=1,Nfrag
      do l=1,Nff(nb)
       pop(ii)=pop(ii)+rn(ii)*XX(j,kk)*cof0(Nf1(nb)+l-1,ii)*
     >          RMM(kt(N1(nb)+l-1,j))
      enddo
      enddo
       j=j+1
 102  continue
      j=j+Nff(na)
 101  continue
c
c calculation corresponding to j belonging to the fragment
      do 800 j=1,Nf
c      ss=0.0D0
      do 800 na=1,Nfrag
      do 800 nb=1,Nfrag
      do 800 l1=1,Nff(na)
      do 800 l2=1,Nff(nb)
       pop(ii)=pop(ii)+rn(ii)*rn(j)*cof0(Nf1(na)+l1-1,ii)*
     >      cof0(Nf1(nb)+l2-1,j)*RMM(kt(N1(na)+l1-1,N1(nb)+l2-1))
c      ss=ss+cof0(Nf1(na)+l1-1,ii)*cof0(Nf1(nb)+l2-1,ii)*
c    >  RMM(kt(N1(na)+l1-1,N1(nb)+l2-1))
 800  continue
c
c
 156  continue
 155  continue
c
      tot=0.0D0
      do jj=1,Nf
      write(*,400) jj,2.0D0*pop(jj)
      tot=tot+pop(jj)
      enddo
c
      write(*,599) 2.0D0*tot
c
 400  format(i3,3x,F7.4)
 599  format('Total population = ',F8.4)
 760  format(I3,9x,I3,6x,F10.4)
 770  format('ATOM #',4x,'ATOM TYPE',4x,'POPULATION')

 999  continue
      return
      end
