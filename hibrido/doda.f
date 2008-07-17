      SUBROUTINE DODA(NATSOL,R,NT,IZ)
C-----CAMBIA DE LA NUMERACION Y UNIDADES MIAS A LAS
C     DE DARIO==> UNIDADES ATOMICAS
      
      INCLUDE 'COMM'
      DIMENSION PCTE(NAT),XTE(NAT),YTE(NAT),ZTE(NAT)
      DIMENSION R(NT,3),IZ(NT),IZTE(NT)
      
c      do i=1+natom,nwat+natom
c      j=i+nwat
c      dist=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
c      dist=dsqrt(dist)/da(1)
c      write(*,*)'DODA 1',i,j,dist
c      enddo 

      DO I=1+NATOM,NPART
      PCTE(I)=PC(I)
      IZTE(I)=IZ(I)
      ENDDO

      DO I=1,NATOM
      R(I,1)=X(I)/A0
      R(I,2)=Y(I)/A0
      R(I,3)=Z(I)/A0
      ENDDO         
      

      NN=NATOM      
      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      R(N,1)=X(NN)/A0
      R(N,2)=Y(NN)/A0
      R(N,3)=Z(NN)/A0
      PC(N)=PCTE(NN)/EE
      IZ(N)=IZTE(NN)
      ENDDO

      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      R(N+1,1)=X(NN)/A0
      R(N+1,2)=Y(NN)/A0
      R(N+1,3)=Z(NN)/A0
      PC(N+1)=PCTE(NN)/EE
      IZ(N+1)=IZTE(NN)
      ENDDO   


      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      R(N+2,1)=X(NN)/A0
      R(N+2,2)=Y(NN)/A0
      R(N+2,3)=Z(NN)/A0
      PC(N+2)=PCTE(NN)/EE
      IZ(N+2)=IZTE(NN)
      ENDDO

c      do i=1+natom,npart,3
c      dist=(r(i,1)-r(i+1,1))**2+(r(i,2)-r(i+1,2))**2+
c     & (r(i,3)-r(i+1,3))**2
c      dist=dsqrt(dist)
c      write(*,*)'DODA 2',i,i+1,dist*0.52918/a0
c      enddo

c      DO I=NATOM+1,NPART
c      XTE(I)=X(I)
c      YTE(I)=Y(I)
c      ZTE(I)=Z(I)
c      PCTE(I)=PC(I)
c      ENDDO

c      DO II=1,NATSOL
c      II1=II-1
c      DO I=NATOM+1,NATOM+NWAT
c      K=(I-(NATOM+1))*(NATSOL-1)  
c      X(I+II1*NWAT)=XTE(I+II1+K)
c      Y(I+II1*NWAT)=YTE(I+II1+K)
c      Z(I+II1*NWAT)=ZTE(I+II1+K)
c      PC(I+II1*NWAT)=PCTE(I+II1+K)
c      write(*,*)' i  j ',i+ii1*nwat,i+ii1+k
c      ENDDO
c      ENDDO

c      do i=1+natom,npart
c      write(*,99)i,x(i),y(i),z(i)

c      DO I=1,NPART
c      R(I,1)=X(I)/A0
c      R(I,2)=Y(I)/A0
c      R(I,3)=Z(I)/A0
c      IF(I.GT.NATOM)THEN
c      PC(I)=PC(I)/EE
c      ENDIF
c      write(*,*)'doda ',i,pc(i),r(i,1),r(i,2),r(i,3)
c      ENDDO
      
    
        
       RETURN
       END
