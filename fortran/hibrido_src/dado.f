      SUBROUTINE DADO(R,NATSOL,IZ)
C-----CAMBIA DE LA NUMERACION Y UNIDADES DE DARIO
C     A LAS MIAS

      use hibrido_common
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PCTE(NAT),XTE(NAT),YTE(NAT),ZTE(NAT)
      DIMENSION R(nt,3),IZ(NT),IZTE(NT)
      
      DO I=1+NATOM,NPART
      PCTE(I)=PC(I)
      IZTE(I)=IZ(I)
      ENDDO
      
      NN=NATOM
      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      X(NN)=R(N,1)*A0
      Y(NN)=R(N,2)*A0
      Z(NN)=R(N,3)*A0
      PC(NN)=PCTE(N)*EE
      IZ(NN)=IZTE(N)
      ENDDO

         
      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      X(NN)=R(N+1,1)*A0
      Y(NN)=R(N+1,2)*A0
      Z(NN)=R(N+1,3)*A0
      PC(NN)=PCTE(N+1)*EE
      IZ(NN)=IZTE(N+1)
      ENDDO

      DO N=NATOM+1,NATOM+NWAT*NATSOL,3
      NN=NN+1
      X(NN)=R(N+2,1)*A0
      Y(NN)=R(N+2,2)*A0
      Z(NN)=R(N+2,3)*A0
      PC(NN)=PCTE(N+2)*EE
      IZ(NN)=IZTE(N+2)
      ENDDO


c      DO I=1,NPART
c      X(I)=R(I,1)*A0
c      Y(I)=R(I,2)*A0
c      Z(I)=R(I,3)*A0
c      IF(I.GT.NATOM)THEN
c      PCTE(I)=PC(I)*EE
c      ENDIF

c      XTE(I)=X(I)
c      YTE(I)=Y(I)
c      ZTE(I)=Z(I)
c      ENDDO
     

c      DO II=1,NATSOL
c      II1=II-1
c      DO I=NATOM+1,NATOM+NWAT
c      K=(I-(NATOM+1))*(NATSOL-1)  
c      X(I+II1*NWAT)=XTE(I+II1+K)
c      Y(I+II1*NWAT)=YTE(I+II1+K)
c      Z(I+II1*NWAT)=ZTE(I+II1+K)
c      PC(I+II1*NWAT)=PCTE(I+II1+K)
c      ENDDO
c      ENDDO

      RETURN
      END
