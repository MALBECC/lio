      SUBROUTINE CONFIG(ANG,NATSOL,ntq,q,HISTO,HISTO2)
      INCLUDE 'COMM'
      dimension q(ntq),HISTO(100,100),HISTO2(100)

      

      IF(ICON.EQ.0) THEN
         IF(NSPECQ.EQ.0) CALL LATHWD(ANG,NATSOL)
          ITEL = -1
         if(ndft.ne.1)then    
c         write(*,*)'config'
         do i=1,natom
         pc(i)=pc1(i)
c         write(*,*)pc(i)/ee,pc1(i)/ee
         enddo     
         endif
          RETURN


      ELSE IF(ABS(ICON).EQ.1)THEN
         READ(8,*) IDUM,X99
         READ(8,*) ITEL,X,Y,Z,X0,Y0,Z0,VX,VY,VZ,VX0,VY0,VZ0
         READ(8,*) ITEL,ST,S0,SD,SD0
         READ(8,*) STQ,S0Q,SDQ,SD0Q
         READ(8,*) ITEL,PC,PC0,VQ
      
         DO I=1+NATOM,NATOM+NWAT
         PC(I) = -PC(I+NWAT)-PC(I+2*NWAT)
         ENDDO
    

         DO I = 1, IDUM
         QR= GAUSSN()
         ENDDO

         IDUM = -1

      ENDIF

      
      IF (ICON.LE.0)THEN
         ITEL = -1
      ELSEIF(ICON.EQ.1)THEN
         READ(8,*) ITEL,KG,KG1,HISTO,HISTO2,QK,AC
      ENDIF



      RETURN
      END




