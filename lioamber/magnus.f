            subroutine magnus(F,rho,rho6,M,N)
! Entrada: Fock(t+(deltat/2)), rho(t)
! Salida:  rho6=rho(t+deltat)
! Ojo: todo entra y sale en base ortonormal!!
!       
       REAL*8 , intent(inout)  :: F(M,M)
       COMPLEX*8, intent(inout) :: rho(M,M),rho6(M,M)
       COMPLEX*8 :: rhonew(M,M),W(M,M)
       complex*8 :: x
       REAL*8 :: factorial(N)
       
       factorial(1)=1.D0
       do i=2,N
        factorial(i)=factorial(i-1)/i
       enddo

       x=cmplx(0.0d0,1.0d0)
! Magnus de orden 1: W = omega1 (JCTC 2011,7,1344-1355)
! F de entrada es F(t+(deltat/2))
       W=-x*F       
! Construimos el propagador de la densidad usando la formula de Baker-Campbell-Hausdorff (BCH) de orden 10 (buscar orden optimo)
       rho6=rho
       do i=1,N
!       rhonew=rho6+factorial(i)*rhonew  
       call conmutcc(W,rho,rhonew,M)
       rho6=rho6+factorial(i)*rhonew
       rho=rhonew
       enddo
       return
       end


























































