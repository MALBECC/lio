!         PROGRAM Fourier Transform

          implicit none

          double precision, allocatable :: mu(:)
          double precision ts, nu, fti, ftr, t, ene, pi, damp, lambda
          integer i, j, k
          integer ns, lmin, lmax

          character(20) input, output


          call get_command_argument(1,input)
          call get_command_argument(2,output)
          open (unit=1,file=input)
          open (unit=2,file=output)
          open (unit=3,file='espectro')
          open (unit=4, file='sq_espectro')
          open (unit=5, file='espectro_nm')

          write(*,*) 'steps, time-step (fs), lambda min, max (nm), damp'
          read(*,*) ns, ts, lmin, lmax, damp
 
          allocate (mu(1:ns))

          pi = 3.1415926535897932384626433832795d0
          ts = ts * 1.E-15 !conversion of time step to seconds
          damp = damp * 1.E15
          do i = 1, ns
            read(1,*) t, mu(i)
          enddo
          do i = 1, ns-1
             mu(i) = mu(i+1) - mu(i)  !takes differences
             mu(i)=mu(i)/ts
          enddo
          t=0
          do i=1,ns
             t=t+ts
             write(77,*) t, mu(i)
          enddo

! Calculation of spectrum in the region 20 - 1000 nm (1.5E+16 - 3E+14 Hz)

         do i = 1, 10*(lmax - lmin)

           lambda = dble((0.1*i) + lmin)
           nu = 3.0E17/lambda ! frequency computed in Hz
           
           fti = 0.0d0
           ftr = 0.0d0
           t=0.0

           do j = 1,ns-1
             t = t + ts
             ftr = ftr + cos(2*pi*t*nu) * mu(j)*exp(-t*damp) ! the real part of the transformation with a damping factor
             fti = fti + sin(2*pi*t*nu) * mu(j)*exp(-t*damp) ! the imaginary part of the transformation with a damping factor
           enddo

           ene = 1.241462463E3/lambda

           write (2,100) lambda, ene, ftr*ts, fti*ts
           write (3,101) ene, ftr*(-1)/(2*pi)
           write (5,101) lambda, ftr*(-1)/(2*pi)
!          write (5,101) lambda, fti*nu
         enddo

100       format (1x,I4,2x,E14.6,2x,E14.6,2x,E14.6)
101       format (1x,E14.6,2x,E14.6)
 
         END
