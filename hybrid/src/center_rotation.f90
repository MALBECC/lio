	subroutine center_rotation(atoms, mass, r, firstcet, In0vec)
!center system and keeps direction of autovectors of inercia tensor
!need more test
!2017 N. foglia
	implicit none
        integer, intent(in):: atoms
	double precision, intent(in), dimension(atoms) :: mass
	double precision, intent(inout), dimension(3,atoms) :: r
	double precision, dimension(3) :: RCM, rtemp
	double precision ::  Tmass
	double precision, dimension(3,3) :: Ini
	double precision, dimension(9) :: Inivec, Intempvec
	double precision, dimension(3,3) :: Inivec2
        double precision, dimension(9), intent(inout) :: In0vec
	double precision :: rfactor
	integer :: firstcet
        double precision, dimension(1000) :: work
        double precision, dimension(3) ::w
        integer, dimension(1000) :: iwork
        integer :: info
	integer :: LWORK, LIWORK, LWMAX

	double precision :: proj1, proj0
	integer :: i, j, j2
	integer :: case12, signcase, caseIn, signcaseIn
	integer :: wi

	LWMAX=1000


!calculo el centro de masa y centro al sistema
        RCM=0.d0
        Tmass=0.d0
        do i=1, atoms
          do j=1, 3
            RCM(j)=RCM(j)+mass(i)*r(j,i)
          end do
          Tmass=Tmass+mass(i)
        end do
        RCM=RCM/Tmass

        do i=1, atoms
          do j=1, 3
            r(j,i)=r(j,i)-RCM(j)
          end do
        end do

!calcula el tensor de inercia inicial
        Ini=0.d0
        do j=1,3
          do j2=1,3
            do i=1, atoms
               rfactor=-r(j,i)*r(j2,i)
               if (j.eq.j2) rfactor=rfactor + r(1,i)*r(1,i) + r(2,i)*r(2,i) + r(3,i)*r(3,i)
               Ini(j,j2)=Ini(j,j2)+mass(i)*rfactor
            end do
          end do
        end do

        wi=0
	Inivec2=0.d0
        do j=1,3
          do j2=j,3
            wi=wi+1
	    Inivec2(j,j2)=Ini(j,j2)
          end do
        end do

!calculo autovectores del tensor de inercia
        LWORK = -1
        LIWORK = -1
        call DSYEVD( 'V', 'U', 3, inivec2, 3, W, WORK, LWORK, IWORK, LIWORK, INFO )
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        LIWORK = MIN( LWMAX, IWORK( 1 ) )
        call DSYEVD( 'V', 'U', 3, inivec2, 3, W, WORK, LWORK, IWORK, LIWORK, INFO )
        if (info .ne. 0) stop "info ne 0, problem in dsyevd"

        wi=0
        do j=1,3
          do j2=1,3
            wi=wi+1
            inivec(wi)=Inivec2(j2,j)
          end do
        end do


        if (firstcet.eq.1) then
           in0vec=inivec
        else !ordeno ejes para q conincidan con el tensor de inercia del paso anterior
           proj0=0.d0
           do case12=1, 6
             do signcase=1, 8
               call switch(inivec, intempvec, case12, signcase)
               proj1=multi3vec(intempvec,in0vec)
               if (proj1.gt.proj0)  then
                 caseIn=case12
                 signcaseIn=signcase
                 proj0=proj1
               end if
             end do
           end do

          call switch(inivec, in0vec, caseIn, signcaseIn)
          inivec=in0vec
	end if

!paso el sistema a la base que diagonaliza el tensor de inercia
        do i=1, atoms
          rtemp=0.d0
          rtemp(1)=r(1,i)*inivec(1)+r(2,i)*inivec(2)+r(3,i)*inivec(3)
          rtemp(2)=r(1,i)*inivec(4)+r(2,i)*inivec(5)+r(3,i)*inivec(6)
          rtemp(3)=r(1,i)*inivec(7)+r(2,i)*inivec(8)+r(3,i)*inivec(9)
          do j=1,3
            r(j,i)=rtemp(j)
          end do
        end do

	write(966,*) atoms
        write(966,*)
        do i=1, atoms
          write(966,*) "1 ", r(1:3,i)
        end do

  345  format(2x, I2,    2x, 3(f10.6,2x))

	contains

	subroutine switch(VV1in, VV2, case12, signcase)
        implicit none
        double precision, dimension(9), intent(in) :: VV1in
	double precision, dimension(9) :: VV1
	double precision, dimension(9), intent(out) :: VV2
	integer, intent(in) :: case12, signcase

	if (signcase .eq. 1) then
	   VV1=VV1in
	elseif (signcase .eq. 2) then
	   VV1(1:3)=VV1in(1:3)
           VV1(4:6)=VV1in(4:6)
           VV1(7:9)=-VV1in(7:9)
        elseif (signcase .eq. 3) then
           VV1(1:3)=VV1in(1:3)
           VV1(4:6)=-VV1in(4:6)
           VV1(7:9)=VV1in(7:9)
        elseif (signcase .eq. 4) then
           VV1(1:3)=VV1in(1:3)
           VV1(4:6)=-VV1in(4:6)
           VV1(7:9)=-VV1in(7:9)
        elseif (signcase .eq. 5) then
           VV1(1:3)=-VV1in(1:3)
           VV1(4:6)=VV1in(4:6)
           VV1(7:9)=VV1in(7:9)
        elseif (signcase .eq. 6) then
           VV1(1:3)=-VV1in(1:3)
           VV1(4:6)=VV1in(4:6)
           VV1(7:9)=-VV1in(7:9)
        elseif (signcase .eq. 7) then
           VV1(1:3)=-VV1in(1:3)
           VV1(4:6)=-VV1in(4:6)
           VV1(7:9)=VV1in(7:9)
        elseif (signcase .eq. 8) then
           VV1(1:3)=-VV1in(1:3)
           VV1(4:6)=-VV1in(4:6)
           VV1(7:9)=-VV1in(7:9)
	else
           stop "bad case selected on sign"
        end if


	if (case12 .eq. 1) then
	   VV2=VV1
	elseif (case12 .eq. 2) then
	   VV2(1:3)=VV1(1:3)
	   VV2(4:6)=VV1(7:9)
	   VV2(7:9)=VV1(4:6)
        elseif (case12 .eq. 3) then
           VV2(1:3)=VV1(4:6)
           VV2(4:6)=VV1(1:3)
           VV2(7:9)=VV1(7:9)
        elseif (case12 .eq. 4) then
           VV2(1:3)=VV1(4:6)
           VV2(4:6)=VV1(7:9)
           VV2(7:9)=VV1(1:3)
        elseif (case12 .eq. 5) then
           VV2(1:3)=VV1(7:9)
           VV2(4:6)=VV1(1:3)
           VV2(7:9)=VV1(4:6)
        elseif (case12 .eq. 6) then
           VV2(1:3)=VV1(7:9)
           VV2(4:6)=VV1(4:6)
           VV2(7:9)=VV1(1:3)
	else
	   stop "bad case selected on switch"
	end if
	return
	end subroutine switch


	double precision function multi3vec(VV1, VV2)
	implicit none
	double precision, dimension(9), intent(in) :: VV1, VV2
	integer :: k
	multi3vec=0.d0
	do k=1, 9
	  multi3vec=multi3vec+VV1(k)*VV2(k)
	end do
	return
	end function multi3vec

	end subroutine center_rotation

