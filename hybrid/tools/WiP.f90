!# Programa para analizar la contribucion atomica al trabajo en la coordenada de reaccion
!# V 1.0 Nick 2017
	Program WiP
	implicit none
	double precision, dimension(:), allocatable :: r, W1, W2
	double precision, dimension(:,:), allocatable :: W, Wcut
	integer :: natom, npic, natomlist
	integer :: i, j, l, trash
	double precision :: rmin, rmax
	double precision :: rminb, rmaxb
	integer, dimension(:), allocatable :: goodlist, catalist
	double precision :: E0


       character(20) input! file to analize
       character(4) rminc, rmaxc,npicc !input variables for lambda min, max and damping factor
       character(8) natomc
       LOGICAL :: existfile


       call get_command_argument(1,input)
       call get_command_argument(2,rminc)
       call get_command_argument(3,rmaxc)
       call get_command_argument(4,natomc)
       call get_command_argument(5,npicc)

       call logo()


	rminb=-1.2
	rmaxb=-0.2

!%%%%%%%%%%%%%%%%%%%%  Lectura de parametros
       if (len_trim(input) .eq. 0) then
        write(*,*) "WiP needs a Work_input, please use:"
        write(*,*) "./Wip Work_input"
        write(*,*) "or"
        write(*,*) "./Wip Work_input rmin rmax natoms npic "
        write(*,*)
        STOP
       end if

       INQUIRE(FILE=INPUT, EXIST=existfile)
       IF ( .NOT. existfile) THEN !verifica la existencia del archivo de parametros
        WRITE(*,*) input, "does not exist, please verify file name"
        STOP
       END IF

       if (len_trim(rminc).gt.0 .and. len_trim(rmaxc).gt.0 .and. len_trim(natomc) .gt. 0 .and. len_trim(npicc) .gt. 0) then
        read(rminc,*) rmin
        read(rmaxc,*) rmax
        read(natomc,*) natom
	read(npicc,*) npic
       else
        write(*,*)
        write(*,*)"please write:"
        write(*,*)'rmin rmax natom npic'
        read(*,*) rmin, rmax, natom, npic
        write(*,*)
       endif
!%%%%%%%%%%%%%%%%%%%%



	allocate(r(npic), W(npic,natom+2), goodlist(natom+2), catalist(natom+2))
	allocate(W1(npic), W2(npic))



	open(unit=20, file="todo.dat")
	do j=1, npic
	  read(20,*) r(j), trash, W(j, 1:natom+2)
	end do
	close(20)
	
	write(*,*) "r", r
!conversion Energia a kcal/mol y E=0 para el primer valor
        E0=W(1,natom+2)*23.06d0
	do j=1, npic
           W(j,natom+2)=W(j,natom+2)*23.06d0
!           if (r(j) .lt. rmin) W(j, 1:natom+2)=0.d0
           W(j,natom+2)=W(j,natom+2)-E0
!           if (r(j) .lt. rmin .or. r(j) .gt. rmax) W(j, 1:natom+2)=0.d0
        end do

!	do j=1, npic
!          write(*,*) r(j), W(j, 1:natom)
!        end do

!Analisis atomos relevantes
	W1=0.d0
!W1 W parte QM
!        do j=1, npic
!	  do i=1,24
!	    if (i.ne.5) then
!            if(i.ne.6) then
!            if(i.ne.9) then
!            if (i.ne.21) then !4 .and. i.ne.5 .and. i.ne. 8 .and. i.ne.20) then
!	      W1(j)=W1(j) + W(j,i)
!	    end if
!            end if
!            end if
!            end if
!
!          end do
!        end do


!Elimino atomos con contribucionMUY pequeña al W
        goodlist=0
	catalist=0
        natomlist=0

	do i=1, natom+2
	  do j=1, npic
	    if (abs(W(j,i)) .gt. 1D-10) then
              goodlist(i)=1
            end if

	    if (r(j) .lt. rmaxb .and. r(j) .gt. rminb) then
	      if (W(j,i) .lt. 0.d0)  then
	        catalist(i)=1
	      end if
	    end if

	  end do
	end do

	do i=1, natom+2
	  natomlist=natomlist+goodlist(i)
	end do

	write(*,*) natomlist, "atomos reelevantes"
	allocate(Wcut(npic,natomlist))



	open(unit=21, file="Analizado-Todo.dat")

	l=0
	do i=1, natom+2
	  if (goodlist(i).eq.1) then
	    if (l .lt. 10) then
	       write(21,996) l,i
	    elseif (l .lt. 100) then
	       write(21,997) l,i
	    elseif (l .lt. 1000) then
	       write(21,998) l,i
	    elseif (l .lt. 10000) then
	       write(21,999) l,i
	    elseif (l .lt. 100000) then
	       write(21,1000) l,i
	    else
	       stop "sistem have more than 100000 atoms"
	    end if
	    l=l+1
	  end if
	end do

	l=0
	do i=1, natom+2
	  if (goodlist(i).eq.1) then
	    if(l .lt. 10) then
               write(21,1101) l
            elseif (l .lt. 100) then
               write(21,1201) l
            elseif (l .lt. 1000) then
               write(21,1301) l
            elseif (l .lt. 10000) then
               write(21,1401) l
            elseif (l .lt. 100000) then
               write(21,1501) l
            else
               stop "sistem have more than 100000 atoms"
            end if

	    write(21,1002)

	    do j=1,npic
             if (r(j) .gt. rmin .and. r(j) .lt. rmax) then
 	      write(21,*) r(j),W(j,i)
	     end if
	    end do
	    write(21,1003)
	    l=l+1
	  end if
	end do
        close(21)





        open(unit=22, file="Analizado-cata.dat")

        l=0
        do i=1, natom+2
          if (catalist(i).eq.1) then
            if (l .lt. 10) then
               write(22,996) l,i
            elseif (l .lt. 100) then
               write(22,997) l,i
            elseif (l .lt. 1000) then
               write(22,998) l,i
            elseif (l .lt. 10000) then
               write(22,999) l,i
            elseif (l .lt. 100000) then
               write(22,1000) l,i
            else
               stop "sistem have more than 100000 atoms"
            end if
            l=l+1
          end if
        end do

        l=0
        do i=1, natom+2
          if (catalist(i).eq.1) then
            if(l .lt. 10) then
               write(22,1101) l
            elseif (l .lt. 100) then
               write(22,1201) l
            elseif (l .lt. 1000) then
               write(22,1301) l
            elseif (l .lt. 10000) then
               write(22,1401) l
            elseif (l .lt. 100000) then
               write(22,1501) l
            else
               stop "sistem have more than 100000 atoms"
            end if

            write(22,1002)

            do j=1,npic
             if (r(j) .gt. rmin .and. r(j) .lt. rmax) then
              write(22,*) r(j),W(j,i)
             end if
            end do
            write(22,1003)
            l=l+1
          end if
        end do
        close(22)







  996 format('@    s',i1,1x,'legend  "',i5,'"')
  997 format('@    s',i2,1x,'legend  "',i5,'"')
  998 format('@    s',i3,1x,'legend  "',i5,'"')
  999 format('@    s',i4,1x,'legend  "',i5,'"')
 1000 format('@    s',i5,1x,'legend  "',i5,'"')
 1101 format('@target G0.S',i1)
 1201 format('@target G0.S',i2)
 1301 format('@target G0.S',i3)
 1401 format('@target G0.S',i4)
 1501 format('@target G0.S',i5)
 1002 format('@type xy')
 1003 format('&')




 939 format(99999(2x,f18.9))

      contains

       subroutine logo()
         write(*,*)
         write(*,1200)
         write(*,1201)
         write(*,1202)
         write(*,1203)
         write(*,1204)
         write(*,1205)
	 write(*,1206)
	 write(*,1207)
	 write(*,1208)
	 write(*,1209)
	 write(*,1210)
	 write(*,1211)
	 write(*,1212)
         write(*,*)


 1200 FORMAT(4x,"█╗    ██╗ ██████╗ ██████╗ ██╗  ██╗     ██╗██╗███╗   ██╗██╗       ")
 1201 FORMAT(4x,"██║    ██║██╔═══██╗██╔══██╗██║ ██╔╝    ██╔╝██║████╗  ██║╚██╗      ")
 1202 FORMAT(4x,"██║ █╗ ██║██║   ██║██████╔╝█████╔╝     ██║ ██║██╔██╗ ██║ ██║      ")
 1203 FORMAT(4x,"██║███╗██║██║   ██║██╔══██╗██╔═██╗     ██║ ██║██║╚██╗██║ ██║      ")
 1204 FORMAT(4x,"╚███╔███╔╝╚██████╔╝██║  ██║██║  ██╗    ╚██╗██║██║ ╚████║██╔╝      ")
 1205 FORMAT(4x," ╚══╝╚══╝  ╚═════╝ ╚═╝  ╚═╝╚═╝  █╚═╝     ╚═╝╚═╝╚═╝  ╚═══╝╚═╝       ")
 1206 FORMAT(4x,"                                                                  ")
 1207 FORMAT(4x,"██████╗ ██████╗  ██████╗  ██████╗ ██████╗ ███████╗███████╗███████╗")
 1208 FORMAT(4x,"██╔══██╗██╔══██╗██╔═══██╗██╔════╝ ██╔══██╗██╔════╝██╔════╝██╔════╝")
 1209 FORMAT(4x,"██████╔╝██████╔╝██║   ██║██║  ███╗██████╔╝█████╗  ███████╗███████╗")
 1210 FORMAT(4x,"██╔═══╝ ██╔══██╗██║   ██║██║   ██║██╔══██╗██╔══╝  ╚════██║╚════██║")
 1211 FORMAT(4x,"██║     ██║  ██║╚██████╔╝╚██████╔╝██║  ██║███████╗███████║███████║")
 1212 FORMAT(4x,"╚═╝     ╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═╝╚══════╝╚══════╝╚══════╝")
       end subroutine logo
	end program
