	subroutine intECP
	use ECP_mod, only : ecpmode, ecptypes, tipeECP, ZlistECP,nECP,bECP, aECP,Zcore, Lmax, expnumbersECP
	implicit none
	integer z,l,j,u,w,t

!        integer, dimension(118,0:5,10) :: nECP
!        double precision, dimension(118,0:5,10) :: bECP, aECP
!        integer, dimension(118,-2:5) :: dataECP
!	write(*,*) "rutina de calculo de integrales"
!	write(*,*) "testeando variables de entrada"
	do z=1,118
!	write(*,*) "carga", z
!	write(*,*) "lmax", dataECP(z,-1)
	do l=0,Lmax(Z)
!	write(*,*) "termminos", dataECP(z,l)
	do j=1,expnumbersECP(Z,l)
!	write(*,9018) z,l,j,nECP(z,l,j), bECP(z,l,j),aECP(z,l,j)
	end do
	end do
	end do
!	write(*,*) dataECP

!	do z=1,118
!		do u=1, dataECP(Z,dataECP(Z,-1))
!			write(*,*) nECP(Z,dataECP(Z,dataECP(Z,-1)),u), bECP(Z,dataECP(Z,dataECP(Z,-1)),u), aECP(Z,dataECP(Z,dataECP(Z,-1)),u)
!		end do
!		do w=0, dataECP(Z,-1)-1
!			do u=1, dataECP(Z,w)
!				write(*,*) nECP(Z,w,u), bECP(Z,w,u), aECP(Z,w,u)
!			end do
!		end do
!	end do

        do z=1,118
                do l=0, Lmax(Z)
                        do t=1, expnumbersECP(Z,l)
                                write(*,9018) Z,l,t,nECP(Z,l,t), bECP(Z,l,t), aECP(Z,l,t)
                        end do
                end do
        end do









	9018 format(/1x,'Z =',i4,2x, 'L =',i4,2x,'coefnumber =',i3,2x, 'n =',i2,2x, 'b =', f15.5,2x,'c =',f15.5)






	end subroutine intECP
