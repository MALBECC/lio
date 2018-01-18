c This subroutine calculates the solute-solvent LJ energy and forces
c----------------------------------------------------------------------------
      subroutine ljef(na_u,nac,natot,r,Em,Rm,f,Es,list)
      implicit none
      double precision
     .  r(3,natot), f(3,natot),
     .  Em(natot), Rm(natot), 
     .  rr(3),flj(3)
      integer list(nac)
c agregado, Nick
       double precision :: Elj, tx, ty, tz, dd, dd2
       double precision :: epsilon, sigma, t1, B, A, fej
       integer, intent(in) :: na_u, nac,natot
       double precision, intent(out) :: Es
       integer :: j2, j1

      Elj=0.0D0
      flj=0.0D0
      Es=0.d0
c----------------------------------------------------------------------------
c loop over MM atoms
      do j2=na_u+1,natot
c interaction only for ListQMMM solvent atoms
        if(list(j2-na_u).eq.0) then
c loop over QM atoms
          do  j1=1,na_u                
            tx=r(1,j1)-r(1,j2)
            ty=r(2,j1)-r(2,j2)
            tz=r(3,j1)-r(3,j2)
            rr(1)=tx
            rr(2)=ty
            rr(3)=tz
            dd2=tx**2d0 + ty**2d0 + tz**2d0
            dd=sqrt(dd2)
c Energy and forces from LJ term
            epsilon=sqrt(Em(j2)*Em(j1))
            sigma=0.50D0*(Rm(j2)+Rm(j1))
            t1=sigma**6d0
            B=4.0D0*epsilon*t1
            A=B*t1
            Elj = Elj+  A/dd**12d0 -B/dd**6d0
            fej = -12.0D0*A/dd**14d0 + 6.0D0*B/dd**8d0
            flj(1)=-2.d0*fej*tx
            flj(2)=-2.d0*fej*ty
            flj(3)=-2.d0*fej*tz

            f(1,j1)=f(1,j1) + flj(1)
            f(2,j1)=f(2,j1) + flj(2)
            f(3,j1)=f(3,j1) + flj(3)
            f(1,j2)=f(1,j2) - flj(1) 
            f(2,j2)=f(2,j2) - flj(2)
            f(3,j2)=f(3,j2) - flj(3)
          enddo !! at st
        endif !! list
      enddo !! at sv
      Es = 2.d0*Elj 
      return
      end subroutine ljef

