!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine lio_finalize()
!--------------------------------------------------------------------!
! DEALLOCATION OF GLOBAL VARIABLES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod
       use ECP_mod, only : ecpmode
       implicit none
       if (idip.eq.1) then
        write(69,8703)
        CLOSE(69)
       end if
!--------------------------------------------------------------------!
       if (allocated(Smat))    deallocate(Smat)
       if (allocated(RealRho)) deallocate(RealRho)
!--------------------------------------------------------------------!
       deallocate(r,v,rqm, Em, Rm, pc,Iz, nnat,
     > af,c,a,cx,ax,cd,ad,B,Nuc,ncont,Nucx,ncontx,Nucd
     > ,ncontd, indexii, indexiid, RMM, X, XX)
c       deallocate(old1,old2,old3)
       deallocate(natomc,nnps,nnpp,nnpd,nns)
       deallocate(nnd,nnp,atmin,jatc,d)
       call g2g_deinit()
       call aint_deinit()
!--------------------------------------------------------------------!
       if (ecpmode) call intECP(4) !desalocatea variables de pseudopotenciales


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Nuevos formatos, Nick
 8703 FORMAT(4x,"╚═══════════════╩",
     >"═══════════════╩═════",
     >"══════════╩══════════",
     >"═════╝")

       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
