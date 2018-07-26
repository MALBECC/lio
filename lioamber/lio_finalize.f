!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine lio_finalize()
!--------------------------------------------------------------------!
! DEALLOCATION OF GLOBAL VARIABLES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod
       use ECP_mod, only : ecpmode
       use fileio_data, only : style
       implicit none
       if (dipole) then
        if (style) write(69,8703)
        CLOSE(69)
       end if

!--------------------------------------------------------------------!
       if (allocated(Smat))      deallocate(Smat)
       if (allocated(RealRho))   deallocate(RealRho)
       if (allocated(sqsm))      deallocate(sqsm)
       if (allocated(Eorbs))     deallocate(Eorbs)
       if (allocated(Eorbs_b))      deallocate(Eorbs_b)
       if (allocated(MO_coef_at))   deallocate(MO_coef_at)
       if (allocated(MO_coef_at_b)) deallocate(MO_coef_at_b)

!--------------------------------------------------------------------!
       deallocate(r,v,rqm, Em, Rm)
       deallocate(pc, Iz, cx, ax, cd, ad, c, a)
      deallocate(Nuc,ncont,Nucx,ncontx,Nucd
     > ,ncontd, indexii, indexiid, RMM, X)
       deallocate(nnat, B, af)
       deallocate(natomc,nnps,nnpp,nnpd,nns)
       deallocate(nnd,nnp,atmin,jatc,d)

       call g2g_timer_summary()
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
