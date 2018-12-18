!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine lio_finalize()
!--------------------------------------------------------------------!
! DEALLOCATION OF GLOBAL VARIABLES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod , only: dipole, Smat, RealRho, sqsm, Eorbs, 
     >                        Eorbs_b, MO_coef_at, MO_coef_at_b, r,
     >                        v, rqm, Em, Rm, pc, Iz, RMM, X, d
       use ECP_mod    , only: ecpmode
       use fileio_data, only: style
       use basis_subs , only: basis_deinit
 
       implicit none
       call basis_deinit()

       ! This should not be here.
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
       deallocate(pc, Iz, RMM, X,d)

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
