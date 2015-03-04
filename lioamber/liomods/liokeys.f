!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module liokeys
!--------------------------------------------------------------------!
       implicit none
       character(len=10) :: runtype

       namelist /lio_keywords/ runtype
       contains
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine start_liokeys(UnitSlot)
!--------------------------------------------------------------------!
       integer,intent(in) :: UnitSlot
       runtype='scf'
!
       read(unit=UnitSlot,nml=lio_keywords)
       end subroutine
       end module
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
