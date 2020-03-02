!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function commutator_dd(MA, MB) result(MC)
   implicit none
   real(kind=8), intent(in)  :: MA(:,:)
   real(kind=8), intent(in)  :: MB(:,:)
   real(kind=8), allocatable :: MC(:,:)
   real(kind=8), allocatable :: MP(:,:), MN(:,:)
   integer                   :: nn

   nn = size(MA, 1)
   allocate(MC(nn,nn), MP(nn,nn), MN(nn,nn))
   MP = matmul(MA, MB)
   MN = matmul(MB, MA)
   MC = MP - MN
end function

!--------------------------------------------------------------------!
function commutator_zd(MA, MB) result(MC)
   implicit none
   complex(kind=8), intent(in)  :: MA(:,:)
   real(kind=8)   , intent(in)  :: MB(:,:)
   complex(kind=8), allocatable :: MC(:,:)
   complex(kind=8), allocatable :: MP(:,:), MN(:,:)
   integer                      :: nn

   nn = size(MA, 1)
   allocate(MC(nn,nn), MP(nn,nn), MN(nn,nn))
   MP = matmul(MA, MB)
   MN = matmul(MB, MA)
   MC = MP - MN
end function

!--------------------------------------------------------------------!
function commutator_dz(MA, MB) result(MC)
   implicit none
   real(kind=8)   , intent(in)  :: MA(:,:)
   complex(kind=8), intent(in)  :: MB(:,:)
   complex(kind=8), allocatable :: MC(:,:)
   complex(kind=8), allocatable :: MP(:,:), MN(:,:)
   integer                      :: nn

   nn = size(MA, 1)
   allocate(MC(nn,nn), MP(nn,nn), MN(nn,nn))
   MP = matmul(MA, MB)
   MN = matmul(MB, MA)
   MC = MP - MN
end function

!--------------------------------------------------------------------!
function commutator_zz(MA, MB) result(MC)
   implicit none
   complex(kind=8), intent(in)  :: MA(:,:)
   complex(kind=8), intent(in)  :: MB(:,:)
   complex(kind=8), allocatable :: MC(:,:)
   complex(kind=8), allocatable :: MP(:,:), MN(:,:)
   integer                      :: nn

   nn = size(MA, 1)
   allocate(MC(nn,nn), MP(nn,nn), MN(nn,nn))
   MP = matmul(MA, MB)
   MN = matmul(MB, MA)
   MC = MP - MN
end function

!--------------------------------------------------------------------!
function commutator_cd(MA, MB) result(MC)
   implicit none
   complex(kind=4), intent(in)  :: MA(:,:)
   real(kind=8)   , intent(in)  :: MB(:,:)
   complex(kind=4), allocatable :: MC(:,:)
   complex(kind=4), allocatable :: MP(:,:), MN(:,:)
   integer                      :: nn

   nn = size(MA, 1)
   allocate(MC(nn,nn), MP(nn,nn), MN(nn,nn))
   MP = matmul(MA,real(MB,4))
   MN = matmul(real(MB,4), MA)
   MC = MP - MN
end function

!--------------------------------------------------------------------!
function commutator_dc(MA, MB) result(MC)
   implicit none
   real(kind=8)   , intent(in)  :: MA(:,:)
   complex(kind=4), intent(in)  :: MB(:,:)
   complex(kind=4), allocatable :: MC(:,:)
   complex(kind=4), allocatable :: MP(:,:), MN(:,:)
   integer                      :: nn

   nn = size(MA, 1)
   allocate(MC(nn,nn), MP(nn,nn), MN(nn,nn))
   MP = matmul(real(MA,4), MB)
   MN = matmul(MB,real(MA,4))
   MC = MP - MN
end function

!--------------------------------------------------------------------!
function commutator_cc(MA, MB) result(MC)
   implicit none
   complex(kind=4), intent(in)  :: MA(:,:)
   complex(kind=4), intent(in)  :: MB(:,:)
   complex(kind=4), allocatable :: MC(:,:)
   complex(kind=4), allocatable :: MP(:,:), MN(:,:)
   integer                      :: nn

   nn = size(MA, 1)
   allocate(MC(nn,nn), MP(nn,nn), MN(nn,nn))
   MP = matmul(MA, MB)
   MN = matmul(MB, MA)
   MC = MP - MN
end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
