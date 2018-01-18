      module fdf2

      private sp, dp

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)

c     Declarations for fdf procedures

      interface fdf_get
        module procedure fdf_int, fdf_dp, fdf_bool,
     $                   fdf_sp, fdf_str, fdf_phys
      end interface

      
      interface

         function fdf_defined(label)
         logical fdf_defined
         character(len=*), intent(in) :: label
         end function fdf_defined

         function fdf_enabled()
         logical fdf_enabled
         end function fdf_enabled

         function fdf_block(label,unit)
         logical fdf_block
         character(len=*), intent(in) :: label
         integer, intent(out)  :: unit
         end function fdf_block

         function fdf_convfac(unit1,unit2)
         real*8 fdf_convfac
         character(len=*), intent(in) :: unit1, unit2
         end function fdf_convfac

         subroutine fdf_init(filein,fileout)
         character(len=*), intent(in) :: filein, fileout
         end subroutine fdf_init

         subroutine fdf_inhibit
         end subroutine fdf_inhibit

      end interface


      contains

         function fdf_int(label,default)
         integer fdf_int
         character(len=*), intent(in) :: label
         integer, intent(in) ::  default
         integer fdf_integer
         external fdf_integer
         fdf_int = fdf_integer(label,default)
         end function fdf_int

         function fdf_dp(label,default)
         real(dp) fdf_dp
         character(len=*), intent(in) :: label
         real(dp), intent(in) ::  default
         real(dp) fdf_double
         external fdf_double
         fdf_dp = fdf_double(label,default)
         end function fdf_dp

         function fdf_sp(label,default)
         real(sp) fdf_sp
         character(len=*), intent(in) :: label
         real(sp), intent(in) ::  default
         real(sp) fdf_single
         external fdf_single
         fdf_sp = fdf_single(label,default)
         end function fdf_sp

         function fdf_phys(label,default,unit)
         real(dp) fdf_phys
         character(len=*), intent(in) :: label, unit
         real(dp), intent(in) ::  default
         real(dp) fdf_physical
         external fdf_physical
         fdf_phys = fdf_physical(label,default,unit)
         end function fdf_phys

         function fdf_bool(label,default)
         logical fdf_bool
         character(len=*), intent(in) :: label
         logical, intent(in) ::  default
         logical fdf_boolean
         external fdf_boolean
         fdf_bool = fdf_boolean(label,default)
         end function fdf_bool

         function fdf_str(label,default)
         character(len=80) fdf_str
         character(len=*), intent(in) :: label
         character(len=*), intent(in) ::  default
         character(len=80) fdf_string
         external fdf_string
         fdf_str =  fdf_string(label,default)
         end function fdf_str


      end module fdf2



