Module mo_model

  use mo_kind   , only: dp, i4

  implicit none

  private

  public :: getrange
  public :: model

  interface model
     module procedure model_1d, model_0d
  end interface model

  ! module variables
  integer(i4) :: npara = 3

contains

  function GetRange()
    implicit none

    real(dp), dimension(npara,2)  :: GetRange   ! min value and max value of parameter
    real(dp), parameter           :: my_pi = 3.141592653589793238462643383279502884197_dp

    GetRange(1,:)  = (/  -my_pi,  my_pi /) ! min & max parameter 1
    GetRange(2,:)  = (/  -my_pi,  my_pi /) ! min & max parameter 2
    GetRange(3,:)  = (/  -my_pi,  my_pi /) ! min & max parameter 3

  end function GetRange

  function model_1d(paraset,x)
    implicit none

    real(dp), dimension(npara),   intent(in)  :: paraset   ! parameter set
    real(dp), dimension(:),       intent(in)  :: x         ! variables
    real(dp), dimension(size(x))              :: model_1d  ! modeled output

    ! local variables
    real(dp), parameter :: my_b = 2.0_dp

    model_1d(:) = sin(paraset(1))  + x(:)*(sin(paraset(2)))**2  + my_b * paraset(3)**4  * sin(paraset(1))

  end function model_1d

  function model_0d(paraset,x)
    implicit none

    real(dp), dimension(npara),   intent(in)  :: paraset   ! parameter set
    real(dp),                     intent(in)  :: x         ! variables
    real(dp)                                  :: model_0d  ! modeled output

    ! local variables
    real(dp), parameter :: my_b = 2.0_dp

    model_0d = sin(paraset(1))  + x*(sin(paraset(2)))**2  + my_b * paraset(3)**4  * sin(paraset(1))
  end function model_0d

end Module mo_model
