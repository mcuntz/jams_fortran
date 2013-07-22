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
  integer(i4) :: npara = 4

contains

  function GetRange()
    implicit none

    real(dp), dimension(npara,2)  :: GetRange   ! min value and max value of parameter

    GetRange(1,:)  = (/  -1.00_dp,  1.00_dp /) ! min & max parameter 1
    GetRange(2,:)  = (/  -1.00_dp,  1.00_dp /) ! min & max parameter 2
    GetRange(3,:)  = (/  -1.00_dp,  1.00_dp /) ! min & max parameter 3
    GetRange(4,:)  = (/  -1.00_dp,  1.00_dp /) ! min & max parameter 4

  end function GetRange

  function model_1d(paraset,x)
    implicit none

    real(dp), dimension(npara),   intent(in)  :: paraset   ! parameter set
    real(dp), dimension(:),       intent(in)  :: x         ! variables
    real(dp), dimension(size(x))              :: model_1d  ! modeled output

    model_1d(:) = paraset(1) * x(:)**3 + paraset(2) * x(:)**2 + paraset(3) * x(:) + paraset(4)
    ! model(:) = 0.0_dp*paraset(1) * x(:)**3 + paraset(2) * x(:)**2 + paraset(3) * x(:) + paraset(4)
  end function model_1d

  function model_0d(paraset,x)
    implicit none

    real(dp), dimension(npara),   intent(in)  :: paraset   ! parameter set
    real(dp),                     intent(in)  :: x         ! variables
    real(dp)                                  :: model_0d  ! modeled output

    model_0d = paraset(1) * x**3 + paraset(2) * x**2 + paraset(3) * x + paraset(4)
  end function model_0d

end Module mo_model
