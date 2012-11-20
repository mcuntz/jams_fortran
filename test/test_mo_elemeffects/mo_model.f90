Module mo_model

  use mo_kind   , only: dp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GetRange
  PUBLIC :: model

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

  function model(paraset,x)
    implicit none

    real(dp), dimension(npara),   intent(in)  :: paraset   ! parameter set
    real(dp), dimension(:),       intent(in)  :: x         ! variables
    real(dp), dimension(size(x))              :: model     ! modeled output

    model(:) = paraset(1) * x(:)**3 + paraset(2) * x(:)**2 + paraset(3) * x(:) + paraset(4)
  end function model

end Module mo_model
