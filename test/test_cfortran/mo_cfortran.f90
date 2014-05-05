module mo_cfortran

  ! This cannot be called mo_c because then mo_c.c and mo_c.f90 produce mo_c.o and overwrite each other

  use mo_kind, only: i4, dp

  implicit none

  public :: fctest

  private

contains

  function fctest(A, n1, n2, c)

    implicit none

    real(dp),    dimension(:,:), intent(in) :: A
    integer(i4),                 intent(in) :: n1, n2
    real(dp),                    intent(in) :: c
    real(dp)                                :: fctest

    real(dp) :: ctest
    external :: ctest

    ! Be aware that C is colunm-major so that one wants to transpose A, probably
    fctest = ctest(A, n1, n2, c)

  end function fctest

end module mo_cfortran
