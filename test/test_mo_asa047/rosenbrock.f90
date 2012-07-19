function rosenbrock(x)

!*****************************************************************************80
!
!! ROSENBROCK evaluates the Rosenbrock parabolic value function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, real(dp) X(2), the argument.
!
!    Output, real(dp) ROSENBROCK, the value of the function.
!
  use mo_kind,   only: i4, dp

  implicit none

  real(dp), dimension(:), intent(in) :: x
  real(dp) :: rosenbrock

  real(dp) :: fx
  real(dp) :: fx1
  real(dp) :: fx2

  fx1 = x(2) - x(1) * x(1)
  fx2 = 1.0_dp - x(1)

  fx = 100.0_dp * fx1 * fx1 + fx2 * fx2

  rosenbrock = fx

  return

end function rosenbrock
