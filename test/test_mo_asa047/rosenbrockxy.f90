function rosenbrockxy(p, xx, yy)

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
!    Input, real(dp) P(2), the argument.
!    Input, real(dp) X(1), the argument 2.
!    Input, real(dp) Y(1), the argument 3.
!
!    Output, real(dp) ROSENBROCK, the value of the function.
!
  use mo_kind,   only: dp

  implicit none

  real(dp), dimension(:), intent(in) :: p
  real(dp), dimension(:), intent(in) :: xx
  real(dp), dimension(:), intent(in) :: yy
  real(dp) :: rosenbrockxy

  real(dp) :: fx
  real(dp) :: fx1
  real(dp) :: fx2

  fx1 = p(2) - p(1) * p(1)
  fx2 = xx(1) - p(1)

  fx = yy(1) * fx1 * fx1 + fx2 * fx2

  rosenbrockxy = fx

end function rosenbrockxy
