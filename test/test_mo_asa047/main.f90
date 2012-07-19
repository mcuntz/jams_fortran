program main

!*****************************************************************************80
!
!! Demonstrates the use of NELMIN on ROSENBROCK.
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
  use mo_kind,   only: i4, dp
  use mo_asa047, only: nelmin

  implicit none

  integer(i4), parameter :: n = 2

  integer(i4) :: i
  integer(i4) :: icount
  integer(i4) :: ifault
  integer(i4) :: kcount
  integer(i4) :: konvge
  integer(i4) :: numres
  real(dp) :: reqmin
  INTERFACE
     FUNCTION rosenbrock(xx)
       USE mo_kind, ONLY: dp
       IMPLICIT NONE
       REAL(dp), DIMENSION(:),           INTENT(IN) :: xx
       REAL(dp) :: rosenbrock
     END FUNCTION rosenbrock
  END INTERFACE
  !real(dp), external :: rosenbrock
  real(dp) :: start(n)
  real(dp) :: step(n)
  real(dp) :: xmin(n)
  real(dp) :: ynewlo

  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_asa04.f90'

  start(1:n) = (/ -1.2_dp, 1.0_dp /)
  reqmin = 1.0E-08_dp
  step(1:n) = (/ 1.0_dp, 1.0_dp /)
  konvge = 10
  kcount = 500
  call nelmin( rosenbrock, start, xmin, ynewlo, reqmin, step, &
    konvge, kcount, icount, numres, ifault )

  isgood = .True.
  isgood = isgood .and. (anint(1000._dp*xmin(1)) == 1001._dp)
  isgood = isgood .and. (anint(1000._dp*xmin(2)) == 1002._dp)

  if (isgood) then
     write(*,*) 'mo_asa047 double precision o.k.'
  else
     write(*,*) 'mo_asa047 double precision failed!'
  endif


end program main
