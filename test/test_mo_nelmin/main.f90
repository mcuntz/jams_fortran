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
  use mo_nelmin, only: nelmin, nelminxy

  implicit none

  integer(i4), parameter :: n = 2

  integer(i4) :: icount
  integer(i4) :: ifault
  integer(i4) :: kcount
  integer(i4) :: konvge
  integer(i4) :: numres
  real(dp) :: reqmin
  INTERFACE
     FUNCTION rosenbrock(pp)
       USE mo_kind, ONLY: dp
       IMPLICIT NONE
       REAL(dp), DIMENSION(:), INTENT(IN) :: pp
       REAL(dp) :: rosenbrock
     END FUNCTION rosenbrock
  END INTERFACE
  INTERFACE
     FUNCTION rosenbrockxy(pp, xx, yy)
       USE mo_kind, ONLY: dp
       IMPLICIT NONE
       REAL(dp), DIMENSION(:), INTENT(IN) :: pp
       REAL(dp), DIMENSION(:), INTENT(IN) :: xx
       REAL(dp), DIMENSION(:), INTENT(IN) :: yy
       REAL(dp) :: rosenbrockxy
     END FUNCTION rosenbrockxy
  END INTERFACE
  real(dp) :: start(n)
  real(dp) :: step(n)
  real(dp) :: xmin(n)
  real(dp) :: ynewlo
  real(dp) :: xx(1)
  real(dp) :: yy(1)
  real(dp), allocatable :: history(:)
  integer(i4) :: i


  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_nelmin.f90'

  isgood = .True.
  start(1:n) = (/ -1.2_dp, 1.0_dp /)
  reqmin = 1.0E-08_dp
  step(1:n) = (/ 1.0_dp, 1.0_dp /)
  konvge = 10
  kcount = 500
  xmin = nelmin( rosenbrock, start, funcmin=ynewlo, varmin=reqmin, step=step, &
       konvge=konvge, maxeval=kcount, neval=icount, numrestart=numres, ierror=ifault, history=history)
  do i=50, icount, 50
     write(*,'(A15,I3,A24,F12.7)') ' Minimum after ', i, ' function evaluations:  ', history(i)
  end do
  write(*,'(A16,I3,A23,F12.7)')    ' Minimum found (',icount,' iterations):            ',ynewlo
  write(*,*) 'Call 1: ', xmin, anint(100._dp*xmin)
  isgood = isgood .and. (nint(100._dp*xmin(1)) == 100)
  isgood = isgood .and. (nint(100._dp*xmin(2)) == 100)

  start(1:n) = (/ -1.2_dp, 1.0_dp /)
  reqmin = 1.0E-08_dp
  step(1:n) = (/ 1.0_dp, 1.0_dp /)
  konvge = 10
  kcount = 500
  xmin = nelmin(rosenbrock, start, varmin=reqmin, step=step, konvge=konvge, maxeval=kcount)
  write(*,*) 'Call 2: ', xmin, anint(100._dp*xmin)
  isgood = isgood .and. (nint(100._dp*xmin(1)) == 100)
  isgood = isgood .and. (nint(100._dp*xmin(2)) == 100)

  start(1:n) = (/ -1.2_dp, 1.0_dp /)
  reqmin = 1.0E-08_dp
  konvge = 10
  kcount = 500
  xmin = nelmin(rosenbrock, start, varmin=reqmin, konvge=konvge, maxeval=kcount)
  write(*,*) 'Call 3: ', xmin, anint(100._dp*xmin)
  isgood = isgood .and. (nint(100._dp*xmin(1)) == 100)
  ! Less precision with new starting points
  isgood = isgood .and. (nint(100._dp*xmin(2)) == 99)

  start(1:n) = (/ -1.2_dp, 1.0_dp /)
  konvge = 10
  kcount = 500
  xmin = nelmin(rosenbrock, start, konvge=konvge, maxeval=kcount)
  write(*,*) 'Call 4: ', xmin, anint(100._dp*xmin)
  isgood = isgood .and. (nint(100._dp*xmin(1)) == 100)
  isgood = isgood .and. (nint(100._dp*xmin(2)) == 100)

  start(1:n) = (/ -1.2_dp, 1.0_dp /)
  xmin = nelmin(rosenbrock, start)
  write(*,*) 'Call 5: ', xmin, anint(100._dp*xmin)
  isgood = isgood .and. (nint(100._dp*xmin(1)) == 100)
  isgood = isgood .and. (nint(100._dp*xmin(2)) == 100)

  start(1:n) = (/ -12.0_dp, 10.0_dp /)
  xmin = nelmin(rosenbrock, start)
  write(*,*) 'Call 6: ', xmin, anint(100._dp*xmin)
  isgood = isgood .and. (nint(100._dp*xmin(1)) == 100)
  isgood = isgood .and. (nint(100._dp*xmin(2)) == 100)

  start(1:n) = (/ -12.0_dp, 10.0_dp /)
  xx(1) = 1.0_dp
  yy(1) = 100.0_dp
  xmin = nelminxy(rosenbrockxy, start, xx, yy, &
                neval=icount, history=history)
  do i=50, icount, 50
     write(*,'(A15,I3,A24,F12.7)') ' Minimum after ', i, ' function evaluations:  ', history(i)
  end do
  write(*,'(A16,I3,A23,F12.7)')           ' Minimum found (',icount,' iterations):            ',ynewlo
  write(*,*) 'Call 7: ', xmin, anint(100._dp*xmin)
  isgood = isgood .and. (nint(100._dp*xmin(1)) == 100)
  isgood = isgood .and. (nint(100._dp*xmin(2)) == 100)

  if (isgood) then
     write(*,*) 'mo_nelmin double precision o.k.'
  else
     write(*,*) 'mo_nelmin double precision failed!'
  endif


end program main
