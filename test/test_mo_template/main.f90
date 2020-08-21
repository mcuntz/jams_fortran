PROGRAM main

  USE mo_kind,     ONLY: dp, sp, i8
  use mo_ansi_colors, only: color, c_red, c_green
  USE mo_template, ONLY: mean, PI_sp, PI_dp

  IMPLICIT NONE

  REAL(dp), DIMENSION(10) :: dat
  REAL(sp), DIMENSION(10) :: sat

  real(sp) :: stmp
  real(dp) :: dtmp

  LOGICAL :: isgood

  Write(*,*) ''
  Write(*,*) 'Test mo_template.f90'

  ! Double precision
  dat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  dtmp = mean(dat)
  isgood = isgood .and. (nint(10._dp*dtmp) == 55)
  dtmp = mean(dat(1:9))
  isgood = isgood .and. (nint(10._dp*dtmp) == 50)
  dtmp = mean(dat, mask=( abs(dat-10._dp) .gt. tiny(1.0_dp) ))
  isgood = isgood .and. (nint(10._dp*dtmp) == 50)

  isgood = isgood .and. (nint(1e15_dp*PI_dp,i8) == 3141592653589793_i8)

  if (isgood) then
     write(*,*) 'mo_template double precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_template double precision ', color('failed!', c_red)
  endif

  ! Double precision
  sat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  stmp = mean(sat)
  isgood = isgood .and. (nint(10._sp*stmp) == 55)
  stmp = mean(sat(1:9))
  isgood = isgood .and. (nint(10._sp*stmp) == 50)
  stmp = mean(sat, mask=( abs(sat-10._sp) .gt. tiny(1.0_sp) ))
  isgood = isgood .and. (nint(10._sp*stmp) == 50)

  isgood = isgood .and. (nint(1e6_sp*PI_sp) == 3141593)

  if (isgood) then
     write(*,*) 'mo_template single precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_template single precision ', color('failed!', c_red)
  endif

END PROGRAM main
