PROGRAM main
  
  USE mo_kind,     ONLY: dp, sp, i8
  use mo_ansi_colors, only: color, c_red, c_green
  USE mo_template, ONLY: mean, PI_sp, PI_dp
  
  IMPLICIT NONE
  
  REAL(dp), DIMENSION(10) :: dat
  REAL(sp), DIMENSION(10) :: sat
  
  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_template.f90'

  ! Double precision
  dat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  isgood = isgood .and. (nint(10._dp*mean(dat)) == 55)
  isgood = isgood .and. (nint(10._dp*mean(dat(1:9))) == 50)
  isgood = isgood .and. (nint(10._dp*mean(dat, mask=( abs(dat-10._dp) .gt. tiny(1.0_dp) ))) == 50)

  isgood = isgood .and. (nint(1e15_dp*PI_dp,i8) == 3141592653589793_i8)

  if (isgood) then
     write(*,*) 'mo_template double precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_template double precision ', color('failed!', c_red)
  endif

  ! Double precision
  sat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  isgood = isgood .and. (nint(10._sp*mean(sat)) == 55)
  isgood = isgood .and. (nint(10._sp*mean(sat(1:9))) == 50)
  isgood = isgood .and. (nint(10._sp*mean(sat, mask=( abs(sat-10._sp) .gt. tiny(1.0_sp) ))) == 50)

  isgood = isgood .and. (nint(1e6_sp*PI_sp) == 3141593)

  if (isgood) then
     write(*,*) 'mo_template single precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_template single precision ', color('failed!', c_red)
  endif

END PROGRAM main
