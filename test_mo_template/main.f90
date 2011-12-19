PROGRAM main
  
  USE mo_kind,     ONLY: dp, sp
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
  isgood = isgood .and. (anint(10._dp*mean(dat)) == 55._dp)
  isgood = isgood .and. (anint(10._dp*mean(dat(1:9))) == 50._dp)
  isgood = isgood .and. (anint(10._dp*mean(dat, mask=(dat /= 10.))) == 50._dp)

  isgood = isgood .and. (anint(1e15_dp*PI_dp) == 3141592653589793._dp)

  if (isgood) then
     write(*,*) 'mo_template double precision o.k.'
  else
     write(*,*) 'mo_template double precision failed!'
  endif

  ! Double precision
  sat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  isgood = isgood .and. (anint(10._sp*mean(sat)) == 55._sp)
  isgood = isgood .and. (anint(10._sp*mean(sat(1:9))) == 50._sp)
  isgood = isgood .and. (anint(10._sp*mean(sat, mask=(sat /= 10.))) == 50._sp)

  isgood = isgood .and. (anint(1e6_sp*PI_sp) == 3141593._dp)

  if (isgood) then
     write(*,*) 'mo_template single precision o.k.'
  else
     write(*,*) 'mo_template single precision failed!'
  endif

END PROGRAM main
