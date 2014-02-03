PROGRAM main
  
  USE mo_kind,    ONLY: dp
  USE mo_boxcox,  ONLY: boxcox, invboxcox, get_boxcox

  IMPLICIT NONE
  
  LOGICAL :: isgood

  real(dp), dimension(20) :: dx  = (/ 4.86, 4.92, 3.37, 3.29, 4.43, 1.53, 3.76, 0.35, &
                                      1.73, 1.98, 4.26, 4.64, 2.88, 2.35, 8.35, 5.6 , &
                                      4.28,13.04, 8.94, 3.94 /)
  real(dp), dimension(20) :: dy, dyx
  real(dp)                :: dl

  Write(*,*) ''
  Write(*,*) 'Test mo_boxcox.f90'
  Write(*,*) ''
  Write(*,*) 'Lambda is  0.39050292804648679 in Python'
  Write(*,*) 'and        0.3889467523247853  in Fortran'

  ! Double precision
  Write(*,*) ''
  isgood = .true.
  dl = get_boxcox(dx)
  write(*,*) 'Dlambda: ', dl
  isgood = isgood .and. ((nint(1000._dp*dl)-389) == 0)
  dy  = boxcox(dx,dl)
  dyx = invboxcox(dy,dl)
  if (any(abs(dx-dyx) > 1000._dp*epsilon(1.0_dp))) isgood = .False.
  
  if (isgood) then
     write(*,*) 'mo_boxcox double precision o.k.'
  else
     write(*,*) 'mo_boxcox double precision failed!'
  endif

  ! Single precision
  Write(*,*) ''
  write(*,*) 'This example is not converging in single precision.'

END PROGRAM main
