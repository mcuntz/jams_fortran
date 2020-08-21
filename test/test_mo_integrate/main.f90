PROGRAM main

  USE mo_kind,        ONLY: dp, sp, i4
  USE mo_integrate,   ONLY: int_regular
  use mo_ansi_colors, only: color, c_red, c_green

  IMPLICIT NONE

  INTEGER(i4), PARAMETER  :: nn = 9
  REAL(dp), DIMENSION(nn) :: dat
  REAL(sp), DIMENSION(nn) :: sat
  INTEGER(i4)             :: ii
  real(sp)                :: stmp
  real(dp)                :: dtmp

  LOGICAL :: isgood

  Write(*,*) ''
  Write(*,*) 'Test mo_integrate.f90'

  ! Double precision
  forall(ii=1:nn) dat(ii) = ii
  isgood = .true.
  dtmp   = int_regular(dat)
  isgood = isgood .and. (nint(dtmp) == 40)

  if (isgood) then
     write(*,*) 'mo_integrate double precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_integrate double precision ', color('failed!', c_red)
  endif

  ! Double precision
  forall(ii=1:nn) sat(ii) = ii
  isgood = .true.
  stmp   = int_regular(sat)
  isgood = isgood .and. (nint(stmp) == 40)

  if (isgood) then
     write(*,*) 'mo_integrate single precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_integrate single precision ', color('failed!', c_red)
  endif

END PROGRAM main
