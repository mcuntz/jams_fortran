PROGRAM main
  
  USE mo_kind, ONLY: dp, sp, i4
  USE mo_integrate,  ONLY: int_regular

  IMPLICIT NONE

  INTEGER(i4), PARAMETER  :: nn = 9
  REAL(dp), DIMENSION(nn) :: dat
  REAL(sp), DIMENSION(nn) :: sat
  INTEGER(i4)             :: ii
  
  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_integrate.f90'

  ! Double precision
  forall(ii=1:nn) dat(ii) = ii
  isgood = .true.
  isgood = isgood .and. (anint(int_regular(dat)) == 40._dp)

  if (isgood) then
     write(*,*) 'mo_integrate double precision o.k.'
  else
     write(*,*) 'mo_integrate double precision failed!'
  endif

  ! Double precision
  forall(ii=1:nn) sat(ii) = ii
  isgood = .true.
  isgood = isgood .and. (anint(int_regular(sat)) == 40._dp)

  if (isgood) then
     write(*,*) 'mo_integrate single precision o.k.'
  else
     write(*,*) 'mo_integrate single precision failed!'
  endif

END PROGRAM main
