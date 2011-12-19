PROGRAM main
  
  USE mo_kind,   ONLY: i4
  USE mo_julian, ONLY: caldat, julday, ndays, ndyin
  
  IMPLICIT NONE
  
  INTEGER(i4) :: dd, mm, yy, julian
  
  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_julian.f90'

  isgood = .true.
  if (julday(01,01,1900) /= 2415021) isgood = .false.
  call caldat(2415021, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.
  if (ndays(01,01,1900) /= 0) isgood = .false.
  call ndyin(0, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.
  if (ndays(1,1,-4713) /= -2415021) isgood = .false.
  call ndyin(-2415021, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= -4713)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_julian o.k.'
  else
     write(*,*) 'mo_julian failed!'
  endif

END PROGRAM main
