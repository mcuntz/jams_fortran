PROGRAM main
  
  USE mo_kind,   ONLY: i4, dp
  USE mo_julian, ONLY: caldat, julday, ndays, ndyin, date2dec, dec2date
  
  IMPLICIT NONE
  
  INTEGER(i4) :: dd, mm, yy
  INTEGER(i4) :: hh, ii, ss
  real(dp)    :: dec
  
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

  if (abs(date2dec(01,01,1900,12,0)-2415021.0_dp) > epsilon(1.0_dp)*2415021._dp*10.0_dp) isgood = .false.
  if (abs(date2dec(01,01,1900)-2415020.5_dp) > epsilon(1.0_dp)*2415020._dp*10.0_dp) isgood = .false.

  call dec2date(date2dec(01,01,1900,12,01,02), dd, mm, yy, hh, ii, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 12) .or. (ii /= 01) .or. (ss /= 02)) isgood = .false.
  call dec2date(date2dec(01,01,1900), dd, mm, yy, hh, ii, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 00) .or. (ii /= 00) .or. (ss /= 00)) isgood = .false.
  call dec2date(date2dec(01,01,1900,23,59,59), dd, mm, yy, hh, ii, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 23) .or. (ii /= 59) .or. (ss /= 59)) isgood = .false.
  call dec2date(date2dec(01,01,01,12,12,12), dd, mm, yy, hh, ii, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 01) .or. (hh /= 12) .or. (ii /= 12) .or. (ss /= 12)) isgood = .false.
  call dec2date(date2dec(01,02,-100,11,11,11), dd, mm, yy, hh, ii, ss)
  if ((dd /= 1) .or. (mm /= 2) .or. (yy /= -100) .or. (hh /= 11) .or. (ii /= 11) .or. (ss /= 11)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_julian o.k.'
  else
     write(*,*) 'mo_julian failed!'
  endif

END PROGRAM main
