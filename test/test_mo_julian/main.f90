PROGRAM main
  
  USE mo_kind,   ONLY: i4, dp
  USE mo_julian, ONLY: caldat, julday, ndays, ndyin, date2dec, dec2date
  
  IMPLICIT NONE
  
  INTEGER(i4) :: dd, mm, yy
  INTEGER(i4) :: hh, ii, ss
  INTEGER(i4), DIMENSION(3) :: dd1, mm1, yy1, hh1, ii1, ss1
  REAL(dp), DIMENSION(3) :: dec1
  
  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_julian.f90'

  isgood = .true.

  ! julday, caldat
  if (julday(01,01,1900) /= 2415021) isgood = .false.
  call caldat(2415021, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.

  ! ndays, ndyin
  if (ndays(01,01,1900) /= 0) isgood = .false.
  call ndyin(0, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.
  if (ndays(1,1,-4713) /= -2415021) isgood = .false.
  call ndyin(-2415021, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= -4713)) isgood = .false.

  ! date2dec
  if (abs(date2dec(01,01,1900,12,0)-2415021.0_dp) > epsilon(1.0_dp)*2415021._dp*10.0_dp) isgood = .false.
  if (abs(date2dec(01,01,1900)-2415020.5_dp) > epsilon(1.0_dp)*2415020._dp*10.0_dp) isgood = .false.

  ! date2dec, dec2date - scalar
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

  ! date2dec, dec2date - vector
  call dec2date(date2dec((/01,10,27/),(/02,12,03/),(/-100,1,2000/),(/11,23,00/),(/11,57,47/),(/11,12,59/)), &
       dd1, mm1, yy1, hh1, ii1, ss1)
  if ((dd1(1) /= 1) .or. (mm1(1) /= 2) .or. (yy1(1) /= -100) .or. &
       (hh1(1) /= 11) .or. (ii1(1) /= 11) .or. (ss1(1) /= 11)) isgood = .false.
  if ((dd1(2) /= 10) .or. (mm1(2) /= 12) .or. (yy1(2) /= 1) .or. &
       (hh1(2) /= 23) .or. (ii1(2) /= 57) .or. (ss1(2) /= 12)) isgood = .false.
  if ((dd1(3) /= 27) .or. (mm1(3) /= 3) .or. (yy1(3) /= 2000) .or. &
       (hh1(3) /= 0) .or. (ii1(3) /= 47) .or. (ss1(3) /= 59)) isgood = .false.

  ! date2dec, dec2date - scalar/vector mix
  ! This does not work with NAG compiler, have to have temporary array
  ! call dec2date(date2dec(01,12,2000,(/11,12,13/),11,59), &
  !      dd1, mm1, yy1, hh1, ii1, ss1)
  dec1 = date2dec(01,12,2000,(/11,12,13/),11,59)
  call dec2date(dec1, dd1, mm1, yy1, hh1, ii1, ss1)
  if ((dd1(1) /= 1) .or. (mm1(1) /= 12) .or. (yy1(1) /= 2000) .or. &
       (hh1(1) /= 11) .or. (ii1(1) /= 11) .or. (ss1(1) /= 59)) isgood = .false.
  if ((dd1(2) /= 1) .or. (mm1(2) /= 12) .or. (yy1(2) /= 2000) .or. &
       (hh1(2) /= 12) .or. (ii1(2) /= 11) .or. (ss1(2) /= 59)) isgood = .false.
  if ((dd1(3) /= 1) .or. (mm1(3) /= 12) .or. (yy1(3) /= 2000) .or. &
       (hh1(3) /= 13) .or. (ii1(3) /= 11) .or. (ss1(3) /= 59)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_julian o.k.'
  else
     write(*,*) 'mo_julian failed!'
  endif

  ! ! Check 59 minutes, 60 seconds problem if eps too small in dec2date
  ! dec1(1) = 2443879.5_dp
  ! write(*,*) '#00: ', dec1(1), 1._dp/24._dp, date2dec(1,1,1990,1,0,0) - date2dec(1,1,1990,0,0,0)
  ! call dec2date(dec1(1), dd, mm, yy, hh, ii, ss)
  ! write(*,*) '#01: ', dec1(1), dd, mm, yy, hh, ii, ss
  ! dec1(1) = dec1(1) + 1._dp/24._dp
  ! call dec2date(dec1(1), dd, mm, yy, hh, ii, ss)
  ! write(*,*) '#02: ', dec1(1), dd, mm, yy, hh, ii, ss
  ! dec1(1) = dec1(1) - 1._dp/24._dp
  ! call dec2date(dec1(1), dd, mm, yy, hh, ii, ss)
  ! write(*,*) '#03: ', dec1(1), dd, mm, yy, hh, ii, ss
  ! dec1(1) = dec1(1) + date2dec(1,1,1990,1,0,0) - date2dec(1,1,1990,0,0,0)
  ! call dec2date(dec1(1), dd, mm, yy, hh, ii, ss)
  ! write(*,*) '#04: ', dec1(1), dd, mm, yy, hh, ii, ss
  
  ! ! Check date2dec with John
  ! write(*,*) '#05: ', date2dec(22,12,2005,15,00,00)
  ! write(*,*) '#06: ', date2dec(22,12,2005,16,00,00)
  ! write(*,*) '#07: ', date2dec(22,12,2005,17,00,00)
  ! write(*,*) '#08: ', date2dec(22,12,2005,17,01,00)
  ! write(*,*) '#09: ', date2dec(22,12,2005,17,00,01)

END PROGRAM main
