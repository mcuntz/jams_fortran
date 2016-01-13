PROGRAM main

  USE mo_kind,   ONLY: i4, dp
  USE mo_julian, ONLY: ndays, ndyin
  USE mo_julian, ONLY: caldat, julday, date2dec, dec2date, setCalendar

  IMPLICIT NONE

  INTEGER(i4) :: dd, mm, yy
  INTEGER(i4) :: hh, nn, ss
  INTEGER(i4) :: dd2, mm2, yy2
  INTEGER(i4) :: hh2, nn2, ss2
  REAL(dp)    :: ff, rr, gg
  INTEGER(i4), DIMENSION(3) :: dd1, mm1, yy1, hh1, nn1, ss1
  INTEGER(i4), DIMENSION(3) :: dd12, mm12, yy12, hh12, nn12, ss12
  REAL(dp),    DIMENSION(3) :: dec1
  REAL(dp),    DIMENSION(3) :: rr1, ff1, gg1
  INTEGER(i4) :: jj

  LOGICAL :: isgood

  Write(*,*) ''
  Write(*,*) 'Test mo_julian.f90'

  isgood = .true.

  
  ! -----------------------
  ! standard calendar tests
  
  ! julday, caldat
  do jj=1, 2524594 ! 01.01.2200
     call caldat(jj,dd,mm,yy)
     ss = julday(dd,mm,yy)
     if (jj /= ss) isgood = .false.
  end do
  if (julday(01,01,1900) /= 2415021) isgood = .false.
  call caldat(2415021, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.

  ! ndays, ndyin
  do jj=0, 109573 ! 01.01.2200
     call ndyin(jj,dd,mm,yy)
     ss = ndays(dd,mm,yy)
     if (jj /= ss) isgood = .false.
  end do
  if (ndays(01,01,1900) /= 0) isgood = .false.
  call ndyin(0, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.
  if (ndays(1,1,-4712) /= -2415021) isgood = .false.
  call ndyin(-2415021, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= -4712)) isgood = .false.

  ! date2dec
  if (abs(date2dec(01,01,1900,12,0)-2415021.0_dp) > epsilon(1.0_dp)*2415021._dp*10.0_dp) isgood = .false.
  if (abs(date2dec(01,01,1900)-2415020.5_dp) > epsilon(1.0_dp)*2415020._dp*10.0_dp) isgood = .false.

  ! date2dec, dec2date - scalar
  do jj=1, 2524594 ! 01.01.2200
     call random_number(ff)
     rr = real(jj) + ff
     call dec2date(rr,dd,mm,yy,hh,nn,ss)
     gg = date2dec(dd,mm,yy,hh,nn,ss)
     call dec2date(gg,dd2,mm2,yy2,hh2,nn2,ss2)
     if ((dd /= dd2) .or. (mm /= mm2) .or. (yy /= yy2) .or. (hh /= hh2) .or. (nn /= nn2) .or. (ss /= ss2)) isgood = .false.
  end do
  call dec2date(date2dec(01,01,1900,12,01,02), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 12) .or. (nn /= 01) .or. (ss /= 02)) isgood = .false.
  call dec2date(date2dec(01,01,1900), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 00) .or. (nn /= 00) .or. (ss /= 00)) isgood = .false.
  call dec2date(date2dec(01,01,1900,23,59,59), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 23) .or. (nn /= 59) .or. (ss /= 59)) isgood = .false.
  call dec2date(date2dec(01,01,01,12,12,12), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 01) .or. (hh /= 12) .or. (nn /= 12) .or. (ss /= 12)) isgood = .false.
  call dec2date(date2dec(01,02,-100,11,11,11), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 2) .or. (yy /= -100) .or. (hh /= 11) .or. (nn /= 11) .or. (ss /= 11)) isgood = .false.

  ! date2dec, dec2date - vector
  do jj=1, 2524594 ! 01.01.2200
     call random_number(ff1)
     rr1 = real(jj) + ff1
     call dec2date(rr1,dd1,mm1,yy1,hh1,nn1,ss1)
     gg1 = date2dec(dd1,mm1,yy1,hh1,nn1,ss1)
     call dec2date(gg1,dd12,mm12,yy12,hh12,nn12,ss12)
     if (any(dd1 /= dd12) .or. any(mm1 /= mm12) .or. any(yy1 /= yy12) .or. &
          any(hh1 /= hh12) .or. any(nn1 /= nn12) .or. any(ss1 /= ss12)) isgood = .false.
  end do
  call dec2date(date2dec((/01,10,27/),(/02,12,03/),(/-100,1,2000/),(/11,23,00/),(/11,57,47/),(/11,12,59/)), &
       dd1, mm1, yy1, hh1, nn1, ss1)
  if ((dd1(1) /= 1) .or. (mm1(1) /= 2) .or. (yy1(1) /= -100) .or. &
       (hh1(1) /= 11) .or. (nn1(1) /= 11) .or. (ss1(1) /= 11)) isgood = .false.
  if ((dd1(2) /= 10) .or. (mm1(2) /= 12) .or. (yy1(2) /= 1) .or. &
       (hh1(2) /= 23) .or. (nn1(2) /= 57) .or. (ss1(2) /= 12)) isgood = .false.
  if ((dd1(3) /= 27) .or. (mm1(3) /= 3) .or. (yy1(3) /= 2000) .or. &
       (hh1(3) /= 0) .or. (nn1(3) /= 47) .or. (ss1(3) /= 59)) isgood = .false.

  ! date2dec, dec2date - scalar/vector mix
  ! This does not work with NAG compiler, have to have temporary array
  ! call dec2date(date2dec(01,12,2000,(/11,12,13/),11,59), &
  !      dd1, mm1, yy1, hh1, nn1, ss1)
  do jj=1, 2524594 ! 01.01.2200
     call random_number(ff1)
     rr1 = real(jj) + ff1
     call dec2date(rr1,dd1,mm1,yy1,hh1,nn1,ss1)
     gg1 = date2dec(dd1(1),mm1(1),yy1(1),hh1,nn1,ss1)
     call dec2date(gg1,dd12,mm12,yy12,hh12,nn12,ss12)
     if (any(dd1(1) /= dd12) .or. any(mm1(1) /= mm12) .or. any(yy1(1) /= yy12) .or. &
          any(hh1 /= hh12) .or. any(nn1 /= nn12) .or. any(ss1 /= ss12)) isgood = .false.
  end do
  dec1 = date2dec(01,12,2000,(/11,12,13/),11,59)
  call dec2date(dec1, dd1, mm1, yy1, hh1, nn1, ss1)
  if ((dd1(1) /= 1) .or. (mm1(1) /= 12) .or. (yy1(1) /= 2000) .or. &
       (hh1(1) /= 11) .or. (nn1(1) /= 11) .or. (ss1(1) /= 59)) isgood = .false.
  if ((dd1(2) /= 1) .or. (mm1(2) /= 12) .or. (yy1(2) /= 2000) .or. &
       (hh1(2) /= 12) .or. (nn1(2) /= 11) .or. (ss1(2) /= 59)) isgood = .false.
  if ((dd1(3) /= 1) .or. (mm1(3) /= 12) .or. (yy1(3) /= 2000) .or. &
       (hh1(3) /= 13) .or. (nn1(3) /= 11) .or. (ss1(3) /= 59)) isgood = .false.

  ! Mix
  do jj=1, 2524594 ! 01.01.2200
     call caldat(jj,dd,mm,yy)
     ss = int(date2dec(dd,mm,yy,12,0,0), i4)
     if (jj /= ss) isgood = .false.
  end do
  
  ! ----------------------
  ! 360days calendar tests

  ! julday360, caldat360
  do jj=0, 792000 ! 01.01.2200
     call caldat(jj, dd, mm, yy, calendar=3)
     ss = julday(dd, mm, yy, calendar=3)
     if (jj /= ss) isgood = .false.
  end do
  if (julday(01, 01, 1900, calendar=3) /= 684000) isgood = .false.
  call caldat(684000, dd, mm, yy, calendar=3)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.

  call setCalendar("360day")

  ! date2dec360
  if (abs(date2dec(01,01,1900,12,0)-684000.0_dp) > epsilon(1.0_dp)*684000.0_dp*10.0_dp) isgood = .false.
  if (abs(date2dec(01,01,1900)-683999.5_dp) > epsilon(1.0_dp)*683999._dp*10.0_dp) isgood = .false.

  ! date2dec360, dec2date360 - scalar
  do jj=1, 109573 ! 01.01.2200
     call random_number(ff)
     rr = real(jj) + ff
     call dec2date(rr,dd,mm,yy,hh,nn,ss)
     gg = date2dec(dd,mm,yy,hh,nn,ss)
     call dec2date(gg,dd2,mm2,yy2,hh2,nn2,ss2)
     if ((dd /= dd2) .or. (mm /= mm2) .or. (yy /= yy2) .or. (hh /= hh2) .or. (nn /= nn2) .or. (ss /= ss2)) isgood = .false.
  end do

  call dec2date(date2dec(01,01,1900,12,01,02), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 12) .or. (nn /= 01) .or. (ss /= 02)) isgood = .false.
  call dec2date(date2dec(01,01,1900), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 00) .or. (nn /= 00) .or. (ss /= 00)) isgood = .false.
  call dec2date(date2dec(01,01,1900,23,59,59), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 23) .or. (nn /= 59) .or. (ss /= 59)) isgood = .false.
  call dec2date(date2dec(01,01,01,12,12,12), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 01) .or. (hh /= 12) .or. (nn /= 12) .or. (ss /= 12)) isgood = .false.
  call dec2date(date2dec(01,02,-100,11,11,11), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 2) .or. (yy /= -100) .or. (hh /= 11) .or. (nn /= 11) .or. (ss /= 11)) isgood = .false.

  ! date2dec360, dec2date360 - vector
  do jj=1, 2524594 ! 01.01.2200
     call random_number(ff1)
     rr1 = real(jj) + ff1
     call dec2date(rr1,dd1,mm1,yy1,hh1,nn1,ss1)
     gg1 = date2dec(dd1,mm1,yy1,hh1,nn1,ss1)
     call dec2date(gg1,dd12,mm12,yy12,hh12,nn12,ss12)
     if (any(dd1 /= dd12) .or. any(mm1 /= mm12) .or. any(yy1 /= yy12) .or. &
          any(hh1 /= hh12) .or. any(nn1 /= nn12) .or. any(ss1 /= ss12)) isgood = .false.
  end do
  call dec2date(date2dec((/01,10,27/),(/02,12,03/),(/-100,1,2000/),(/11,23,00/),(/11,57,47/),(/11,12,59/)), &
       dd1, mm1, yy1, hh1, nn1, ss1)
  if ((dd1(1) /= 1) .or. (mm1(1) /= 2) .or. (yy1(1) /= -100) .or. &
       (hh1(1) /= 11) .or. (nn1(1) /= 11) .or. (ss1(1) /= 11)) isgood = .false.
  if ((dd1(2) /= 10) .or. (mm1(2) /= 12) .or. (yy1(2) /= 1) .or. &
       (hh1(2) /= 23) .or. (nn1(2) /= 57) .or. (ss1(2) /= 12)) isgood = .false.
  if ((dd1(3) /= 27) .or. (mm1(3) /= 3) .or. (yy1(3) /= 2000) .or. &
       (hh1(3) /= 0) .or. (nn1(3) /= 47) .or. (ss1(3) /= 59)) isgood = .false.

  ! date2dec360, dec2date360 - scalar/vector mix
  ! This does not work with NAG compiler, have to have temporary array
  ! call dec2date360(date2dec360(01,12,2000,(/11,12,13/),11,59), &
  !      dd1, mm1, yy1, hh1, nn1, ss1)
  do jj=1, 2524594 ! 01.01.2200
     call random_number(ff1)
     rr1 = real(jj) + ff1
     call dec2date(rr1,dd1,mm1,yy1,hh1,nn1,ss1)
     gg1 = date2dec(dd1(1),mm1(1),yy1(1),hh1,nn1,ss1)
     call dec2date(gg1,dd12,mm12,yy12,hh12,nn12,ss12)
     if (any(dd1(1) /= dd12) .or. any(mm1(1) /= mm12) .or. any(yy1(1) /= yy12) .or. &
          any(hh1 /= hh12) .or. any(nn1 /= nn12) .or. any(ss1 /= ss12)) isgood = .false.
  end do
  dec1 = date2dec(01,12,2000,(/11,12,13/),11,59)
  call dec2date(dec1, dd1, mm1, yy1, hh1, nn1, ss1)
  if ((dd1(1) /= 1) .or. (mm1(1) /= 12) .or. (yy1(1) /= 2000) .or. &
       (hh1(1) /= 11) .or. (nn1(1) /= 11) .or. (ss1(1) /= 59)) isgood = .false.
  if ((dd1(2) /= 1) .or. (mm1(2) /= 12) .or. (yy1(2) /= 2000) .or. &
       (hh1(2) /= 12) .or. (nn1(2) /= 11) .or. (ss1(2) /= 59)) isgood = .false.
  if ((dd1(3) /= 1) .or. (mm1(3) /= 12) .or. (yy1(3) /= 2000) .or. &
       (hh1(3) /= 13) .or. (nn1(3) /= 11) .or. (ss1(3) /= 59)) isgood = .false.

  ! Mix
  do jj=1, 79200 ! 01.01.2200
     call caldat(jj,dd,mm,yy)
     ss = int(date2dec(dd,mm,yy,12,0,0), i4)
     if (jj /= ss) isgood = .false.
  end do

  ! ---------------------
  ! 365day calendar tests

  call setCalendar("365day")
  
  ! julday365, caldat365
  do jj=0, 803000 ! 01.01.2200
     call caldat(jj,dd,mm,yy)
     ss = julday(dd,mm,yy)
     if (jj /= ss) isgood = .false.
  end do
  if (julday(01,01,1900) /= 693500) isgood = .false.
  call caldat(693500, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.


  ! date2dec365
  if (abs(date2dec(01,01,1900,12,0)-693500.0_dp) > epsilon(1.0_dp)*693500.0_dp*10.0_dp) isgood = .false.
  if (abs(date2dec(01,01,1900)-693499.5_dp) > epsilon(1.0_dp)*693499._dp*10.0_dp) isgood = .false.

  ! date2dec365, dec2date365 - scalar
  do jj=1, 803000 ! 01.01.2200
     call random_number(ff)
     rr = real(jj) + ff
     call dec2date(rr,dd,mm,yy,hh,nn,ss)
     gg = date2dec(dd,mm,yy,hh,nn,ss)
     call dec2date(gg,dd2,mm2,yy2,hh2,nn2,ss2)
     if ((dd /= dd2) .or. (mm /= mm2) .or. (yy /= yy2) .or. (hh /= hh2) .or. (nn /= nn2) .or. (ss /= ss2)) isgood = .false.
  end do

  call dec2date(date2dec(01,01,1900,12,01,02), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 12) .or. (nn /= 01) .or. (ss /= 02)) isgood = .false.
  call dec2date(date2dec(01,01,1900), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 00) .or. (nn /= 00) .or. (ss /= 00)) isgood = .false.
  call dec2date(date2dec(01,01,1900,23,59,59), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 23) .or. (nn /= 59) .or. (ss /= 59)) isgood = .false.
  call dec2date(date2dec(01,01,01,12,12,12), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 01) .or. (hh /= 12) .or. (nn /= 12) .or. (ss /= 12)) isgood = .false.
  call dec2date(date2dec(01,02,-100,11,11,11), dd, mm, yy, hh, nn, ss)
  if ((dd /= 1) .or. (mm /= 2) .or. (yy /= -100) .or. (hh /= 11) .or. (nn /= 11) .or. (ss /= 11)) isgood = .false.

  call setCalendar("julian")
  
  ! date2dec365, dec2date365 - vector
  do jj=1, 803000 ! 01.01.2200
     call random_number(ff1)
     rr1 = real(jj) + ff1
     call dec2date(rr1, dd1, mm1, yy1, hh1, nn1, ss1, calendar=2)
     gg1 = date2dec(dd1, mm1, yy1, hh1, nn1, ss1, calendar=2)
     call dec2date(gg1, dd12, mm12, yy12, hh12, nn12, ss12, calendar=2)
     if (any(dd1 /= dd12) .or. any(mm1 /= mm12) .or. any(yy1 /= yy12) .or. &
          any(hh1 /= hh12) .or. any(nn1 /= nn12) .or. any(ss1 /= ss12)) isgood = .false.
  end do

  call dec2date(date2dec((/01,10,27/),(/02,12,03/),(/-100,1,2000/),(/11,23,00/),(/11,57,47/),(/11,12,59/), calendar=2), &
       dd1, mm1, yy1, hh1, nn1, ss1, calendar=2)
  if ((dd1(1) /= 1) .or. (mm1(1) /= 2) .or. (yy1(1) /= -100) .or. &
       (hh1(1) /= 11) .or. (nn1(1) /= 11) .or. (ss1(1) /= 11)) isgood = .false.
  if ((dd1(2) /= 10) .or. (mm1(2) /= 12) .or. (yy1(2) /= 1) .or. &
       (hh1(2) /= 23) .or. (nn1(2) /= 57) .or. (ss1(2) /= 12)) isgood = .false.
  if ((dd1(3) /= 27) .or. (mm1(3) /= 3) .or. (yy1(3) /= 2000) .or. &
       (hh1(3) /= 0) .or. (nn1(3) /= 47) .or. (ss1(3) /= 59)) isgood = .false.


  ! date2dec365, dec2date365 - scalar/vector mix
  ! This does not work with NAG compiler, have to have temporary array
  ! call dec2date365(date2dec365(01,12,2000,(/11,12,13/),11,59), &
  !      dd1, mm1, yy1, hh1, nn1, ss1)
  do jj=1, 803000 ! 01.01.2200
     call random_number(ff1)
     rr1 = real(jj) + ff1
     call dec2date(rr1, dd1, mm1, yy1, hh1, nn1, ss1, calendar=2)
     gg1 = date2dec(dd1(1), mm1(1), yy1(1), hh1, nn1, ss1, calendar=2)
     call dec2date(gg1, dd12, mm12, yy12, hh12, nn12, ss12, calendar=2)
     if (any(dd1(1) /= dd12) .or. any(mm1(1) /= mm12) .or. any(yy1(1) /= yy12) .or. &
          any(hh1 /= hh12) .or. any(nn1 /= nn12) .or. any(ss1 /= ss12)) isgood = .false.
  end do
  dec1 = date2dec(01, 12, 2000, (/11,12,13/), 11, 59, calendar=2)
  call dec2date(dec1, dd1, mm1, yy1, hh1, nn1, ss1, calendar=2)
  if ((dd1(1) /= 1) .or. (mm1(1) /= 12) .or. (yy1(1) /= 2000) .or. &
       (hh1(1) /= 11) .or. (nn1(1) /= 11) .or. (ss1(1) /= 59)) isgood = .false.
  if ((dd1(2) /= 1) .or. (mm1(2) /= 12) .or. (yy1(2) /= 2000) .or. &
       (hh1(2) /= 12) .or. (nn1(2) /= 11) .or. (ss1(2) /= 59)) isgood = .false.
  if ((dd1(3) /= 1) .or. (mm1(3) /= 12) .or. (yy1(3) /= 2000) .or. &
       (hh1(3) /= 13) .or. (nn1(3) /= 11) .or. (ss1(3) /= 59)) isgood = .false.

  ! Mix
  do jj=1, 803000 ! 01.01.2200
     call caldat(jj, dd, mm, yy, calendar=2)
     ss = int(date2dec(dd, mm, yy, 12, 0, 0, calendar=2), i4)
     if (jj /= ss) isgood = .false.
  end do

  if (isgood) then
     write(*,*) 'mo_julian o.k.'
  else
     write(*,*) 'mo_julian failed!'
  endif

  ! !  ! ! Check 59 minutes, 60 seconds problem if eps too small in dec2date
  ! !  ! dec1(1) = 2443879.5_dp
  ! !  ! write(*,*) '#00: ', dec1(1), 1._dp/24._dp, date2dec(1,1,1990,1,0,0) - date2dec(1,1,1990,0,0,0)
  ! !  ! call dec2date(dec1(1), dd, mm, yy, hh, nn, ss)
  ! !  ! write(*,*) '#01: ', dec1(1), dd, mm, yy, hh, nn, ss
  ! !  ! dec1(1) = dec1(1) + 1._dp/24._dp
  ! !  ! call dec2date(dec1(1), dd, mm, yy, hh, nn, ss)
  ! !  ! write(*,*) '#02: ', dec1(1), dd, mm, yy, hh, nn, ss
  ! !  ! dec1(1) = dec1(1) - 1._dp/24._dp
  ! !  ! call dec2date(dec1(1), dd, mm, yy, hh, nn, ss)
  ! !  ! write(*,*) '#03: ', dec1(1), dd, mm, yy, hh, nn, ss
  ! !  ! dec1(1) = dec1(1) + date2dec(1,1,1990,1,0,0) - date2dec(1,1,1990,0,0,0)
  ! !  ! call dec2date(dec1(1), dd, mm, yy, hh, nn, ss)
  ! !  ! write(*,*) '#04: ', dec1(1), dd, mm, yy, hh, nn, ss

  ! !  ! ! Check date2dec with John
  ! !  ! write(*,*) '#05: ', date2dec(22,12,2005,15,00,00)
  ! !  ! write(*,*) '#06: ', date2dec(22,12,2005,16,00,00)
  ! !  ! write(*,*) '#07: ', date2dec(22,12,2005,17,00,00)
  ! !  ! write(*,*) '#08: ', date2dec(22,12,2005,17,01,00)
  ! !  ! write(*,*) '#09: ', date2dec(22,12,2005,17,00,01)

END PROGRAM main
