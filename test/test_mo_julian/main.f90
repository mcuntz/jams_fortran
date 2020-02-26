program main

  use mo_kind,   only: i4, dp
  use mo_ansi_colors, only: color, c_red, c_green
  use mo_julian, only: ndays, ndyin
  use mo_julian, only: caldat, julday, date2dec, dec2date
  use mo_utils,  only: ne

  implicit none

  integer(i4) :: dd, mm, yy
  integer(i4) :: hh, nn, ss
  integer(i4) :: dd2, mm2, yy2
  integer(i4) :: hh2, nn2, ss2
  real(dp)    :: ff, rr, gg, d0
  integer(i4), dimension(3) :: dd1, mm1, yy1, hh1, nn1, ss1
  integer(i4), dimension(3) :: dd12, mm12, yy12, hh12, nn12, ss12
  real(dp),    dimension(3) :: dec1
  real(dp),    dimension(3) :: rr1, ff1, gg1, fd1, fd2
  integer(i4) :: jj

  logical :: isgood, allgood

  integer(i4),      parameter :: ncal = 10
  character(len=8), dimension(ncal) :: calendars
  integer(i4),      dimension(ncal) :: is1900
  integer(i4)       :: ical, dmin, dmax
  character(len=20) :: iscal
  real(dp)          :: eps

  write(*,*) ''
  write(*,*) 'Test mo_julian.f90'
  
  allgood = .true.

  calendars(1)  = 'standard'
  calendars(2)  = '360day'
  calendars(3)  = '365day'
  calendars(4)  = '366day'
  calendars(5)  = '360_day'
  calendars(6)  = '365_day'
  calendars(7)  = '366_day'
  calendars(8)  = 'noleap'
  calendars(9)  = 'all_leap'
  calendars(10) = 'lilian'

  is1900(1)  = 2415021
  is1900(2)  = 684000
  is1900(3)  = 693500
  is1900(4)  = 695400
  is1900(5)  = 2380320
  is1900(6)  = 2413380
  is1900(7)  = 2419992
  is1900(8)  = 2413380
  is1900(9)  = 2419992
  is1900(10) = 115861

  dmin = 2415021 - 366 ! < 1899
  dmax = 2524594       ! 01.01.2200 in julian days
  
  ! ----------------------
  ! IMSL tests

  isgood = .true.

  ! ndays, ndyin
  do jj=0, 109573 ! 01.01.2200
     call ndyin(jj,dd,mm,yy)
     ss = ndays(dd,mm,yy)
     if (jj /= ss) isgood = .false.
  end do
  if (ndays(01,01,1900) /= 0) isgood = .false.
  call ndyin(0, dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.
  if (ndays(1,1,-4712) /= -is1900(1)) isgood = .false.
  call ndyin(-is1900(1), dd, mm, yy)
  if ((dd /= 1) .or. (mm /= 1) .or. (yy /= -4712)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_julian ndyin/ndays ', color('o.k.', c_green)
  else
     write(*,*) 'mo_julian ndyin/ndays ', color('failed!', c_red)
  endif

  ! ----------------------
  ! calendar tests
  
  do ical=1, ncal

     iscal = calendars(ical)
     isgood = .true.

     ! julday, caldat
     do jj=dmin, dmax ! 01.01.2200
        call caldat(jj, dd, mm, yy, calendar=trim(iscal))
        ss = julday(dd, mm, yy, calendar=trim(iscal))
        if (jj /= ss) isgood = .false.
     end do
     if (julday(01, 01, 1900, calendar=trim(iscal)) /= is1900(ical)) isgood = .false.
     call caldat(is1900(ical), dd, mm, yy, calendar=trim(iscal))
     if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900)) isgood = .false.

     ! date2dec
     
     if (trim(iscal)=='lilian') then
        hh = 0
        ff = 0._dp
        gg = 0.5_dp
        eps = epsilon(1.0_dp)*(real(is1900(ical),dp)+2299159.5_dp)
     else
        hh = 12
        ff = 0.5_dp
        gg = 0._dp
        eps = epsilon(1.0_dp)*real(is1900(ical),dp)
     endif
     if (abs(date2dec(01,01,1900,hh,0, calendar=trim(iscal)) - real(is1900(ical),dp)) > eps*10.0_dp) isgood = .false.
     if (abs(date2dec(01,01,1900, calendar=trim(iscal)) - real(is1900(ical),dp) + ff) > eps*10.0_dp) isgood = .false.
     if (abs(date2dec(01,01,1900, fracday=0.5_dp, calendar=trim(iscal)) - real(is1900(ical),dp) - gg) > &
             eps*10.0_dp) isgood = .false.

     ! date2dec, dec2date - scalar
     do jj=dmin, dmax ! 01.01.2200
        call random_number(ff)
        rr = real(jj,dp) + ff
        call dec2date(rr,dd,mm,yy,hh,nn,ss, calendar=trim(iscal))
        gg = date2dec(dd,mm,yy,hh,nn,ss, calendar=trim(iscal))
        call dec2date(gg,dd2,mm2,yy2,hh2,nn2,ss2, calendar=trim(iscal))
        if ((dd /= dd2) .or. (mm /= mm2) .or. (yy /= yy2) .or. (hh /= hh2) .or. (nn /= nn2) .or. (ss /= ss2)) &
             isgood = .false.
     end do

     call dec2date(date2dec(01,01,1900,12,01,02, calendar=trim(iscal)), dd, mm, yy, hh, nn, ss, calendar=trim(iscal))
     if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 12) .or. (nn /= 01) .or. (ss /= 02)) isgood = .false.
     call dec2date(date2dec(01,01,1900, calendar=trim(iscal)), dd, mm, yy, hh, nn, ss, calendar=trim(iscal))
     if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 00) .or. (nn /= 00) .or. (ss /= 00)) isgood = .false.
     call dec2date(date2dec(01,01,1900,23,59,59, calendar=trim(iscal)), dd, mm, yy, hh, nn, ss, calendar=trim(iscal))
     if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 23) .or. (nn /= 59) .or. (ss /= 59)) isgood = .false.
     call dec2date(date2dec(01,01,01,12,12,12, calendar=trim(iscal)), dd, mm, yy, hh, nn, ss, calendar=trim(iscal))
     if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 01) .or. (hh /= 12) .or. (nn /= 12) .or. (ss /= 12)) isgood = .false.
     call dec2date(date2dec(01,02,-100,11,11,11, calendar=trim(iscal)), dd, mm, yy, hh, nn, ss, calendar=trim(iscal))
     if ((dd /= 1) .or. (mm /= 2) .or. (yy /= -100) .or. (hh /= 11) .or. (nn /= 11) .or. (ss /= 11)) isgood = .false.

     ! date2dec, dec2date - vector
     do jj=dmin, dmax ! 01.01.2200
        call random_number(ff1)
        rr1 = real(jj,dp) + ff1
        call dec2date(rr1,dd1,mm1,yy1,hh1,nn1,ss1,calendar=trim(iscal))
        gg1 = date2dec(dd1,mm1,yy1,hh1,nn1,ss1,calendar=trim(iscal))
        call dec2date(gg1,dd12,mm12,yy12,hh12,nn12,ss12,calendar=trim(iscal))
        if (any(dd1 /= dd12) .or. any(mm1 /= mm12) .or. any(yy1 /= yy12) .or. &
             any(hh1 /= hh12) .or. any(nn1 /= nn12) .or. any(ss1 /= ss12)) isgood = .false.
     end do

     call dec2date(date2dec((/01,10,27/),(/02,12,03/),(/-100,1,2000/),(/11,23,00/),(/11,57,47/),(/11,12,59/), &
          calendar=trim(iscal)), dd1, mm1, yy1, hh1, nn1, ss1, calendar=trim(iscal))
     if ((dd1(1) /= 1) .or. (mm1(1) /= 2) .or. (yy1(1) /= -100) .or. &
          (hh1(1) /= 11) .or. (nn1(1) /= 11) .or. (ss1(1) /= 11)) isgood = .false.
     if ((dd1(2) /= 10) .or. (mm1(2) /= 12) .or. (yy1(2) /= 1) .or. &
          (hh1(2) /= 23) .or. (nn1(2) /= 57) .or. (ss1(2) /= 12)) isgood = .false.
     if ((dd1(3) /= 27) .or. (mm1(3) /= 3) .or. (yy1(3) /= 2000) .or. &
          (hh1(3) /= 0) .or. (nn1(3) /= 47) .or. (ss1(3) /= 59)) isgood = .false.

     ! date2dec, dec2date - scalar/vector mix
     ! This does not work with all NAG compilers. Have to have temporary array.
     ! call dec2date360(date2dec360(01,12,2000,(/11,12,13/),11,59), dd1, mm1, yy1, hh1, nn1, ss1)
     do jj=dmin, dmax ! 01.01.2200
        call random_number(ff1)
        rr1 = real(jj,dp) + ff1
        call dec2date(rr1,dd1,mm1,yy1,hh1,nn1,ss1, calendar=trim(iscal))
        gg1 = date2dec(dd1(1),mm1(1),yy1(1),hh1,nn1,ss1, calendar=trim(iscal))
        call dec2date(gg1,dd12,mm12,yy12,hh12,nn12,ss12, calendar=trim(iscal))
        if (any(dd1(1) /= dd12) .or. any(mm1(1) /= mm12) .or. any(yy1(1) /= yy12) .or. &
             any(hh1 /= hh12) .or. any(nn1 /= nn12) .or. any(ss1 /= ss12)) isgood = .false.
     end do
     
     dec1 = date2dec(01,12,2000,(/11,12,13/),11,59, calendar=trim(iscal))
     call dec2date(dec1, dd1, mm1, yy1, hh1, nn1, ss1, calendar=trim(iscal))
     if ((dd1(1) /= 1) .or. (mm1(1) /= 12) .or. (yy1(1) /= 2000) .or. &
          (hh1(1) /= 11) .or. (nn1(1) /= 11) .or. (ss1(1) /= 59)) isgood = .false.
     if ((dd1(2) /= 1) .or. (mm1(2) /= 12) .or. (yy1(2) /= 2000) .or. &
          (hh1(2) /= 12) .or. (nn1(2) /= 11) .or. (ss1(2) /= 59)) isgood = .false.
     if ((dd1(3) /= 1) .or. (mm1(3) /= 12) .or. (yy1(3) /= 2000) .or. &
          (hh1(3) /= 13) .or. (nn1(3) /= 11) .or. (ss1(3) /= 59)) isgood = .false.

     ! Mix
     do jj=dmin, dmax ! 01.01.2200
        call caldat(jj,dd,mm,yy, calendar=trim(iscal))
        ss = int(date2dec(dd,mm,yy,12,0,0, calendar=trim(iscal)), i4)
        if (jj /= ss) isgood = .false.
     end do

     ! fracday
     do jj=dmin, dmax ! 01.01.2200
        call random_number(ff1)
        rr1 = real(jj,dp) + ff1
        call dec2date(rr1,dd1,mm1,yy1,fracday=fd1, calendar=trim(iscal))
        gg1 = date2dec(dd1(1),mm1(1),yy1(1),fracday=fd1, calendar=trim(iscal))
        call dec2date(gg1,dd12,mm12,yy12,fracday=fd2, calendar=trim(iscal))
        if (any(dd1(1) /= dd12) .or. any(mm1(1) /= mm12) .or. any(yy1(1) /= yy12) .or. any(ne(fd1,fd2))) &
             isgood = .false.
     end do

     ! units
     d0 = date2dec(1,1,1900,calendar=trim(iscal))
     do jj=is1900(ical), is1900(ical)+dmax-dmin
        call random_number(ff)
        rr = real(jj,dp) + ff
        call dec2date(rr,dd,mm,yy,hh,nn,ss,calendar=trim(iscal))
        rr = date2dec(dd,mm,yy,hh,nn,ss,calendar=trim(iscal)) ! include eps
        gg = rr - d0     ! days since 01.01.1900 00:00:00
        call dec2date(gg,dd2,mm2,yy2,hh2,nn2,ss2,units='days since 1900-01-01 00:00:00',calendar=trim(iscal))
        if ((dd /= dd2) .or. (mm /= mm2) .or. (yy /= yy2) .or. (hh /= hh2) .or. (nn /= nn2) .or. (ss /= ss2)) &
             isgood = .false.
        gg = gg * 24._dp ! hours since 01.01.1900 00:00:00
        call dec2date(gg,dd2,mm2,yy2,hh2,nn2,ss2,units='hours since 1900-1-01 00:00:00',calendar=trim(iscal))
        if ((dd /= dd2) .or. (mm /= mm2) .or. (yy /= yy2) .or. (hh /= hh2) .or. (nn /= nn2) .or. (ss /= ss2)) &
             isgood = .false.
        gg = gg * 60._dp ! minutes since 01.01.1900 00:00:00
        call dec2date(gg,dd2,mm2,yy2,hh2,nn2,ss2,units='minutes since 1900-01-1 00:00',calendar=trim(iscal))
        if ((dd /= dd2) .or. (mm /= mm2) .or. (yy /= yy2) .or. (hh /= hh2) .or. (nn /= nn2) .or. (ss /= ss2)) &
             isgood = .false.
        gg = gg * 60._dp ! seconds since 01.01.1900 00:00:00
        call dec2date(gg,dd2,mm2,yy2,hh2,nn2,ss2,units='seconds since 1900-1-1 00',calendar=trim(iscal))
        if ((dd /= dd2) .or. (mm /= mm2) .or. (yy /= yy2) .or. (hh /= hh2) .or. (nn /= nn2) .or. (ss /= ss2)) &
             isgood = .false.
     end do

     call dec2date(0.5_dp, dd, mm, yy, hh, nn, ss, units='days since 1900-1-1',calendar=trim(iscal))
     if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 12) .or. (nn /= 0) .or. (ss /= 0)) isgood = .false.
     call dec2date(12.0_dp, dd, mm, yy, hh, nn, ss, units='hours since 1900-01-01 12:00',calendar=trim(iscal))
     if ((dd /= 2) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 0) .or. (nn /= 0) .or. (ss /= 0)) isgood = .false.
     call dec2date(60.0_dp, dd, mm, yy, hh, nn, ss, units='minutes since 1900-01-01 12:00',calendar=trim(iscal))
     if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 13) .or. (nn /= 0) .or. (ss /= 0)) isgood = .false.
     call dec2date(59.0_dp, dd, mm, yy, hh, nn, ss, units='seconds since 1900-01-01 12:00',calendar=trim(iscal))
     if ((dd /= 1) .or. (mm /= 1) .or. (yy /= 1900) .or. (hh /= 12) .or. (nn /= 0) .or. (ss /= 59)) isgood = .false.

     allgood = allgood .and. isgood
     
     if (isgood) then
        write(*,*) 'mo_julian ', trim(iscal), ' ', color('o.k.', c_green)
     else
        write(*,*) 'mo_julian ', trim(iscal), ' ', color('failed!', c_red)
     endif
     
  end do

  if (allgood) then
     write(*,*) 'mo_julian ', color('o.k.', c_green)
  else
     write(*,*) 'mo_julian ', color('failed!', c_red)
  endif

  ! ! Check 59 minutes, 60 seconds problem if eps too small in dec2date
  ! dec1(1) = 2443879.5_dp
  ! write(*,*) '#00: ', dec1(1), 1._dp/24._dp, date2dec(1,1,1990,1,0,0) - date2dec(1,1,1990,0,0,0)
  ! call dec2date(dec1(1), dd, mm, yy, hh, nn, ss)
  ! write(*,*) '#01: ', dec1(1), dd, mm, yy, hh, nn, ss
  ! dec1(1) = dec1(1) + 1._dp/24._dp
  ! call dec2date(dec1(1), dd, mm, yy, hh, nn, ss)
  ! write(*,*) '#02: ', dec1(1), dd, mm, yy, hh, nn, ss
  ! dec1(1) = dec1(1) - 1._dp/24._dp
  ! call dec2date(dec1(1), dd, mm, yy, hh, nn, ss)
  ! write(*,*) '#03: ', dec1(1), dd, mm, yy, hh, nn, ss
  ! dec1(1) = dec1(1) + date2dec(1,1,1990,1,0,0) - date2dec(1,1,1990,0,0,0)
  ! call dec2date(dec1(1), dd, mm, yy, hh, nn, ss)
  ! write(*,*) '#04: ', dec1(1), dd, mm, yy, hh, nn, ss

  ! ! Check date2dec with John
  ! write(*,*) '#05: ', date2dec(22,12,2005,15,00,00)
  ! write(*,*) '#06: ', date2dec(22,12,2005,16,00,00)
  ! write(*,*) '#07: ', date2dec(22,12,2005,17,00,00)
  ! write(*,*) '#08: ', date2dec(22,12,2005,17,01,00)
  ! write(*,*) '#09: ', date2dec(22,12,2005,17,00,01)

  ! ! check python
  ! dd2 = 16
  ! mm2 = 1
  ! yy2 = 2004
  ! hh2 = 12
  ! nn2 = 0
  ! ss2 = 0
  ! gg = date2dec(16,01,2004,12,0,calendar='360day') - date2dec(01,01,1850,calendar='360day')
  ! print*, 'days of 16.01.2004 12:00 from 01.01.1850 00:00: ', gg

END PROGRAM main
