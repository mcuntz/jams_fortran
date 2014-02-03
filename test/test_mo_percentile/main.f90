PROGRAM main
  
  USE mo_kind,        ONLY: dp, sp, i4
  USE mo_percentile,  ONLY: ksmallest, median, percentile, qmedian
  USE mo_sort,        ONLY: sort
  
  IMPLICIT NONE
  
  REAL(dp), DIMENSION(10) :: dat
  REAL(dp), DIMENSION(2)  :: dqua
  REAL(sp), DIMENSION(10) :: sat
  REAL(sp), DIMENSION(2)  :: squa

#ifndef ABSOFT
  integer,  parameter :: nele = 10000000
  real(dp), dimension(nele) :: big, buf
  real(dp) :: med
  integer  :: i, istart, istop
#endif

  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_percentile.f90'

  ! Double precision
  isgood = .true.
  dat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = isgood .and. (abs(median(dat) - 5.5_dp) .lt. tiny(1.0_dp))
  isgood = isgood .and. (abs(median(dat,mask=(abs(dat - 10._dp) .gt. tiny(1.0_dp))) - 5._dp) .lt. tiny(1.0_dp))
  isgood = isgood .and. (abs(ksmallest(dat,4) - 4._dp) .lt. tiny(1.0_dp))
  isgood = isgood .and. (abs(percentile(dat,95._dp) - 10._dp) .lt. tiny(1.0_dp))
  dqua = percentile(dat,(/50._dp,95._dp/))
  isgood = isgood .and. (abs(dqua(1) - 5._dp) .lt. tiny(1.0_dp))
  isgood = isgood .and. (abs(dqua(2) - 10._dp) .lt. tiny(1.0_dp))
  dqua = percentile(dat,(/50._dp,75._dp/),mode_in=1_i4)
  dqua = percentile(dat,(/50._dp,75._dp/),mode_in=3_i4)
  dqua = percentile(dat,(/50._dp,75._dp/),mode_in=4_i4)
  dqua = percentile(dat,(/50._dp,75._dp/),mode_in=5_i4)
  dqua = percentile(dat,(/50._dp,75._dp/),mode_in=6_i4)
  dqua = percentile(dat,(/50._dp,75._dp/),mode_in=7_i4)
  dqua = percentile(dat,(/50._dp,75._dp/),mode_in=8_i4)
  dqua = percentile(dat,(/50._dp,75._dp/),mode_in=2_i4)
  isgood = isgood .and. (abs(dqua(1) - 5._dp) .lt. tiny(1.0_dp))
  isgood = isgood .and. (abs(dqua(2) - 7.5_dp) .lt. tiny(1.0_dp))

  if (isgood) then
     write(*,*) 'mo_percentile double precision o.k.'
  else
     write(*,*) 'mo_percentile double precision failed!'
  endif

  ! Single precision
  isgood = .true.
  sat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = isgood .and. (abs(median(sat) - 5.5_sp) .lt. tiny(1.0_sp))
  isgood = isgood .and. (abs(median(sat,mask=(abs(sat - 10._dp) .gt. tiny(1.0_dp))) - 5._sp) .lt. tiny(1.0_sp))
  isgood = isgood .and. (abs(ksmallest(sat,4) - 4._sp) .lt. tiny(1.0_sp))
  isgood = isgood .and. (abs(percentile(sat,95._sp) - 10._sp) .lt. tiny(1.0_sp))
  squa = percentile(sat,(/50._sp,95._sp/))
  isgood = isgood .and. (abs(squa(1) - 5._sp) .lt. tiny(1.0_sp))
  isgood = isgood .and. (abs(squa(2) - 10._sp) .lt. tiny(1.0_sp))
  squa = percentile(sat,(/50._sp,75._sp/),mode_in=2_i4)
  isgood = isgood .and. (abs(squa(1) - 5._sp) .lt. tiny(1.0_sp))
  isgood = isgood .and. (abs(squa(2) - 7.5_sp) .lt. tiny(1.0_sp))

  if (isgood) then
     write(*,*) 'mo_percentile single precision o.k.'
  else
     write(*,*) 'mo_percentile single precision failed!'
  endif

! Absoft gives segmentation fault on all of the following test, do not know why
#ifndef ABSOFT
  ! Test speed
  do i = 1, nele
     call random_number(big(i))
  enddo
  ! forall(i=1:nele) big(i) = real(i,dp)
  ! forall(i=1:nele/2) big(2*i) = real(i,dp)

  buf = big
  call system_clock(istart)
  med =  median(buf)
  call system_clock(istop)
  write(*,*) "median: ", med, istop - istart

  buf = big
  call system_clock(istart)
  med =  ksmallest(buf,nele/2+1)
  call system_clock(istop)
  write(*,*) "ksmallest: ", med, istop - istart

  buf = big
  call system_clock(istart)
  med =  median(buf)
  call system_clock(istop)
  write(*,*) "median: ", med, istop - istart

  buf = big
  call system_clock(istart)
  med = qmedian(buf)
  call system_clock(istop)
  write(*,*) "qmedian: ", med, istop - istart

  buf = big
  call system_clock(istart)
  call sort(buf)
  med = buf(nele/2+1)
  call system_clock(istop)
  write(*,*) "sort: ", med, istop - istart
#endif

END PROGRAM main
