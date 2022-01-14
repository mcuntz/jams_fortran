PROGRAM main

  USE mo_kind,        ONLY: dp, sp, i4, i8
  use mo_ansi_colors, only: color, c_red, c_green
  USE mo_percentile,  ONLY: median, percentile, qmedian, n_element
  USE mo_orderpack,   ONLY: sort
  USE mo_utils,       ONLY: eq, ne

  IMPLICIT NONE

  REAL(dp), DIMENSION(10) :: dat
  REAL(dp), DIMENSION(2)  :: dqua
  REAL(sp), DIMENSION(10) :: sat
  REAL(sp), DIMENSION(2)  :: squa

  real(sp) :: stmp
  real(dp) :: dtmp

#ifndef __ABSOFT__
#ifdef __GFORTRAN__
  ! If arrays re too big, the flag -frecursive gives either
  !   Illegal instruction: 4
  ! or
  !   Segmentation fault: 11
  ! The large numbers work with gfortran if not -frecursive. Then it states
  ! Warning: Array 'big' at (1) is larger than limit set by '-fmax-stack-var-size=',
  ! moved from stack to static storage. This makes the procedure unsafe when called
  ! recursively, or concurrently from multiple threads. Consider using '-frecursive',
  ! or increase the -fmax-stack-var-size=' limit, or change the code to use an
  ! ALLOCATABLE array. [-Wsurprising]
  integer,  parameter :: nele = 100000
#else
  integer,  parameter :: nele = 10000000
#endif
  real(dp), dimension(nele) :: big, buf
  real(dp) :: med
  integer  :: i
  integer(i8) :: istart, istop
#endif

  LOGICAL :: isgood

  write(*,*) ''
  write(*,*) 'Test mo_percentile.f90'

  ! Double precision
  isgood = .true.
  dat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  dtmp = median(dat)
  isgood = isgood .and. eq(dtmp, 5.5_dp)
  dtmp = median(dat, mask=ne(dat, 10._dp))
  isgood = isgood .and. eq(dtmp, 5._dp)
  dtmp = n_element(dat, 4)
  isgood = isgood .and. eq(dtmp, 4._dp)
  dtmp = percentile(dat, 95._dp)
  isgood = isgood .and. eq(dtmp, 10._dp)
  dqua = percentile(dat, (/50._dp, 95._dp/))
  isgood = isgood .and. eq(dqua(1), 5._dp)
  isgood = isgood .and. eq(dqua(2), 10._dp)
  dqua = percentile(dat,(/50._dp, 75._dp/), mode_in=1_i4)
  dqua = percentile(dat,(/50._dp, 75._dp/), mode_in=3_i4)
  dqua = percentile(dat,(/50._dp, 75._dp/), mode_in=4_i4)
  dqua = percentile(dat,(/50._dp, 75._dp/), mode_in=5_i4)
  dqua = percentile(dat,(/50._dp, 75._dp/), mode_in=6_i4)
  dqua = percentile(dat,(/50._dp, 75._dp/), mode_in=7_i4)
  dqua = percentile(dat,(/50._dp, 75._dp/), mode_in=8_i4)
  dqua = percentile(dat,(/50._dp, 75._dp/), mode_in=2_i4)
  isgood = isgood .and. eq(dqua(1), 5._dp)
  isgood = isgood .and. eq(dqua(2), 7.5_dp)

  if (isgood) then
     write(*,*) 'mo_percentile double precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_percentile double precision ', color('failed!', c_red)
  endif

  ! Single precision
  isgood = .true.
  sat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  stmp = median(sat)
  isgood = isgood .and. eq(stmp, 5.5_sp)
  stmp = median(sat, mask=ne(sat, 10._sp))
  isgood = isgood .and. eq(stmp, 5._sp)
  stmp = n_element(sat,4)
  isgood = isgood .and. eq(stmp, 4._sp)
  stmp = percentile(sat,95._sp)
  isgood = isgood .and. eq(stmp, 10._sp)
  squa = percentile(sat,(/50._sp, 95._sp/))
  isgood = isgood .and. eq(squa(1), 5._sp)
  isgood = isgood .and. eq(squa(2), 10._sp)
  squa = percentile(sat,(/50._sp, 75._sp/), mode_in=2_i4)
  isgood = isgood .and. eq(squa(1), 5._sp)
  isgood = isgood .and. eq(squa(2), 7.5_sp)

  if (isgood) then
     write(*,*) 'mo_percentile single precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_percentile single precision ', color('failed!', c_red)
  endif

! Absoft gives segmentation fault on all of the following test, do not know why
#ifndef __ABSOFT__
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
  med =  n_element(buf,nele/2+1)
  call system_clock(istop)
  write(*,*) "n_element: ", med, istop - istart

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
