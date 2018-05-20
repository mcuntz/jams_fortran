program test_kernel

  use mo_kind,   only: i4, dp
  use mo_kernel, only: kernel_regression_h, kernel_regression
  use mo_kernel, only: kernel_density_h, kernel_density
  use mo_kernel, only: kernel_cumdensity
  !$ USE omp_lib,    only: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
  ! use IEEE_exceptions

  implicit none

  integer(i4)                            :: ii, nn
  real(dp),    dimension(2)              :: h
  real(dp),    dimension(10,2)           :: x
  real(dp),    dimension(10)             :: y
  real(dp),    dimension(10)             :: ysmooth
  logical                                :: isgood
  integer(i4), dimension(10)             :: test1

  real(dp),    dimension(:), allocatable :: xin_1d      ! measured
  real(dp),    dimension(:), allocatable :: yin_1d      ! smoothed at x_in

  ! for OMP
  !$  integer(i4)                           :: n_threads

  !$  write(*,*) '--------------------------------------------------'
  !$  write(*,*) ' This program is parallel.'
  !$OMP parallel
  !$    n_threads = OMP_GET_NUM_THREADS()
  !$OMP end parallel
  !$  write(*,*) ' It will run using ',n_threads,' threads'
  !$  write(*,*) '--------------------------------------------------'

  nn = size(x,1)
  forall(ii=1:nn) x(ii,1) = real(ii-1,dp)/9.0_dp
  forall(ii=1:nn) x(ii,2) = 1.0_dp / (real(ii-1,dp)/9.0_dp + 0.1_dp)
  y(:) = 1.0_dp + x(:,1)**2 - sin(x(:,2))** 2

  isgood = .true.
  write(*,*) '  '

  ! ---------------------------------------------------------------------------------
  write(*,*) 'Testing approximation of bandwith h for regression'

  h = kernel_regression_h(x,y,silverman=.true.)
  write(*,*) '   Bandwith h (Silvermans rule)  : ', h
  if ( int(h(1)*10000.0_dp,i4) .ne.  2291_i4 ) isgood = .false.
  if ( int(h(2)*10000.0_dp,i4) .ne. 19033_i4 ) isgood = .false.

  h = kernel_regression_h(x,y,silverman=.false.)
  write(*,*) '   Bandwith h (Cross-validation) : ', h
  if ( int(h(1)*10000.0_dp,i4) .ne.  1726_i4 ) isgood = .false.
  if ( int(h(2)*10000.0_dp,i4) .ne. 95169_i4 ) isgood = .false.
  write(*,*) ' '

  ! ---------------------------------------------------------------------------------
  write(*,*) 'Testing kernel regression using bandwith h approximated by cross-validation'
  ysmooth = kernel_regression(x,y,silverman=.false.)

  write(*,*) '   x1        x2        y approx.   y input  '
  do ii=1,nn
     write(*,'(A4,F7.4,A3,F7.4,A3,F7.4,A5,F7.4)') '    ',x(ii,1), '   ',x(ii,2), '   ', ysmooth(ii), '     ', y(ii)
  end do
  write(*,*) ' '

  test1 = (/ 5224_i4, 5256_i4, 5417_i4, 5178_i4, 4764_i4, 4923_i4, 6034_i4, 7774_i4, 9545_i4, 10960_i4 /)
  do ii=1,10
     if (int(ysmooth(ii)*10000.0_dp,i4) .ne. test1(ii)) isgood = .false.
  end do

  ! ---------------------------------------------------------------------------------
  ! Test data: 20 temperature data of Guayaquil
  nn = 20_i4
  allocate(xin_1d(nn))
  allocate(yin_1d(nn))
  xin_1d = (/ 26.1_dp, 24.5_dp, 24.8_dp, 24.5_dp, 24.1_dp, &
       24.3_dp, 26.4_dp, 24.9_dp, 23.7_dp, 23.5_dp, &
       24.0_dp, 24.1_dp, 23.7_dp, 24.3_dp, 26.6_dp, &
       24.6_dp, 24.8_dp, 24.4_dp, 26.8_dp, 25.2_dp /)

  ! ---------------------------------------------------------------------------------
  write(*,*)'Testing approximation of bandwith h for density estimation'

  h(1)        = kernel_density_h(xin_1d(:),silverman=.false.)
  write(*,*) '   estimated h: ',h(1)

  if (int(h(1)*10000.0_dp,i4) .ne. 4582_i4) isgood = .false.

  write(*,*) ' '

  ! ---------------------------------------------------------------------------------
  write(*,*) 'Testing kernel prob. density (PDF) estimation'

  yin_1d = kernel_density(xin_1d, h=h(1), silverman=.false.)

  write(*,'(A4,A7,A4,A7)') '    ','x      ','    ','pdf(x) '
  do ii=1,size(yin_1d)
     write(*,'(A4,F7.4,A4,F8.5)') '    ', xin_1d(ii), '    ',yin_1d(ii)
  end do

  test1 = (/ 1260_i4, 4713_i4, 3900_i4, 4713_i4, 4443_i4, 4773_i4, 1498_i4, 3494_i4, 3035_i4, 2214_i4 /)
  do ii=1,10
     if (int(yin_1d(ii)*10000.0_dp,i4) .ne. test1(ii)) isgood = .false.
  end do

  write(*,*) ' '

  ! ---------------------------------------------------------------------------------
  write(*,*) 'Testing kernel cum. density (CDF) estimation'

  yin_1d = kernel_cumdensity(xin_1d, h=h(1), silverman=.false., romberg=.false.)

  write(*,'(A4,A7,A4,A7)') '    ','x      ','    ','cdf(x) '
  do ii=1,size(yin_1d)
     write(*,'(A4,F7.4,A4,F8.5)') '    ', xin_1d(ii), '    ', yin_1d(ii)
  end do

  ! With normalisation: test1 = (/ 8847_i4, 4544_i4, 6062_i4, 4544_i4, 2363_i4, 3437_i4, 9332_i4, 6490_i4, 0607_i4, 0_i4 /)
  test1 = (/ 8461_i4, 4744_i4, 6055_i4, 4744_i4, 2861_i4, 3789_i4, 8880_i4, 6425_i4, 1344_i4, 819_i4 /)
  do ii=1,10
     if (int(yin_1d(ii)*10000.0_dp,i4) .ne. test1(ii)) isgood = .false.
  end do

  write(*,*) ' '

  ! ---------------------------------------------------------------------------------
  if (isgood) then
     write(*,*) 'mo_kernel o.k.'
  else
     write(*,*) 'mo_kernel failed !'
  end if

end program test_kernel
