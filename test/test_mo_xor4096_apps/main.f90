!*******************************************************
!
!   TO TEST xor4096_apps
!
!*******************************************************
program xor4096_apps

  use mo_kind,         only: i4, i8, dp
  use mo_ansi_colors, only: color, c_red, c_green
  use mo_xor4096,      only: xor4096,xor4096g, n_save_state
  use mo_xor4096_apps, only: xor4096_array, xor4096_range, xor4096g_mvn

  implicit none

  integer(i8)                            :: seed1
  real(dp),    dimension(20)             :: rn1
  integer(i8), dimension(n_save_state)   :: save_state1

  integer(i8), dimension(3)              :: seed2
  real(dp),    dimension(20,3)           :: rn2
  integer(i8), dimension(3,n_save_state) :: save_state2

  integer(i8), dimension(2)              :: seed3
  real(dp),    dimension(20,2)           :: rn_mvn
  integer(i8), dimension(2,n_save_state) :: save_state3
  real(dp), dimension(2,2)               :: cov  ! covariance matrix
  real(dp), dimension(2,2)               :: lcho ! lower cholesky factor

  integer(i4) :: i
  logical     :: isgood
  integer     :: info

  external DPOTRF   ! Lapack routine for calculating the cholesky factor

  isgood = .true.

  ! *********************************************************************
  ! TEST XOR4096_ARRAY
  ! *********************************************************************

  write(*,*) ''
  write(*,*) 'Testing xor4096_array: one stream'

  ! initialize a stream and save its state
  seed1 = 123456_i8
  call xor4096(seed1, rn1(1), save_state=save_state1)
  seed1 = 0_i8
  ! now ready for usage

  call xor4096_array(rn1, save_state=save_state1)
  do i=1,20
     write(*,*) i, rn1(i)
  end do

  if ( nint(rn1(20)*10000.0_dp, i4) /= 8767_i4 ) isgood = .false.

  write(*,*) ''
  write(*,*) 'Testing xor4096_array: multiple streams'

  ! initialize multiple streams and save their states
  seed2 = (/ 123456_i8, 567362_i8, 4572417_i8 /)
  call xor4096(seed2, rn2(1,:), save_state=save_state2)
  seed2 = 0_i8
  ! now ready for usage

  call xor4096_array(rn2, save_state=save_state2)
  do i=1,20
     write(*,*) i, rn2(i,:)
  end do

  if ( nint(rn2(20,3)*10000.0_dp, i4) /= 8125_i4 ) isgood = .false.

  ! *********************************************************************
  ! TEST XOR4096_MVN
  ! *********************************************************************

  write(*,*) ''
  write(*,*) 'Testing xor4096_mvn: time series'

  ! initialize a stream and save its state
  seed3 = (/ 123456_i8, 154725_i8 /)

  call xor4096g(seed3, rn_mvn(1,:), save_state=save_state3)
  seed3 = 0_i8
  ! now ready for usage

  ! calculating lcho: Lower Cholesky factor of covariance matrix cov
  cov(1,:) = (/ 1.0_dp, 0.25_dp /)
  cov(2,:) = (/ 0.25_dp, 1.0_dp /)
  lcho = cov
  call DPOTRF('L', size(lcho,2), lcho, size(lcho,1), info) ! changes lower triangular (incl. main diagonal) of Lcho
  ! DONT forget to set upper triangular matrix to zero
  lcho(1,2) = 0._dp

  call xor4096g_mvn(lcho, rn_mvn, save_state=save_state3)
  do i=1,20
     write(*,*) i, rn_mvn(i,:)
  end do

  if ( nint(rn_mvn(20,2)*10000.0_dp, i4) /= 771_i4 ) isgood = .false.

  ! *********************************************************************
  ! TEST XOR4096_RANGE
  ! *********************************************************************

  write(*,*) ''
  write(*,*) 'Testing xor4096_range: one stream'

  ! initialize a stream and save its state
  seed1 = 123456_i8
  call xor4096(seed1, rn1(1), save_state=save_state1)
  seed1 = 0_i8
  ! now ready for usage

  call xor4096_range((/ 1.0_dp, 5.0_dp /), rn1, save_state=save_state1)
  do i=1,20
     write(*,*) i, rn1(i)
  end do

  if ( nint(rn1(20)*10000.0_dp, i4) /= 45068_i4 ) isgood = .false.

  write(*,*) ''
  write(*,*) 'Testing xor4096_range: multiple streams'

  ! initialize multiple streams and save their states
  seed2 = (/ 123456_i8, 567362_i8, 4572417_i8 /)
  call xor4096(seed2, rn2(1,:), save_state=save_state2)
  seed2 = 0_i8
  ! now ready for usage

  call xor4096_range((/ 1.0_dp, 5.0_dp /), rn2, save_state=save_state2)
  do i=1,20
     write(*,*) i, rn2(i,:)
  end do

  if ( nint(rn2(20,3)*10000.0_dp, i4) /= 42501_i4 ) isgood = .false.

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_xor4096_apps ', color('o.k.', c_green)
  else
     write(*,*) 'mo_xor4096_apps ', color('failed', c_red)
  end if

end program xor4096_apps
