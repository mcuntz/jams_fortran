program test_mo_random_field

  use mo_kind,         only: i4, i8, dp
  use mo_random_field, only: random_velocity_field_gauss

  implicit none

  real(dp), dimension(:,:), allocatable  :: coord
  real(dp), dimension(2)                 :: corr_length
  real(dp), dimension(:,:), allocatable  :: velocity
  real(dp), dimension(:),   allocatable  :: potential
  real(dp), dimension(:,:), allocatable  :: cosine_mode

  integer(i4) :: i,n
  logical     :: isgood

  isgood = .true.

  n = 16
  allocate(coord(n,2))
  allocate(velocity(size(coord,1),size(coord,2)))
  allocate(potential(size(coord,1)))
  coord(:,1)  = (/ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
       0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, &
       0.2_dp, 0.2_dp, 0.2_dp, 0.2_dp, &
       0.3_dp, 0.3_dp, 0.3_dp, 0.3_dp /)
  coord(:,2)  = (/ 0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, &
       0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, &
       0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, &
       0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp /)
  velocity  = -9999.0_dp
  potential = -9999.0_dp

  corr_length = (/ 0.05_dp, 0.01_dp /)

  ! random field, i.e. cosine_modes, is determined and
  ! velocity (all dimensions) and potential (scalar) are calculated for 3 points
  velocity(1:3,:) = random_velocity_field_gauss(coord(1:3,:), corr_length, &
       sigma2=0.001_dp, seed=13546_i8, &
       cosine_modes_out=cosine_mode, &
       potential=potential(1:3) )

  ! previously determined modes are used to calculate velocity and potential 
  ! for additional points 4 to 10
  velocity(4:16,:) = random_velocity_field_gauss(coord(4:16,:), corr_length, &
       sigma2=0.001_dp, &
       cosine_modes_in=cosine_mode, &
       potential=potential(4:16) )
  
  write(*,*) ''
  write(*,*) ' x      y             velocity_x    velocity_y    potential'
  write(*,*) '--------------------------------------------------------------'
  do i=1, size(coord,1)
     write(*, '(F7.3,2X,F7.3,2X,F12.5,2X,F12.5,2X,F12.5)') coord(i,:), velocity(i,:), potential(i)
  end do
  write(*,*) ''

  ! write(*,*) 'velocity ', int(velocity(:,1)*100000.0_dp)
  if ( any(int(velocity(:,1)*100000.0_dp) .ne. (/-1286, 1541, -709,  173, 1188, -497, 1077, &
       2775, -2006, 2563, -369, -499, -1108,  992,  262,  235/)) ) isgood = .false.
  ! write(*,*) 'velocity ', int(velocity(:,2)*100000.0_dp)
  if ( any(int(velocity(:,2)*100000.0_dp) .ne. (/-13247, 2117, 2852, -2973, -6919, -3199, &
       -10483, 13322, -24272, 15308, -9172, -264, -2682, -6343, -13550, -9434/)) ) isgood = .false.
  ! write(*,*) 'potential ', int(potential(:)*100000.0_dp)
  if ( any(int(potential(:)*100000.0_dp) .ne. (/0, 71, -126, -215, -40, 58, 26, 127, 11, &
       -37, 70, 69, 103, 68, -57, -21/)) ) isgood = .false.
  write(*,*) ''

  if (isgood) then
     write(*,*) 'mo_random_field o.k.'
  else
     write(*,*) 'mo_random_field failed!'
  endif

end program test_mo_random_field
