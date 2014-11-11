program test_mo_select_distant_entries

  use mo_kind,                   only: i4, dp
  use mo_select_distant_entries, only: select_distant_entries

  implicit none

  real(dp), dimension(:,:),  allocatable :: distance_matrix        ! square matrix with distance between element i and j 
  integer(i4)                            :: n_total                ! total number of elements
  integer(i4)                            :: n_select               ! number of elements to select
  integer(i4), dimension(:), allocatable :: max_entries
  logical                                :: isgood
  
  isgood = .true. 
  n_total  = 5
  n_select = 3

  allocate( distance_matrix(n_total, n_total) )
  allocate( max_entries(n_select) )

  !             1                   2
  !             
  !             
  !                       5
  !             
  !     
  !             4                   3

  distance_matrix(1,:) = (/ 0.0_dp, 5.0_dp, 7.0_dp, 5.0_dp, 3.5_dp /)
  distance_matrix(2,:) = (/ 5.0_dp, 0.0_dp, 5.0_dp, 7.0_dp, 3.5_dp /)
  distance_matrix(3,:) = (/ 7.0_dp, 5.0_dp, 0.0_dp, 5.0_dp, 3.5_dp /)
  distance_matrix(4,:) = (/ 5.0_dp, 7.0_dp, 5.0_dp, 0.0_dp, 3.5_dp /)
  distance_matrix(5,:) = (/ 3.5_dp, 3.5_dp, 3.5_dp, 3.5_dp, 0.0_dp /)

  ! selects first pair as the one with overall largest distance
  max_entries = select_distant_entries( distance_matrix, n_select, first=0 )
  write(*,*) 'Elements selected (first pair as the one with overall largest distance): ', max_entries
  if ( .not. all(max_entries .eq. (/1,2,3/)) ) isgood = .false.

  ! selects first pair randomly
  max_entries = select_distant_entries( distance_matrix, n_select, first=1, seed=12345_i4 )
  write(*,*) 'Elements selected (first pair randomly):                                 ', max_entries
  if ( .not. all(max_entries .eq. (/2,4,5/)) ) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_select_distant_entries: o.k.'
  else
     write(*,*) 'mo_select_distant_entries: failed!'
  endif

end program test_mo_select_distant_entries
