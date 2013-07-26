program test
 
  USE mo_kind , ONLY: i4, i8, sp, dp
  USE mo_sobol, ONLY: sobol, sobol_array
  !
  IMPLICIT NONE
  !
  INTEGER(i4)                           :: ii
  INTEGER(i4)                           :: dims_i4
  INTEGER(i4)                           :: seed_i4
  REAL(sp), DIMENSION(:), ALLOCATABLE   :: sobol_set_sp
  INTEGER(i8)                           :: dims_i8
  INTEGER(i8)                           :: seed_i8
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: sobol_set_dp
  INTEGER(I4)                           :: nPoints_i4
  INTEGER(I4)                           :: skip_i4
  REAL(sp), DIMENSION(:,:), ALLOCATABLE :: sobol_database_sp
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: sobol_database_dp
  LOGICAL                               :: isgood

  write(*,*) ''
  write(*,*) 'Test mo_sobol.f90'
  isgood = .true.

  ! -----------------------------------------------
  ! Testing single precision sobol vector
  ! -----------------------------------------------
  dims_i4 = 3_i4
  allocate(sobol_set_sp(dims_i4))
  seed_i4 = 0_i4
  write(*,*) ''
  do ii=1,10
     write(*,*)  'call sobol( dims=', dims_i4,', seed=', seed_i4,', sobol_set_sp )'
     call sobol( dims_i4, seed_i4, sobol_set_sp )
     write(*,*) '  ', ii, sobol_set_sp
  end do
  if (any( nint(sobol_set_sp*10000,i4) .ne. (/ 6875_i4, 8125_i4, 8125_i4 /) )) isgood = .false.

  ! -----------------------------------------------
  ! Testing double precision sobol vector
  ! -----------------------------------------------
  dims_i8 = 3_i8
  allocate(sobol_set_dp(dims_i8))
  seed_i8 = 0_i8
  write(*,*) ''
  do ii=1,10
     write(*,*)  'call sobol( dims=', dims_i8,', seed=', seed_i8,', sobol_set_dp )'
     call sobol ( dims_i8, seed_i8, sobol_set_dp )
     write(*,*) '  ', ii, sobol_set_dp
  end do
  if (any( nint(sobol_set_dp*10000,i4) .ne. (/ 6875_i4, 8125_i4, 8125_i4 /) )) isgood = .false.

  ! -----------------------------------------------
  ! Testing single precision sobol database
  ! -----------------------------------------------
  dims_i4    = 3_i4
  nPoints_i4 = 4_i4
  skip_i4    = 5_i4
  allocate(sobol_database_sp(nPoints_i4, dims_i4))
  call sobol_array( dims_i4, nPoints_i4, skip_i4, sobol_database_sp )

  write(*,*) ''
  write(*,*)  'call sobol_array( dims=', dims_i4,', nPoints=', nPoints_i4,', skip=', skip_i4, ', sobol_database_sp )'
  do ii=1,nPoints_i4
     write(*,*) '  ', ii, sobol_database_sp(ii,:)
  end do
  if (any( nint(sobol_database_sp(nPoints_i4,:)*10000,i4) .ne. (/ 1875_i4, 3125_i4, 3125_i4 /) )) isgood = .false.

  ! -----------------------------------------------
  ! Testing double precision sobol database
  ! -----------------------------------------------
  dims_i4    = 3_i4
  nPoints_i4 = 4_i4
  skip_i4    = 5_i4
  allocate(sobol_database_dp(nPoints_i4, dims_i4))
  call sobol_array( dims_i4, nPoints_i4, skip_i4, sobol_database_dp )
  if (any( nint(sobol_database_dp(nPoints_i4,:)*10000,i4) .ne. (/ 1875_i4, 3125_i4, 3125_i4 /) )) isgood = .false.

  write(*,*) ''
  write(*,*) 'call sobol_array( dims=', dims_i4,', nPoints=', nPoints_i4,', skip=', skip_i4, ', sobol_database_dp )'
  do ii=1,nPoints_i4
     write(*,*) '  ', ii, sobol_database_dp(ii,:)
  end do

  ! deallocation
  deallocate(sobol_set_sp)
  deallocate(sobol_set_dp)
  deallocate(sobol_database_sp)
  deallocate(sobol_database_dp)

  if (isgood) then
     write(*,*) 'mo_sobol o.k.'
  else
     write(*,*) 'mo_sobol failed'
  endif
  !
end program test
