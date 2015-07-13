PROGRAM main
  
  USE mo_kind,           ONLY: i4, i8, dp, sp
  USE mo_utils,          ONLY: ne
  USE mo_linear_algebra, ONLY: diag, inverse, solve_linear_equations, solve_linear_equations_svd
  
  IMPLICIT NONE
  
  ! data
  INTEGER(i4), PARAMETER :: nn = 100
  REAL(dp),    DIMENSION(nn,nn) :: dat
  REAL(sp),    DIMENSION(nn,nn) :: sat
  INTEGER(i4), DIMENSION(nn,nn) :: i4at
  INTEGER(i8), DIMENSION(nn,nn) :: i8at
  ! diag
  REAL(dp),    DIMENSION(:), allocatable :: ddat
  REAL(sp),    DIMENSION(:), allocatable :: dsat
  INTEGER(i4), DIMENSION(:), allocatable :: di4at
  INTEGER(i8), DIMENSION(:), allocatable :: di8at
  ! inverse
  REAL(dp),    DIMENSION(:,:), allocatable :: idat
  REAL(sp),    DIMENSION(:,:), allocatable :: isat
  ! lin eq
  REAL(dp),    DIMENSION(nn) :: bdat
  REAL(sp),    DIMENSION(nn) :: bsat
  REAL(dp),    DIMENSION(:), allocatable :: xdat
  REAL(sp),    DIMENSION(:), allocatable :: xsat
  
  INTEGER(i4) :: i, nseed
  INTEGER(i4), DIMENSION(:), allocatable :: iseed
  LOGICAL :: isgood, allgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_linear_algebra.f90'

  ! reproducible results because routines can easily fail if badly conditioned
  call random_seed(size=nseed)
  allocate(iseed(nseed))
  forall(i=1:nseed) iseed(i) = i*10
  call random_seed(put=iseed)
  deallocate(iseed)

  call random_number(dat)
  call random_number(sat)
  i4at = int(sat, i4)
  i8at = int(dat, i8)
  allgood = .true.

  ! diag
  isgood = .true.
#ifdef pgiFortran
  allocate(ddat(nn))
  allocate(dsat(nn))
  allocate(di4at(nn))
  allocate(di8at(nn))
#endif
  ddat  = diag(dat)
  dsat  = diag(sat)
  di4at = diag(i4at)
  di8at = diag(i8at)
  if (size(ddat) /= nn)  isgood = .false.
  if (size(dsat) /= nn)  isgood = .false.
  if (size(di4at) /= nn) isgood = .false.
  if (size(di8at) /= nn) isgood = .false.
  do i=1, nn
     if (ne(ddat(i),dat(i,i)))  isgood = .false.
     if (ne(dsat(i),sat(i,i)))  isgood = .false.
     if (di4at(i) /= i4at(i,i)) isgood = .false.
     if (di8at(i) /= i8at(i,i)) isgood = .false.
  end do

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra diag failed!'

  ! inverse
  isgood = .true.
#ifdef pgiFortran
  allocate(idat(nn,nn))
  allocate(isat(nn,nn))
#endif
  idat = inverse(dat)
  isat = inverse(sat)
  ! allow eps in each element of diag -> very close but can be slightly higher -> *100
  ! print*, 'Inv ', abs(sum(diag(matmul(dat,idat))) - real(size(dat,1),dp)), (real(size(dat,1),dp)*epsilon(1.0_dp))*100._dp
  ! print*, 'Inv ', abs(sum(diag(matmul(sat,isat))) - real(size(sat,1),sp)), (real(size(sat,1),sp)*epsilon(1.0_sp))*100._dp
  ! print*, 'Inv ', abs(sum(matmul(dat,idat)) - real(size(dat,1),dp)), (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp
  ! print*, 'Inv ', abs(sum(matmul(sat,isat)) - real(size(sat,1),sp)), (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._dp
  if (abs(sum(diag(matmul(dat,idat))) - real(size(dat,1),dp)) > (real(size(dat,1),dp)*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (abs(sum(diag(matmul(sat,isat))) - real(size(sat,1),sp)) > (real(size(sat,1),sp)*epsilon(1.0_sp))*100._sp) isgood = .false.
  ! allow eps in each element of matrix -> very close but can be slightly higher -> *10
  if (abs(sum(matmul(dat,idat)) - real(size(dat,1),dp)) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (abs(sum(matmul(sat,isat)) - real(size(sat,1),sp)) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra inverse failed!'

  ! linear eq with LU decomposition
  isgood = .true.
  bdat = matmul(dat,ddat)
  bsat = matmul(sat,dsat)
#ifdef pgiFortran
  allocate(xdat(nn))
  allocate(xsat(nn))
#endif
  xdat = solve_linear_equations(dat, bdat)
  xsat = solve_linear_equations(sat, bsat)
  ! allow eps in each element of matrix -> very close but can be slightly higher -> *100
  ! print*, 'LU  ', sum(abs(xdat-ddat)), (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp
  ! print*, 'LU  ', sum(abs(xsat-dsat)), (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp
  ! print*, 'LU  ', sum(matmul(dat,xdat)-bdat), (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp
  ! print*, 'LU  ', sum(matmul(sat,xsat)-bsat), (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp
  if (sum(abs(xdat-ddat)) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(abs(xsat-dsat)) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.
  if (sum(matmul(dat,xdat)-bdat) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(matmul(sat,xsat)-bsat) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra solve_linear_equations failed!'

  ! linear eq with SVD
  isgood = .true.
  bdat = matmul(dat,ddat)
  bsat = matmul(sat,dsat)
  xdat = solve_linear_equations_svd(dat, bdat)!, condition=.false.)
  xsat = solve_linear_equations_svd(sat, bsat)!, condition=.false.)
  ! allow eps in each element of matrix -> very close but can be slightly higher -> *100
  ! print*, 'SVD ', sum(abs(xdat-ddat)), (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp
  ! print*, 'SVD ', sum(abs(xsat-dsat)), (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp
  ! print*, 'SVD ', sum(matmul(dat,xdat)-bdat), (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp
  ! print*, 'SVD ', sum(matmul(sat,xsat)-bsat), (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp
  if (sum(abs(xdat-ddat)) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(abs(xsat-dsat)) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.
  if (sum(matmul(dat,xdat)-bdat) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(matmul(sat,xsat)-bsat) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra solve_linear_equations_svd failed!'


  allgood = allgood .and. isgood
  if (allgood) then
     write(*,*) 'mo_linear_algebra o.k.'
  else
     write(*,*) 'mo_linear_algebra failed!'
  endif

END PROGRAM main
