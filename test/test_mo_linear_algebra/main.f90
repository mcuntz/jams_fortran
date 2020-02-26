PROGRAM main

  USE mo_kind,           ONLY: i4, i8, dp, sp
  USE mo_utils,          ONLY: ne
  use mo_ansi_colors, only: color, c_red, c_green
  USE mo_linear_algebra, ONLY: diag, min_diag, banded, inverse, solve_linear_equations, solve_linear_equations_svd, &
                               solve_linear_equations_band

  IMPLICIT NONE

  ! data
  ! Borders so that single precision versions do not work anymore are
  !   70/80 for inverse
  !   120/130 for solve_linear_equations
  !   680/690 for solve_linear_equations_svd
  ! Now, single precision versions are computed via double precision.
  ! This improves inverse and solve_linear_equations but not solve_linear_equations_svd.
  ! Now, solve_linear_equations (LU factorisation) seems more stable then SVD version.
  INTEGER(i4), PARAMETER :: nn = 100
  INTEGER(i4), PARAMETER :: l  = 6
  INTEGER(i4), PARAMETER :: u  = 4
  INTEGER(i4), PARAMETER :: mm = l+u+1 ! for the banded-algorithms
  REAL(dp),    DIMENSION(nn,nn) :: dat
  REAL(sp),    DIMENSION(nn,nn) :: sat
  REAL(dp),    DIMENSION(nn,nn) :: dat2
  REAL(sp),    DIMENSION(nn,nn) :: sat2
  INTEGER(i4), DIMENSION(nn,nn) :: i4at
  INTEGER(i8), DIMENSION(nn,nn) :: i8at
  ! diag
  REAL(dp),    DIMENSION(:), allocatable :: ddat
  REAL(sp),    DIMENSION(:), allocatable :: dsat
  INTEGER(i4), DIMENSION(:), allocatable :: di4at
  INTEGER(i8), DIMENSION(:), allocatable :: di8at
  ! inverse/banded
  REAL(dp),    DIMENSION(:,:), allocatable :: idat, b_dat!, tdat ! test invers
  REAL(sp),    DIMENSION(:,:), allocatable :: isat, b_sat!, tsat
  ! lin eq
  REAL(dp),    DIMENSION(nn) :: bdat!, tbdat ! test linear system
  REAL(sp),    DIMENSION(nn) :: bsat!, tbsat
  REAL(dp),    DIMENSION(:), allocatable :: xdat
  REAL(sp),    DIMENSION(:), allocatable :: xsat

  INTEGER(i4) :: i, ii, nseed
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
  allocate(ddat(nn))
  allocate(dsat(nn))
  allocate(di4at(nn))
  allocate(di8at(nn))
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
  if (.not. isgood) write(*,*) 'mo_linear_algebra diag ', color('failed!', c_red)

  !min_diag
  isgood = .true.
  do ii=-nn+1,nn-1
    deallocate(ddat)
    deallocate(dsat)
    deallocate(di4at)
    deallocate(di8at)
    allocate(ddat(nn-abs(ii)))
    allocate(dsat(nn-abs(ii)))
    allocate(di4at(nn-abs(ii)))
    allocate(di8at(nn-abs(ii)))
    ddat  = min_diag(dat,ii)
    dsat  = min_diag(sat,ii)
    di4at = min_diag(i4at,ii)
    di8at = min_diag(i8at,ii)
    do i=1, nn-abs(ii)
      if (ne(ddat(i),  dat(i-min(ii,0),i+max(ii,0)))) isgood = .false.
      if (ne(dsat(i),  sat(i-min(ii,0),i+max(ii,0)))) isgood = .false.
      if (di4at(i) /= i4at(i-min(ii,0),i+max(ii,0)))  isgood = .false.
      if (di8at(i) /= i8at(i-min(ii,0),i+max(ii,0)))  isgood = .false.
    end do
  end do

  deallocate(ddat)
  deallocate(dsat)
  allocate(ddat(nn)) ! use later
  allocate(dsat(nn))
  ddat = diag(dat)
  dsat = diag(sat)

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra min_diag ', color('failed!', c_red)

  !banded
  isgood = .true.
  allocate(b_dat(mm,nn))
  allocate(b_sat(mm,nn))
  b_dat = banded(dat,l,u)
  b_sat = banded(sat,l,u)
  do i=1,u+1
    if (any(ne(b_dat(i,u+2-i:nn),  min_diag(dat,u+1-i)))) isgood = .false.
    if (any(ne(b_sat(i,u+2-i:nn),  min_diag(sat,u+1-i)))) isgood = .false.
  end do
  do i=1,l
    if (any(ne(b_dat(u+1+i,1:nn-i),  min_diag(dat,-i)))) isgood = .false.
    if (any(ne(b_sat(u+1+i,1:nn-i),  min_diag(sat,-i)))) isgood = .false.
 end do

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra banded ', color('failed!', c_red)

  ! inverse
  isgood = .true.
  allocate(idat(nn,nn))
  allocate(isat(nn,nn))
  idat = inverse(dat)
  isat = inverse(sat)
  ! allow eps in each element of diag -> very close but can be slightly higher -> *100
  if (abs(sum(diag(matmul(dat,idat))) - real(size(dat,1),dp)) > (real(size(dat,1),dp)*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (abs(sum(diag(matmul(sat,isat))) - real(size(sat,1),sp)) > (real(size(sat,1),sp)*epsilon(1.0_sp))*100._sp) isgood = .false.
  ! allow eps in each element of matrix -> very close but can be slightly higher -> *10
  if (abs(sum(matmul(dat,idat)) - real(size(dat,1),dp)) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (abs(sum(matmul(sat,isat)) - real(size(sat,1),sp)) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.
  ! tdat = matmul(dat,idat)
  ! do i=2, nn-1
  !    print*, tdat(i-1,i), tdat(i,i), tdat(i+1,i)
  ! enddo
  ! tsat = matmul(sat,isat)
  ! do i=2, nn-1
  !    print*, tsat(i-1,i), tsat(i,i), tsat(i+1,i)
  ! enddo

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra inverse ', color('failed!', c_red)

  ! linear eq with LU decomposition
  isgood = .true.
  bdat = matmul(dat,ddat)
  bsat = matmul(sat,dsat)
  allocate(xdat(nn))
  allocate(xsat(nn))
  xdat = solve_linear_equations(dat, bdat)
  xsat = solve_linear_equations(sat, bsat)
  ! allow eps in each element of matrix -> very close but can be slightly higher -> *100
  if (sum(abs(xdat-ddat)) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(abs(xsat-dsat)) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.
  if (sum(matmul(dat,xdat)-bdat) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(matmul(sat,xsat)-bsat) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.
  ! tbdat = matmul(dat,xdat)-bdat
  ! print*, tbdat
  ! tbsat = matmul(sat,xsat)-bsat
  ! print*, tbsat

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra solve_linear_equations ', color('failed!', c_red)

  ! linear eq with SVD
  isgood = .true.
  bdat = matmul(dat,ddat)
  bsat = matmul(sat,dsat)
  xdat = solve_linear_equations_svd(dat, bdat)!, condition=.false.)
  xsat = solve_linear_equations_svd(sat, bsat)!, condition=.false.)
  ! allow eps in each element of matrix -> very close but can be slightly higher -> *100
  if (sum(abs(xdat-ddat)) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(abs(xsat-dsat)) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.
  if (sum(matmul(dat,xdat)-bdat) > (real(size(dat,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(matmul(sat,xsat)-bsat) > (real(size(sat,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.
  ! tbdat = matmul(dat,xdat)-bdat
  ! print*, tbdat
  ! tbsat = matmul(sat,xsat)-bsat
  ! print*, tbsat

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra solve_linear_equations_svd ', color('failed!', c_red)

  ! linear eq with banded
  isgood = .true.

  dat2 = 0._dp
  sat2 = 0._dp
  do ii=-l,u
    do i=1,nn-abs(ii)
      dat2(i-min(ii,0),i+max(ii,0)) = dat(i-min(ii,0),i+max(ii,0))
      sat2(i-min(ii,0),i+max(ii,0)) = sat(i-min(ii,0),i+max(ii,0))
    end do
  end do
  bdat = matmul(dat2,ddat)
  bsat = matmul(sat2,dsat)
  xdat = solve_linear_equations_band(banded(dat2,l,u), bdat, l, u)!, condition=.false.)
  xsat = solve_linear_equations_band(banded(sat2,l,u), bsat, l, u)!, condition=.false.)
  ! allow eps in each element of matrix -> very close but can be slightly higher -> *100
  if (sum(abs(xdat-ddat)) > (real(size(dat2,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(abs(xsat-dsat)) > (real(size(sat2,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.
  if (sum(matmul(dat2,xdat)-bdat) > (real(size(dat2,1),dp)**2*epsilon(1.0_dp))*100._dp) isgood = .false.
  if (sum(matmul(sat2,xsat)-bsat) > (real(size(sat2,1),sp)**2*epsilon(1.0_sp))*100._sp) isgood = .false.

  allgood = allgood .and. isgood
  if (.not. isgood) write(*,*) 'mo_linear_algebra solve_linear_equations_band ', color('failed!', c_red)


  if (allgood) then
     write(*,*) 'mo_linear_algebra ', color('o.k.', c_green)
  else
     write(*,*) 'mo_linear_algebra ', color('failed!', c_red)
  endif

END PROGRAM main
