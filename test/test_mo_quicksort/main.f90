PROGRAM main
  
  USE mo_kind,      ONLY: i4, sp, dp
  use mo_timer,     only: timers_init, timer_start, timer_stop, timer_print, timer_clear
  USE mo_sort,      ONLY: sort, sort_index
  USE mo_orderpack, ONLY: mrgrnk, refsor
  use mo_orderpack, only: osort => sort, osort_index => sort_index
  use mo_quicksort, only: qsort, qsort_index, qsortc, quick_sort, qsortb, qsortmp
  use mo_quicksort, only: isort => sort, isort_index => sort_index, sortmp, sortmp_index
  use mo_utils,     only: ne

  IMPLICIT NONE
  
  integer(i4), parameter :: nn = 1000
  integer(i4), parameter :: nloop = 100

  REAL(dp),    DIMENSION(nn) :: dat, dat2
  INTEGER(i4), DIMENSION(nn) :: ii, jj
  integer(i4) :: i
  REAL(dp),    DIMENSION(nn) :: sat, sat2
  integer(i4), DIMENSION(nn) :: iat, iat2

  LOGICAL :: isgood

  
  Write(*,*) ''
  Write(*,*) 'Test mo_quicksort.f90'

  call timers_init()

  isgood = .true.
  Write(*,*) ''
  Write(*,*) 'Numerical Recipes'
  call timer_start(1)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     call sort(sat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     ii = sort_index(sat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  sat = sat(ii)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  sat = real(new_data(nn),sp)
  sat2 = sat
  ii = sort_index(sat)
  call sort(sat)
  if (any(ne(sat2(ii),sat))) isgood = .false.

  Write(*,*) 'Orderpack'
  call timer_start(1)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     call osort(sat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     ii = osort_index(sat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  sat = sat(ii)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  sat = real(new_data(nn),sp)
  sat2 = sat
  ii = osort_index(sat)
  call osort(sat)
  if (any(ne(sat2(ii),sat))) isgood = .false.

  Write(*,*) 'Wikibooks'
  call timer_start(1)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     call qsort(sat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     ii = qsort_index(sat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  sat = sat(ii)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  sat = real(new_data(nn),sp)
  ii = qsort_index(sat)
  call qsort(sat, jj)
  if (any(ii /= jj)) isgood = .false.

  sat = real(new_data(nn),sp)
  sat2 = sat
  ii = qsort_index(sat)
  call qsort(sat)
  if (any(ne(sat2(ii),sat))) isgood = .false.

  call timer_start(1)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     call isort(sat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     ii = isort_index(sat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  sat = sat(ii)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  sat = real(new_data(nn),sp)
  sat2 = sat
  ii = isort_index(sat)
  call isort(sat)
  if (any(ne(sat2(ii),sat))) isgood = .false.

  Write(*,*) 'David Bal - OpenMP'
  call timer_start(1)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     call qsortmp(sat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     call qsortmp(sat,ii)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  sat = real(new_data(nn),sp)
  ii = qsort_index(sat)
  call qsortmp(sat, jj)
  if (any(ii /= jj)) isgood = .false.

  sat = real(new_data(nn),sp)
  sat2 = sat
  call qsortmp(sat, ii)
  if (any(ne(sat2(ii),sat))) isgood = .false.

  call timer_start(1)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     call sortmp(sat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     sat = real(new_data(nn),sp)
     ii = sortmp_index(sat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  sat = sat(ii)
  do i=1, nn-1
     if (sat(i) > sat(i+1)) isgood= .false.
  end do

  sat = real(new_data(nn),sp)
  sat2 = sat
  ii = sortmp_index(sat)
  call sortmp(sat)
  if (any(ne(sat2(ii),sat))) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_quicksort single precision o.k.'
  else
     write(*,*) 'mo_quicksort single precision failed!'
  endif


  isgood = .true.
  Write(*,*) ''
  Write(*,*) 'Numerical Recipes'
  call timer_start(1)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     call sort(iat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     ii = sort_index(iat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  iat = iat(ii)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  iat = int(new_data(nn)*real(nn,dp),i4)
  iat2 = iat
  ii = sort_index(iat)
  call sort(iat)
  if (any(iat2(ii) /= iat)) isgood = .false.

  Write(*,*) 'Orderpack'
  call timer_start(1)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     call osort(iat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     ii = osort_index(iat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  iat = iat(ii)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  iat = int(new_data(nn)*real(nn,dp),i4)
  iat2 = iat
  ii = osort_index(iat)
  call osort(iat)
  if (any(iat2(ii) /= iat)) isgood = .false.

  Write(*,*) 'Wikibooks'
  call timer_start(1)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     call qsort(iat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     ii = qsort_index(iat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  iat = iat(ii)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  iat = int(new_data(nn)*real(nn,dp),i4)
  ii = qsort_index(iat)
  call qsort(iat, jj)
  if (any(ii /= jj)) isgood = .false.

  iat = int(new_data(nn)*real(nn,dp),i4)
  iat2 = iat
  ii = qsort_index(iat)
  call qsort(iat)
  if (any(iat2(ii) /= iat)) isgood = .false.

  call timer_start(1)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     call isort(iat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     ii = isort_index(iat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  iat = iat(ii)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  iat = int(new_data(nn)*real(nn,dp),i4)
  iat2 = iat
  ii = isort_index(iat)
  call isort(iat)
  if (any(iat2(ii) /= iat)) isgood = .false.

  Write(*,*) 'David Bal - OpenMP'
  call timer_start(1)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     call qsortmp(iat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     call qsortmp(iat,ii)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  iat = int(new_data(nn)*real(nn,dp),i4)
  iat2 = iat
  call qsortmp(iat, ii)
  if (any(iat2(ii) /= iat)) isgood = .false.

  call timer_start(1)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     call sortmp(iat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     iat = int(new_data(nn)*real(nn,dp),i4)
     ii = sortmp_index(iat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  iat = iat(ii)
  do i=1, nn-1
     if (iat(i) > iat(i+1)) isgood= .false.
  end do

  iat = int(new_data(nn)*real(nn,dp),i4)
  iat2 = iat
  ii = sortmp_index(iat)
  call sortmp(iat)
  if (any(iat2(ii) /= iat)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_quicksort integer o.k.'
  else
     write(*,*) 'mo_quicksort integer failed!'
  endif


  isgood = .true.
  Write(*,*) ''
  Write(*,*) 'Numerical Recipes'
  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call sort(dat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     dat = new_data(nn)
     ii = sort_index(dat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  dat = dat(ii)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  dat = new_data(nn)
  dat2 = dat
  ii = sort_index(dat)
  call sort(dat)
  if (any(ne(dat2(ii),dat))) isgood = .false.

  Write(*,*) 'Orderpack'
  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call refsor(dat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     dat = new_data(nn)
     call mrgrnk(dat, ii)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  dat = dat(ii)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  dat = new_data(nn)
  dat2 = dat
  call mrgrnk(dat, ii)
  call refsor(dat)
  if (any(ne(dat2(ii),dat))) isgood = .false.

  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call osort(dat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     dat = new_data(nn)
     ii = osort_index(dat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  dat = dat(ii)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  dat = new_data(nn)
  dat2 = dat
  ii = osort_index(dat)
  call osort(dat)
  if (any(ne(dat2(ii),dat))) isgood = .false.

  Write(*,*) 'Wikibooks'
  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call qsort(dat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     dat = new_data(nn)
     call qsort(dat,ii)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  call timer_start(3)
  do i=1, nloop
     dat = new_data(nn)
     ii = qsort_index(dat)
  end do
  call timer_stop(3)
  call timer_print(3)
  call timer_clear(3)
  dat = dat(ii)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  dat = new_data(nn)
  ii = qsort_index(dat)
  call qsort(dat, jj)
  if (any(ii /= jj)) isgood = .false.

  dat = new_data(nn)
  dat2 = dat
  ii = qsort_index(dat)
  call qsort(dat)
  if (any(ne(dat2(ii),dat))) isgood = .false.

  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call isort(dat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     dat = new_data(nn)
     ii = isort_index(dat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  dat = dat(ii)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  dat = new_data(nn)
  dat2 = dat
  ii = isort_index(dat)
  call isort(dat)
  if (any(ne(dat2(ii),dat))) isgood = .false.

  ! Fails with sorted arrays
  Write(*,*) 'Cormen et al. (1997)'
  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call qsortc(dat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  Write(*,*) 'Alan Miller'
  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call quick_sort(dat, ii)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do


  Write(*,*) 'David Bal'
  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call qsortb(dat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  Write(*,*) 'David Bal - OpenMP'
  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call qsortmp(dat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     dat = new_data(nn)
     call qsortmp(dat,ii)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  dat = real(new_data(nn),sp)
  ii = qsort_index(dat)
  call qsortmp(dat, jj)
  if (any(ii /= jj)) isgood = .false.

  dat = real(new_data(nn),sp)
  dat2 = dat
  call qsortmp(dat, ii)
  if (any(ne(dat2(ii),dat))) isgood = .false.

  call timer_start(1)
  do i=1, nloop
     dat = new_data(nn)
     call sortmp(dat)
  end do
  call timer_stop(1)
  call timer_print(1)
  call timer_clear(1)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  call timer_start(2)
  do i=1, nloop
     dat = new_data(nn)
     ii = sortmp_index(dat)
  end do
  call timer_stop(2)
  call timer_print(2)
  call timer_clear(2)
  dat = dat(ii)
  do i=1, nn-1
     if (dat(i) > dat(i+1)) isgood= .false.
  end do

  dat = new_data(nn)
  dat2 = dat
  ii = sortmp_index(dat)
  call sortmp(dat)
  if (any(ne(dat2(ii),dat))) isgood = .false.


  if (isgood) then
     write(*,*) 'mo_quicksort double precision o.k.'
  else
     write(*,*) 'mo_quicksort double precision failed!'
  endif

contains

  function new_data(n)

    implicit none

    integer(i4) :: n
    real(dp), dimension(n) :: new_data

    ! integer(i4) :: i
    ! forall(i=1:n) new_data(i) = real(n-i,dp)
    call random_number(new_data)

  end function new_data

END PROGRAM main
