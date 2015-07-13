PROGRAM main

  USE mo_kind,      ONLY: i4, dp, sp
  USE mo_utils,     ONLY: eq
  USE mo_orderpack, ONLY: refsor, mrgrnk, omedian, sort, sort_index

  IMPLICIT NONE

  REAL(dp),    DIMENSION(10) :: dat1, dat2
  REAL(sp),    DIMENSION(10) :: sat1, sat2
  INTEGER(i4), DIMENSION(10) :: ii

  LOGICAL :: isgood

  Write(*,*) ''
  Write(*,*) 'Test mo_orderpack.f90'

  ! Double precision
  ! sort
  dat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  dat2 = dat1(10:1:-1)
  call refsor(dat2)
  if (any(abs(dat1-dat2) > 0._dp)) isgood = .false.
  dat2 = dat1(10:1:-1)
  call mrgrnk(dat2, ii)
  if (any(abs(dat1-dat2(ii)) > 0._dp)) isgood = .false.
  ! sort aliases
  dat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  dat2 = dat1(10:1:-1)
  call sort(dat2)
  if (any(abs(dat1-dat2) > 0._dp)) isgood = .false.
  dat2 = dat1(10:1:-1)
  ii = sort_index(dat2)
  if (any(abs(dat1-dat2(ii)) > 0._dp)) isgood = .false.

  ! median
  isgood = isgood .and. eq(omedian(dat1),5.5_dp)

  if (isgood) then
     write(*,*) 'mo_orderpack double precision o.k.'
  else
     write(*,*) 'mo_orderpack double precision failed!'
  endif

  ! Single precision
  ! sort
  sat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  sat2 = sat1(10:1:-1)
  call refsor(sat2)
  if (any(abs(sat1-sat2) > 0._sp)) isgood = .false.
  sat2 = sat1(10:1:-1)
  call mrgrnk(sat2, ii)
  if (any(abs(sat1-sat2(ii)) > 0._sp)) isgood = .false.
  ! sort aliases
  sat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  sat2 = sat1(10:1:-1)
  call sort(sat2)
  if (any(abs(sat1-sat2) > 0._dp)) isgood = .false.
  sat2 = sat1(10:1:-1)
  ii = sort_index(sat2)
  if (any(abs(sat1-sat2(ii)) > 0._dp)) isgood = .false.

  ! median
  isgood = isgood .and. eq(omedian(sat1),5.5_sp)

  if (isgood) then
     write(*,*) 'mo_orderpack single precision o.k.'
  else
     write(*,*) 'mo_orderpack single precision failed!'
  endif

END PROGRAM main
