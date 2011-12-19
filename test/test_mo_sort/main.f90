PROGRAM main
  
  USE mo_kind, ONLY: i4, dp, sp
  USE mo_sort, ONLY: sort, sort_index
  
  IMPLICIT NONE
  
  REAL(dp),    DIMENSION(10) :: dat1, dat2
  REAL(sp),    DIMENSION(10) :: sat1, sat2
  INTEGER(i4), DIMENSION(10) :: ii

  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_sort.f90'

  ! Double precision
  dat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  dat2 = dat1(10:1:-1)
  call sort(dat2)
  if (any(abs(dat1-dat2) > 0._dp)) isgood = .false.
  dat2 = dat1(10:1:-1)
  ii = sort_index(dat2)
  if (any(abs(dat1-dat2(ii)) > 0._dp)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_sort double precision o.k.'
  else
     write(*,*) 'mo_sort double precision failed!'
  endif

  ! Single precision
  sat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  sat2 = sat1(10:1:-1)
  call sort(sat2)
  if (any(abs(sat1-sat2) > 0._sp)) isgood = .false.
  sat2 = sat1(10:1:-1)
  ii = sort_index(sat2)
  if (any(abs(sat1-sat2(ii)) > 0._sp)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_sort single precision o.k.'
  else
     write(*,*) 'mo_sort single precision failed!'
  endif

END PROGRAM main
