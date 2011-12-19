PROGRAM main
  
  USE mo_kind,     ONLY: dp, sp
  USE mo_interpol, ONLY: interpol
  
  IMPLICIT NONE
  
  REAL(dp), DIMENSION(10) :: dyin, dyout, dx, dsoll
  LOGICAL,  DIMENSION(10) :: mask
  REAL(sp), DIMENSION(10) :: syin, syout, sx, ssoll

  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_interpol.f90'

  ! Double precision
  isgood = .true.
  dsoll   = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  dx      = dsoll
  dyin    = dsoll
  mask(:) = .true.
  dyin(5) = -9.
  mask(5) = .false.
  dyout = interpol(pack(dyin,mask), pack(dx,mask), dx)
  if (any(abs(dsoll-dyout) > 0._dp)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_interpol double precision o.k.'
  else
     write(*,*) 'mo_interpol double precision failed.'
  endif

  ! Single precision
  isgood = .true.
  ssoll   = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  sx      = ssoll
  syin    = ssoll
  mask(:) = .true.
  syin(5) = -9.
  mask(5) = .false.
  syout = interpol(pack(syin,mask), pack(sx,mask), sx)
  if (any(abs(ssoll-syout) > 0._sp)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_interpol single precision o.k.'
  else
     write(*,*) 'mo_interpol single precision failed!'
  endif

END PROGRAM main
