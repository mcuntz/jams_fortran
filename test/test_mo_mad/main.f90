PROGRAM main

  USE mo_kind, ONLY: dp, sp
  use mo_ansi_colors, only: color, c_red, c_green
  USE mo_mad,  ONLY: mad

  IMPLICIT NONE

  REAL(dp), DIMENSION(26) :: dat
  REAL(sp), DIMENSION(26) :: sat
  LOGICAL,  DIMENSION(26) :: mask, sollmask, mask2

  LOGICAL :: isgood

  Write(*,*) ''
  Write(*,*) 'Test mo_mad.f90'

  ! Double precision
  isgood = .true.
  dat = (/ -0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62, &
       2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30, &
       4.64,5.34,5.42,8.01 /)
  sollmask(:) = .true.
  mask        = mad(dat)
  if (any(mask .neqv. sollmask)) isgood = .false.
  sollmask(26) = .false.
  mask         = mad(dat,z=4._dp)
  if (any(mask .neqv. sollmask)) isgood = .false.
  sollmask(25) = .false.
  sollmask(1)  = .false.
  mask         = mad(dat,z=3._dp)
  if (any(mask .neqv. sollmask)) isgood = .false.
  sollmask(:)  = .true.
  sollmask(25) = .false.
  mask         = mad(dat,z=4._dp,deriv=2)
  if (any(mask .neqv. sollmask)) isgood = .false.
  sollmask(:)  = .true.
  sollmask(1)  = .false.
  sollmask(25) = .false.
  sollmask(26) = .false.
  mask(:)      = .true.
  mask(1)      = .false.
  mask2        = mad(dat,z=3._dp)
  mask         = mask .and. mask2
  if (any(mask .neqv. sollmask)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_mad double precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_mad double precision ', color('failed!', c_red)
  endif

  ! Single precision
  isgood = .true.
  sat = (/ -0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62, &
       2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30, &
       4.64,5.34,5.42,8.01 /)
  sollmask(:) = .true.
  mask        = mad(sat)
  if (.not. all(mask .eqv. sollmask)) isgood = .false.
  sollmask(26) = .false.
  mask         = mad(sat,z=4._sp)
  if (.not. all(mask .eqv. sollmask)) isgood = .false.
  sollmask(25) = .false.
  sollmask(1)  = .false.
  mask         = mad(sat,z=3._sp)
  if (.not. all(mask .eqv. sollmask)) isgood = .false.
  sollmask(:)  = .true.
  sollmask(25) = .false.
  mask         = mad(sat,z=4._sp,deriv=2)
  if (.not. all(mask .eqv. sollmask)) isgood = .false.
  sollmask(:)  = .true.
  sollmask(1)  = .false.
  sollmask(25) = .false.
  sollmask(26) = .false.
  mask(:)      = .true.
  mask(1)      = .false.
  mask2        = mad(sat, z=3._sp)
  mask         = mask .and. mask2
  if (.not. all(mask .eqv. sollmask)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_mad single precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_mad single precision ', color('failed!', c_red)
  endif

END PROGRAM main
