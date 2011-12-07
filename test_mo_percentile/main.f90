PROGRAM main
  
  USE mo_kind,        ONLY: dp, sp
  USE mo_percentile,  ONLY: ksmallest, median, percentile, empqua
  
  IMPLICIT NONE
  
  REAL(dp), DIMENSION(10) :: dat
  REAL(dp), DIMENSION(2)  :: dqua
  REAL(sp), DIMENSION(10) :: sat
  REAL(sp), DIMENSION(2)  :: squa

  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_percentile.f90'

  ! Double precision
  isgood = .true.
  dat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  call empqua(dat,(/50._dp,95._dp/),dqua)
  isgood = isgood .and. (median(dat) == 5.5_dp)
  isgood = isgood .and. (median(dat,mask=(dat /= 10._dp)) == 5._dp)
  isgood = isgood .and. (ksmallest(dat,4) == 4._dp)
  isgood = isgood .and. (percentile(dat,95._dp) == 9._dp)
  isgood = isgood .and. (dqua(1) == 5._dp)
  isgood = isgood .and. (dqua(2) == 9._dp)

  if (isgood) then
     write(*,*) 'Double precision o.k.'
  else
     write(*,*) 'Double precision failed'
  endif

  ! Single precision
  isgood = .true.
  sat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  call empqua(sat,(/50._sp,95._sp/),squa)
  isgood = isgood .and. (median(sat) == 5.5_sp)
  isgood = isgood .and. (median(sat,mask=(sat /= 10._sp)) == 5._sp)
  isgood = isgood .and. (ksmallest(sat,4) == 4._sp)
  isgood = isgood .and. (percentile(sat,95._sp) == 9._sp)
  isgood = isgood .and. (squa(1) == 5._sp)
  isgood = isgood .and. (squa(2) == 9._sp)

  if (isgood) then
     write(*,*) 'Single precision o.k.'
  else
     write(*,*) 'Single precision failed'
  endif

END PROGRAM main
