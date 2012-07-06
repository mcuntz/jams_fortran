PROGRAM main
  
  USE mo_kind,        ONLY: dp, sp, i4
  USE mo_percentile,  ONLY: ksmallest, median, percentile
  
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
  isgood = isgood .and. (median(dat) == 5.5_dp)
  isgood = isgood .and. (median(dat,mask=(dat /= 10._dp)) == 5._dp)
  isgood = isgood .and. (ksmallest(dat,4) == 4._dp)
  isgood = isgood .and. (percentile(dat,95._dp) == 10._dp)
  dqua = percentile(dat,(/50._dp,95._dp/))
  isgood = isgood .and. (dqua(1) == 5._dp)
  isgood = isgood .and. (dqua(2) == 10._dp)
  dqua = percentile(dat,(/50._dp,75._dp/),mode_in=2_i4)
  isgood = isgood .and. (dqua(1) == 5._dp)
  isgood = isgood .and. (dqua(2) == 7.5_dp)

  if (isgood) then
     write(*,*) 'mo_percentile double precision o.k.'
  else
     write(*,*) 'mo_percentile double precision failed!'
  endif

  ! Single precision
  isgood = .true.
  sat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = isgood .and. (median(sat) == 5.5_sp)
  isgood = isgood .and. (median(sat,mask=(sat /= 10._sp)) == 5._sp)
  isgood = isgood .and. (ksmallest(sat,4) == 4._sp)
  isgood = isgood .and. (percentile(sat,95._sp) == 10._sp)
  squa = percentile(sat,(/50._sp,95._sp/))
  isgood = isgood .and. (squa(1) == 5._sp)
  isgood = isgood .and. (squa(2) == 10._sp)
  squa = percentile(sat,(/50._sp,75._sp/),mode_in=2_i4)
  isgood = isgood .and. (squa(1) == 5._sp)
  isgood = isgood .and. (squa(2) == 7.5_sp)

  if (isgood) then
     write(*,*) 'mo_percentile single precision o.k.'
  else
     write(*,*) 'mo_percentile single precision failed!'
  endif

END PROGRAM main
