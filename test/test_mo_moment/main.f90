PROGRAM main
  
  USE mo_kind,   ONLY: dp, sp
  USE mo_moment, ONLY: absdev, average, central_moment, central_moment_var, correlation, &
       covariance, kurtosis, mean, mixed_central_moment, mixed_central_moment_var, &
       moment, skewness, stddev, variance
  USE mo_utils, only: ne
  
  IMPLICIT NONE
  
  REAL(dp), DIMENSION(10) :: dat1, dat2
  REAL(dp) :: a_dp, v_dp, s_dp, k_dp, m_dp, std_dp, ad_dp

  REAL(sp), DIMENSION(10) :: sat1, sat2
  REAL(sp) :: a_sp, v_sp, s_sp, k_sp, m_sp, std_sp, ad_sp

  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_moment.f90'

  ! Double precision
  dat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  ! Test standard
  isgood = isgood .and. (nint(10._dp*mean(dat1)) == 55)
  isgood = isgood .and. (nint(10._dp*average(dat1)) == 55)
  isgood = isgood .and. (nint(10._dp*variance(dat1)) == 92)
  isgood = isgood .and. (nint(10._dp*stddev(dat1)) == 30)
  isgood = isgood .and. (nint(10._dp*skewness(dat1)) == 00)
  isgood = isgood .and. (nint(10._dp*kurtosis(dat1)) == -16)
  isgood = isgood .and. (nint(10._dp*absdev(dat1)) == 25)
  isgood = isgood .and. (nint(100._dp*central_moment(dat1,2)) == 825)
  isgood = isgood .and. (nint(10._dp*central_moment_var(dat1,2)) == 53)
  ! Test the mask
  isgood = isgood .and. (nint(10._dp*mean(dat1(1:9))) == 50)
  isgood = isgood .and. (nint(10._dp*average(dat1(1:9))) == 50)
  isgood = isgood .and. (nint(10._dp*variance(dat1(1:9))) == 75)
  isgood = isgood .and. (nint(10._dp*stddev(dat1(1:9))) == 27)
  isgood = isgood .and. (nint(10._dp*skewness(dat1(1:9))) == 00)
  isgood = isgood .and. (nint(10._dp*kurtosis(dat1(1:9))) == -16)
  isgood = isgood .and. (nint(10._dp*absdev(dat1(1:9))) == 22)
  isgood = isgood .and. (nint(10._dp*central_moment(dat1(1:9),2)) == 67)
  isgood = isgood .and. (nint(10._dp*central_moment_var(dat1(1:9),2)) == 38)
  isgood = isgood .and. (nint(10._dp*mean(dat1, mask=ne(dat1,10._dp))) == 50)
  isgood = isgood .and. (nint(10._dp*average(dat1, mask=ne(dat1,10._dp))) == 50)
  isgood = isgood .and. (nint(10._dp*variance(dat1, mask=ne(dat1,10._dp))) == 75)
  isgood = isgood .and. (nint(10._dp*stddev(dat1, mask=ne(dat1,10._dp))) == 27)
  isgood = isgood .and. (nint(10._dp*skewness(dat1, mask=ne(dat1,10._dp))) == 00)
  isgood = isgood .and. (nint(10._dp*kurtosis(dat1, mask=ne(dat1,10._dp))) == -16)
  isgood = isgood .and. (nint(10._dp*absdev(dat1, mask=ne(dat1,10._dp))) == 22)
  isgood = isgood .and. (nint(10._dp*central_moment(dat1, 2, mask=ne(dat1,10._dp))) == 67)
  isgood = isgood .and. (nint(10._dp*central_moment_var(dat1, 2, mask=ne(dat1,10._dp))) == 38)
  ! Test moment
  call moment(dat1, average=a_dp, variance=v_dp, skewness=s_dp, kurtosis=k_dp, mean=m_dp, &
       stddev=std_dp, absdev=ad_dp, mask=ne(dat1,10._dp))
  isgood = isgood .and. (nint(10._dp*m_dp) == 50)
  isgood = isgood .and. (nint(10._dp*a_dp) == 50)
  isgood = isgood .and. (nint(10._dp*v_dp) == 75)
  isgood = isgood .and. (nint(10._dp*std_dp) == 27)
  isgood = isgood .and. (nint(10._dp*s_dp) == 00)
  isgood = isgood .and. (nint(10._dp*k_dp) == -16)
  isgood = isgood .and. (nint(10._dp*ad_dp) == 22)
  ! Test the double input functions
  dat2 = (/ 1.3, 2.7, 3.9, 5.5, 7., 8.8, 11., 13., 15., 16. /)
  isgood = isgood .and. (nint(10._dp*correlation(dat1,dat2)) == 9)
  isgood = isgood .and. (nint(10._dp*covariance(dat1,dat2)) == 141)
  isgood = isgood .and. (nint(10._dp*mixed_central_moment(dat1,dat2,1,1)) == 141)
  isgood = isgood .and. (nint(10._dp*mixed_central_moment_var(dat1,dat2,1,1)) == 142)
  isgood = isgood .and. (nint(10._dp*correlation(dat1(1:9),dat2(1:9))) == 9)
  isgood = isgood .and. (nint(10._dp*covariance(dat1(1:9),dat2(1:9))) == 115)
  isgood = isgood .and. (nint(10._dp*mixed_central_moment(dat1(1:9),dat2(1:9),1,1)) == 115)
  isgood = isgood .and. (nint(10._dp*mixed_central_moment_var(dat1(1:9),dat2(1:9),1,1)) == 113)
  isgood = isgood .and. (nint(10._dp*correlation(dat1,dat2, mask=(ne(dat1,10._dp).and.ne(dat2,10._dp)))) == 9)
  isgood = isgood .and. (nint(10._dp*covariance(dat1,dat2, mask=(ne(dat1,10._dp).and.ne(dat2,10._dp)))) == 115)
  isgood = isgood .and. (nint(10._dp*mixed_central_moment(dat1,dat2,1,1, mask=(ne(dat1,10._dp).and.ne(dat2,10._dp)))) == 115)
  isgood = isgood .and. (nint(10._dp*mixed_central_moment_var(dat1,dat2,1,1, mask=(ne(dat1,10._dp).and.ne(dat2,10._dp)))) == 113)

  if (isgood) then
     write(*,*) 'mo_moment double precision o.k.'
  else
     write(*,*) 'mo_moment double precision failed!'
  endif

  ! Single precision
  sat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  isgood = .true.
  ! Test standard
  isgood = isgood .and. (nint(10._sp*mean(sat1)) == 55)
  isgood = isgood .and. (nint(10._sp*average(sat1)) == 55)
  isgood = isgood .and. (nint(10._sp*variance(sat1)) == 92)
  isgood = isgood .and. (nint(10._sp*stddev(sat1)) == 30)
  isgood = isgood .and. (nint(10._sp*skewness(sat1)) == 00)
  isgood = isgood .and. (nint(10._sp*kurtosis(sat1)) == -16)
  isgood = isgood .and. (nint(10._sp*absdev(sat1)) == 25)
  isgood = isgood .and. (nint(100._sp*central_moment(sat1,2)) == 825)
  isgood = isgood .and. (nint(10._sp*central_moment_var(sat1,2)) == 53)
  ! Test the mask
  isgood = isgood .and. (nint(10._sp*mean(sat1(1:9))) == 50)
  isgood = isgood .and. (nint(10._sp*average(sat1(1:9))) == 50)
  isgood = isgood .and. (nint(10._sp*variance(sat1(1:9))) == 75)
  isgood = isgood .and. (nint(10._sp*stddev(sat1(1:9))) == 27)
  isgood = isgood .and. (nint(10._sp*skewness(sat1(1:9))) == 00)
  isgood = isgood .and. (nint(10._sp*kurtosis(sat1(1:9))) == -16)
  isgood = isgood .and. (nint(10._sp*absdev(sat1(1:9))) == 22)
  isgood = isgood .and. (nint(10._sp*central_moment(sat1(1:9),2)) == 67)
  isgood = isgood .and. (nint(10._sp*central_moment_var(sat1(1:9),2)) == 38)
  isgood = isgood .and. (nint(10._sp*mean(sat1, mask=ne(sat1,10._sp))) == 50)
  isgood = isgood .and. (nint(10._sp*average(sat1, mask=ne(sat1,10._sp))) == 50)
  isgood = isgood .and. (nint(10._sp*variance(sat1, mask=ne(sat1,10._sp))) == 75)
  isgood = isgood .and. (nint(10._sp*stddev(sat1, mask=ne(sat1,10._sp))) == 27)
  isgood = isgood .and. (nint(10._sp*skewness(sat1, mask=ne(sat1,10._sp))) == 00)
  isgood = isgood .and. (nint(10._sp*kurtosis(sat1, mask=ne(sat1,10._sp))) == -16)
  isgood = isgood .and. (nint(10._sp*absdev(sat1, mask=ne(sat1,10._sp))) == 22)
  isgood = isgood .and. (nint(10._sp*central_moment(sat1, 2, mask=ne(sat1,10._sp))) == 67)
  isgood = isgood .and. (nint(10._sp*central_moment_var(sat1, 2, mask=ne(sat1,10._sp))) == 38)
  ! Test moment
  call moment(sat1, average=a_sp, variance=v_sp, skewness=s_sp, kurtosis=k_sp, mean=m_sp, &
       stddev=std_sp, absdev=ad_sp, mask=ne(sat1,10._sp))
  isgood = isgood .and. (nint(10._sp*m_sp) == 50)
  isgood = isgood .and. (nint(10._sp*a_sp) == 50)
  isgood = isgood .and. (nint(10._sp*v_sp) == 75)
  isgood = isgood .and. (nint(10._sp*std_sp) == 27)
  isgood = isgood .and. (nint(10._sp*s_sp) == 00)
  isgood = isgood .and. (nint(10._sp*k_sp) == -16)
  isgood = isgood .and. (nint(10._sp*ad_sp) == 22)
  ! Test the double input functions
  sat2 = (/ 1.3, 2.7, 3.9, 5.5, 7., 8.8, 11., 13., 15., 16. /)
  isgood = isgood .and. (nint(10._sp*correlation(sat1,sat2)) == 9)
  isgood = isgood .and. (nint(10._sp*covariance(sat1,sat2)) == 141)
  isgood = isgood .and. (nint(10._sp*mixed_central_moment(sat1,sat2,1,1)) == 141)
  isgood = isgood .and. (nint(10._sp*mixed_central_moment_var(sat1,sat2,1,1)) == 142)
  isgood = isgood .and. (nint(10._sp*correlation(sat1(1:9),sat2(1:9))) == 9)
  isgood = isgood .and. (nint(10._sp*covariance(sat1(1:9),sat2(1:9))) == 115)
  isgood = isgood .and. (nint(10._sp*mixed_central_moment(sat1(1:9),sat2(1:9),1,1)) == 115)
  isgood = isgood .and. (nint(10._sp*mixed_central_moment_var(sat1(1:9),sat2(1:9),1,1)) == 113)
  isgood = isgood .and. (nint(10._sp*correlation(sat1,sat2, mask=(ne(sat1,10._sp).and.ne(sat2,10._sp)))) == 9)
  isgood = isgood .and. (nint(10._sp*covariance(sat1,sat2, mask=(ne(sat1,10._sp).and.ne(sat2,10._sp)))) == 115)
  isgood = isgood .and. (nint(10._sp*mixed_central_moment(sat1,sat2,1,1, mask=(ne(sat1,10._sp).and.ne(sat2,10._sp)))) == 115)
  isgood = isgood .and. (nint(10._sp*mixed_central_moment_var(sat1,sat2,1,1, mask=(ne(sat1,10._sp).and.ne(sat2,10._sp)))) == 113)

  if (isgood) then
     write(*,*) 'mo_moment single precision o.k.'
  else
     write(*,*) 'mo_moment single precision failed!'
  endif

END PROGRAM main
