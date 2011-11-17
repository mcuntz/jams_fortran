PROGRAM main
  
  USE mo_kind,   ONLY: dp, sp
  USE mo_moment, ONLY: absdev, average, central_moment, central_moment_var, correlation, &
       covariance, kurtosis, mean, mixed_central_moment, mixed_central_moment_var, &
       moment, skewness, stddev, variance
  
  IMPLICIT NONE
  
  REAL(dp), DIMENSION(10) :: dat1, dat2
  REAL(dp) :: a_dp, v_dp, s_dp, k_dp, m_dp, std_dp, ad_dp
  REAL(sp), DIMENSION(10) :: sat1, sat2
  REAL(sp) :: a_sp, v_sp, s_sp, k_sp, m_sp, std_sp, ad_sp
  
  Write(*,*) ''
  Write(*,*) 'Test mo_moment.f90'

  ! Double precision
  Write(*,*) ''
  Write(*,*) '  Test double precision'

  ! Test the single input functions
  dat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  Write(*,'(A,10f5.1)') '     dat1 is: ', dat1
  Write(*,'(A,f5.1)') '     Mean       of dat1 (target=5.5):  ', mean(dat1)
  Write(*,'(A,f5.1)') '     Average    of dat1 (target=5.5):  ', average(dat1)
  Write(*,'(A,f5.1)') '     Variance   of dat1 (target=9.2):  ', variance(dat1)
  Write(*,'(A,f5.1)') '     Stddev     of dat1 (target=3.0):  ', stddev(dat1)
  Write(*,'(A,f5.1)') '     Skewness   of dat1 (target=0.0):  ', skewness(dat1)
  Write(*,'(A,f5.1)') '     Kurtosis   of dat1 (target=-1.6): ', kurtosis(dat1)
  Write(*,'(A,f5.1)') '     MeanAbsDev of dat1 (target=2.5):  ', absdev(dat1)
  Write(*,'(A,f5.1)') '     2nd Central Moment             of dat1 (target=8.2): ', central_moment(dat1,2)
  Write(*,'(A,f5.1)') '     Variance of 2nd central moment of dat1 (target=5.3): ', central_moment_var(dat1,2)
  ! Test the mask
  Write(*,*) '    Without last'
  Write(*,'(A,f5.1)') '     Mean       of dat1 (target=5.0):  ', mean(dat1(1:9))
  Write(*,'(A,f5.1)') '     Average    of dat1 (target=5.0):  ', average(dat1(1:9))
  Write(*,'(A,f5.1)') '     Variance   of dat1 (target=7.5):  ', variance(dat1(1:9))
  Write(*,'(A,f5.1)') '     Stddev     of dat1 (target=2.7):  ', stddev(dat1(1:9))
  Write(*,'(A,f5.1)') '     Skewness   of dat1 (target=0.0):  ', skewness(dat1(1:9))
  Write(*,'(A,f5.1)') '     Kurtosis   of dat1 (target=-1.6): ', kurtosis(dat1(1:9))
  Write(*,'(A,f5.1)') '     MeanAbsDev of dat1 (target=2.2):  ', absdev(dat1(1:9))
  Write(*,'(A,f5.1)') '     2nd Central Moment             of dat1 (target=6.7): ', central_moment(dat1(1:9),2)
  Write(*,'(A,f5.1)') '     Variance of 2nd central moment of dat1 (target=3.8): ', central_moment_var(dat1(1:9),2)
  Write(*,*) '    With 10 masked'
  Write(*,'(A,f5.1)') '     Mean       of dat1 (target=5.0):  ', mean(dat1, mask=(dat1 /= 10.))
  Write(*,'(A,f5.1)') '     Average    of dat1 (target=5.5):  ', average(dat1, mask=(dat1 /= 10.))
  Write(*,'(A,f5.1)') '     Variance   of dat1 (target=7.5):  ', variance(dat1, mask=(dat1 /= 10.))
  Write(*,'(A,f5.1)') '     Stddev     of dat1 (target=2.7):  ', stddev(dat1, mask=(dat1 /= 10.))
  Write(*,'(A,f5.1)') '     Skewness   of dat1 (target=0.0):  ', skewness(dat1, mask=(dat1 /= 10.))
  Write(*,'(A,f5.1)') '     Kurtosis   of dat1 (target=-1.6): ', kurtosis(dat1, mask=(dat1 /= 10.))
  Write(*,'(A,f5.1)') '     MeanAbsDev of dat1 (target=2.2):  ', absdev(dat1, mask=(dat1 /= 10.))
  Write(*,'(A,f5.1)') '     2nd Central Moment             of dat1 (target=6.7): ', &
       central_moment(dat1,2, mask=(dat1 /= 10.))
  Write(*,'(A,f5.1)') '     Variance of 2nd central moment of dat1 (target=3.8): ', &
       central_moment_var(dat1,2, mask=(dat1 /= 10.))
  call moment(dat1, average=a_dp, variance=v_dp, skewness=s_dp, kurtosis=k_dp, mean=m_dp, &
       stddev=std_dp, absdev=ad_dp, mask=(dat1 /= 10.))
  Write(*,'(A,f5.1,1X,f5.1,1X,f5.1,1X,f5.1,1X,f5.1,1X,f5.1,1X,f5.1,1X)') '     Moments of a are: ', &
       m_dp, a_dp, v_dp, std_dp, s_dp, k_dp, ad_dp
  ! 
  ! Test the double input functions
  Write(*,*) ''
  dat2 = (/ 1.3, 2.7, 3.9, 5.5, 7., 8.8, 11., 13., 15., 16. /)
  Write(*,'(A,10f5.1)') '     dat2 is: ', dat2
  Write(*,'(A,f5.2)') '     Correlation                          of dat1 and dat2 (target=0.90):  ', &
       correlation(dat1,dat2)
  Write(*,'(A,f5.2)') '     Covariance                           of dat1 and dat2 (target=14.11): ', &
       covariance(dat1,dat2)
  Write(*,'(A,f5.2)') '     Mixed central moment 1,1             of dat1 and dat2 (target=14.11): ', &
       mixed_central_moment(dat1,dat2,1,1)
  Write(*,'(A,f5.2)') '     Variance of mixed central moment 1,1 of dat1 and dat2 (target=14.24): ', &
       mixed_central_moment_var(dat1,dat2,1,1)
  Write(*,*) '    Without last'
  Write(*,'(A,f5.2)') '     Correlation                          of dat1 and dat2 (target=0.88):  ', &
       correlation(dat1(1:9),dat2(1:9))
  Write(*,'(A,f5.2)') '     Covariance                           of dat1 and dat2 (target=11.47): ', &
       covariance(dat1(1:9),dat2(1:9))
  Write(*,'(A,f5.2)') '     Mixed central moment 1,1             of dat1 and dat2 (target=11.47): ', &
       mixed_central_moment(dat1(1:9),dat2(1:9),1,1)
  Write(*,'(A,f5.2)') '     Variance of mixed central moment 1,1 of dat1 and dat2 (target=11.29): ', &
       mixed_central_moment_var(dat1(1:9),dat2(1:9),1,1)
  Write(*,*) '    With 10 masked'
  Write(*,'(A,f5.2)') '     Correlation                          of dat1 and dat2 (target=0.88):  ', &
       correlation(dat1,dat2, mask=((dat1 /= 10.) .and. (dat2 /= 10.)))
  Write(*,'(A,f5.2)') '     Covariance                           of dat1 and dat2 (target=11.47): ', &
       covariance(dat1,dat2, mask=((dat1 /= 10.) .and. (dat2 /= 10.)))
  Write(*,'(A,f5.2)') '     Mixed central moment 1,1             of dat1 and dat2 (target=11.47): ', &
       mixed_central_moment(dat1,dat2,1,1, mask=((dat1 /= 10.) .and. (dat2 /= 10.)))
  Write(*,'(A,f5.2)') '     Variance of mixed central moment 1,1 of dat1 and dat2 (target=11.29): ', &
       mixed_central_moment_var(dat1,dat2,1,1, mask=((dat1 /= 10.) .and. (dat2 /= 10.)))
  

  ! Single precision
  Write(*,*) ''
  Write(*,*) '  Test single precision'

  ! Test the single input functions
  sat1 = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  Write(*,'(A,10f5.1)') '     sat1 is: ', sat1
  Write(*,'(A,f5.1)') '     Mean       of sat1 (target=5.5):  ', mean(sat1)
  Write(*,'(A,f5.1)') '     Average    of sat1 (target=5.5):  ', average(sat1)
  Write(*,'(A,f5.1)') '     Variance   of sat1 (target=9.2):  ', variance(sat1)
  Write(*,'(A,f5.1)') '     Stddev     of sat1 (target=3.0):  ', stddev(sat1)
  Write(*,'(A,f5.1)') '     Skewness   of sat1 (target=0.0):  ', skewness(sat1)
  Write(*,'(A,f5.1)') '     Kurtosis   of sat1 (target=-1.6): ', kurtosis(sat1)
  Write(*,'(A,f5.1)') '     MeanAbsDev of sat1 (target=2.5):  ', absdev(sat1)
  Write(*,'(A,f5.1)') '     2nd Central Moment             of sat1 (target=8.2): ', central_moment(sat1,2)
  Write(*,'(A,f5.1)') '     Variance of 2nd central moment of sat1 (target=5.3): ', central_moment_var(sat1,2)
  ! Test the mask
  Write(*,*) '    Without last'
  Write(*,'(A,f5.1)') '     Mean       of sat1 (target=5.0):  ', mean(sat1(1:9))
  Write(*,'(A,f5.1)') '     Average    of sat1 (target=5.0):  ', average(sat1(1:9))
  Write(*,'(A,f5.1)') '     Variance   of sat1 (target=7.5):  ', variance(sat1(1:9))
  Write(*,'(A,f5.1)') '     Stddev     of sat1 (target=2.7):  ', stddev(sat1(1:9))
  Write(*,'(A,f5.1)') '     Skewness   of sat1 (target=0.0):  ', skewness(sat1(1:9))
  Write(*,'(A,f5.1)') '     Kurtosis   of sat1 (target=-1.6): ', kurtosis(sat1(1:9))
  Write(*,'(A,f5.1)') '     MeanAbsDev of sat1 (target=2.2):  ', absdev(sat1(1:9))
  Write(*,'(A,f5.1)') '     2nd Central Moment             of sat1 (target=6.7): ', central_moment(sat1(1:9),2)
  Write(*,'(A,f5.1)') '     Variance of 2nd central moment of sat1 (target=3.8): ', central_moment_var(sat1(1:9),2)
  Write(*,*) '    With 10 masked'
  Write(*,'(A,f5.1)') '     Mean       of sat1 (target=5.0):  ', mean(sat1, mask=(sat1 /= 10.))
  Write(*,'(A,f5.1)') '     Average    of sat1 (target=5.5):  ', average(sat1, mask=(sat1 /= 10.))
  Write(*,'(A,f5.1)') '     Variance   of sat1 (target=7.5):  ', variance(sat1, mask=(sat1 /= 10.))
  Write(*,'(A,f5.1)') '     Stddev     of sat1 (target=2.7):  ', stddev(sat1, mask=(sat1 /= 10.))
  Write(*,'(A,f5.1)') '     Skewness   of sat1 (target=0.0):  ', skewness(sat1, mask=(sat1 /= 10.))
  Write(*,'(A,f5.1)') '     Kurtosis   of sat1 (target=-1.6): ', kurtosis(sat1, mask=(sat1 /= 10.))
  Write(*,'(A,f5.1)') '     MeanAbsDev of sat1 (target=2.2):  ', absdev(sat1, mask=(sat1 /= 10.))
  Write(*,'(A,f5.1)') '     2nd Central Moment             of sat1 (target=6.7): ', &
       central_moment(sat1,2, mask=(sat1 /= 10.))
  Write(*,'(A,f5.1)') '     Variance of 2nd central moment of sat1 (target=3.8): ', &
       central_moment_var(sat1,2, mask=(sat1 /= 10.))
  call moment(sat1, average=a_sp, variance=v_sp, skewness=s_sp, kurtosis=k_sp, mean=m_sp, &
       stddev=std_sp, absdev=ad_sp, mask=(sat1 /= 10.))
  Write(*,'(A,f5.1,1X,f5.1,1X,f5.1,1X,f5.1,1X,f5.1,1X,f5.1,1X,f5.1,1X)') '     Moments of a are: ', &
       m_sp, a_sp, v_sp, std_sp, s_sp, k_sp, ad_sp
  ! 
  ! Test the double input functions
  Write(*,*) ''
  sat2 = (/ 1.3, 2.7, 3.9, 5.5, 7., 8.8, 11., 13., 15., 16. /)
  Write(*,'(A,10f5.1)') '     sat2 is: ', sat2
  Write(*,'(A,f5.2)') '     Correlation                          of sat1 and sat2 (target=0.90):  ', &
       correlation(sat1,sat2)
  Write(*,'(A,f5.2)') '     Covariance                           of sat1 and sat2 (target=14.11): ', &
       covariance(sat1,sat2)
  Write(*,'(A,f5.2)') '     Mixed central moment 1,1             of sat1 and sat2 (target=14.11): ', &
       mixed_central_moment(sat1,sat2,1,1)
  Write(*,'(A,f5.2)') '     Variance of mixed central moment 1,1 of sat1 and sat2 (target=14.24): ', &
       mixed_central_moment_var(sat1,sat2,1,1)
  Write(*,*) '    Without last'
  Write(*,'(A,f5.2)') '     Correlation                          of sat1 and sat2 (target=0.88):  ', &
       correlation(sat1(1:9),sat2(1:9))
  Write(*,'(A,f5.2)') '     Covariance                           of sat1 and sat2 (target=11.47): ', &
       covariance(sat1(1:9),sat2(1:9))
  Write(*,'(A,f5.2)') '     Mixed central moment 1,1             of sat1 and sat2 (target=11.47): ', &
       mixed_central_moment(sat1(1:9),sat2(1:9),1,1)
  Write(*,'(A,f5.2)') '     Variance of mixed central moment 1,1 of sat1 and sat2 (target=11.29): ', &
       mixed_central_moment_var(sat1(1:9),sat2(1:9),1,1)
  Write(*,*) '    With 10 masked'
  Write(*,'(A,f5.2)') '     Correlation                          of sat1 and sat2 (target=0.88):  ', &
       correlation(sat1,sat2, mask=((sat1 /= 10.) .and. (sat2 /= 10.)))
  Write(*,'(A,f5.2)') '     Covariance                           of sat1 and sat2 (target=11.47): ', &
       covariance(sat1,sat2, mask=((sat1 /= 10.) .and. (sat2 /= 10.)))
  Write(*,'(A,f5.2)') '     Mixed central moment 1,1             of sat1 and sat2 (target=11.47): ', &
       mixed_central_moment(sat1,sat2,1,1, mask=((sat1 /= 10.) .and. (sat2 /= 10.)))
  Write(*,'(A,f5.2)') '     Variance of mixed central moment 1,1 of sat1 and sat2 (target=11.29): ', &
       mixed_central_moment_var(sat1,sat2,1,1, mask=((sat1 /= 10.) .and. (sat2 /= 10.)))


END PROGRAM main
