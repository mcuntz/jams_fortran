PROGRAM main
  
  USE mo_kind,   ONLY: dp, sp
  USE mo_corr,   ONLY: autocoeffk, autocorr, corr, crosscoeffk, crosscorr
  USE mo_moment, ONLY: mean, variance, covariance

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: n = 1500
  REAL(dp), DIMENSION(n) :: dat1, dat2, dc
  REAL(sp), DIMENSION(n) :: sat1, sat2, sc
  INTEGER :: i, nout

  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_corr.f90'

  ! Double precision
  forall(i=1:n) dat1(i) = real(i,dp)
  dat2 = dat1 + sin(dat1)
  isgood = .true.
  if (abs(variance(dat2)*real(n-1,dp)/real(n,dp) - autocoeffk(dat2,0)) > epsilon(1.0_dp)) &
       isgood =.false.
  if (abs(autocorr(dat2,0) - 1.0_dp) > epsilon(1.0_dp)) isgood =.false.
  if (abs(covariance(dat1,dat2) - crosscoeffk(dat1,dat2,0))> epsilon(1.0_dp)) isgood =.false.
  if (abs(autocorr(dat1,10) - crosscorr(dat1,dat1,10))> epsilon(1.0_dp)) isgood =.false.
  ! Covariance function: covariance(x,y) = corr(x,y)/n
  call random_number(dat1)
  dat2 = dat1 + sin(dat1)
  dc = corr(dat1-mean(dat1), dat2-mean(dat2), nadjust=nout)
  if (abs(dc(1)/real(nout,dp) - covariance(dat1,dat2)) > 0.01_dp) isgood =.false.

  if (isgood) then
     write(*,*) 'mo_corr double precision o.k.'
  else
     write(*,*) 'mo_corr double precision failed!'
  endif

  ! Single precision
  forall(i=1:n) sat1(i) = real(i,sp)
  sat2 = sat1 + sin(sat1)
  isgood = .true.
  if (abs(variance(sat2)*real(n-1,sp)/real(n,sp) - autocoeffk(sat2,0)) > epsilon(1.0_sp)) &
       isgood =.false.
  if (abs(autocorr(sat2,0) - 1.0_sp) > epsilon(1.0_sp)) isgood =.false.
  if (abs(covariance(sat1,sat2) - crosscoeffk(sat1,sat2,0))> epsilon(1.0_sp)) isgood =.false.
  if (abs(autocorr(sat1,10) - crosscorr(sat1,sat1,10))> epsilon(1.0_sp)) isgood =.false.
  ! Covariance function: covariance(x,y) = corr(x,y)/n
  call random_number(sat1)
  sat2 = sat1 + sin(sat1)
  sc = corr(sat1-mean(sat1), sat2-mean(sat2), nadjust=nout)
  if (abs(sc(1)/real(nout,sp) - covariance(sat1,sat2)) > 0.01_sp) isgood =.false.

  if (isgood) then
     write(*,*) 'mo_corr single precision o.k.'
  else
     write(*,*) 'mo_corr single precision failed!'
  endif

END PROGRAM main
