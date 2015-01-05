PROGRAM main
  
  USE mo_kind,      ONLY: i4, dp, sp
  USE mo_linfit,    ONLY: linfit
  
  IMPLICIT NONE
  
  ! data
  INTEGER(i4), PARAMETER :: nn = 100

  REAL(dp), DIMENSION(nn) :: xx, yy, yout
  REAL(dp) :: a, b, a2, b2, siga, sigb, chi2

  REAL(sp), DIMENSION(nn) :: sxx, syy, syout
  REAL(sp) :: sa, sb, sa2, sb2, ssiga, ssigb, schi2

  LOGICAL  :: model2
  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_linfit.f90'

  isgood = .true.

  ! Double precision
  a = 2._dp
  b = 3._dp
  call random_number(xx)
  yy = a + b*xx

  ! model I
  model2 = .false.
  ! fit should be perfect
  yout = linfit(xx, yy, a=a2, b=b2, siga=siga, sigb=sigb, chi2=chi2, model2=model2)
  yout = yout
  if (abs(a-a2) > 100._dp*epsilon(1.0_dp)) isgood = .false.
  if (siga > 100._dp*epsilon(1.0_dp)) isgood = .false.
  if (sigb > 100._dp*epsilon(1.0_dp)) isgood = .false.
  if (chi2 > epsilon(1.0_dp)) isgood = .false.
  
  ! model II
  model2 = .true.
  yout = linfit(xx, yy, a=a2, b=b2, siga=siga, sigb=sigb, chi2=chi2, model2=model2)
  yout = yout
  if (abs(a-a2) > 100._dp*epsilon(1.0_dp)) isgood = .false.
  if (siga > 1.0e-4*a) isgood = .false.
  if (sigb > 1.0e-4*b) isgood = .false.
  if (chi2 > epsilon(1.0_dp)) isgood = .false.


  ! Single precision
  sa = 2._sp
  sb = 3._sp
  call random_number(sxx)
  syy = sa + sb*sxx

  ! model I
  model2 = .false.
  ! fit should be perfect
  syout = linfit(sxx, syy, a=sa2, b=sb2, siga=ssiga, sigb=ssigb, chi2=schi2, model2=model2)
  syout = syout
  if (abs(sa-sa2) > 100._sp*epsilon(1.0_sp)) isgood = .false.
  if (ssiga > 100._sp*epsilon(1.0_sp)) isgood = .false.
  if (ssigb > 100._sp*epsilon(1.0_sp)) isgood = .false.
  if (schi2 > epsilon(1.0_sp)) isgood = .false.
  
  ! model II
  model2 = .true.
  syout = linfit(sxx, syy, a=sa2, b=sb2, siga=ssiga, sigb=ssigb, chi2=schi2, model2=model2)
  syout = syout
  if (abs(sa-sa2) > 100._sp*epsilon(1.0_sp)) isgood = .false.
  if (ssiga > 1.0e-4*sa) isgood = .false.
  if (ssigb > 1.0e-4*sb) isgood = .false.
  if (schi2 > epsilon(1.0_sp)) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_linfit o.k.'
  else
     write(*,*) 'mo_linfit failed!'
  endif

END PROGRAM main
