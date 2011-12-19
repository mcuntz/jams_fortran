PROGRAM main
  
  USE mo_kind, ONLY: dp, sp
  USE mo_fit,  ONLY: fitfun, fpoly_dp, fpoly_sp, polyfit, svdfit
  USE mo_fit,  ONLY: linfit

  IMPLICIT NONE
  
  LOGICAL :: isgood
  
  integer, parameter      :: degree = 2
  integer                 :: i

  real(dp), dimension(11) :: dx = (/ (i,i=0,10) /)
  real(dp), dimension(11) :: dy = (/ 1,   6,  17,  34, 57,  86, 121, 162, 209, 262, 321 /)
  real(dp), dimension(11) :: dsig, dout
  real(dp), dimension(degree+1) :: da, dw
  real(dp), dimension(degree+1,degree+1) :: dv
  real(dp) :: dchisq
  real(dp), dimension(2) :: dla
  real(dp) :: da1, db1

  real(sp), dimension(11) :: sx = (/ (i,i=0,10) /)
  real(sp), dimension(11) :: sy = (/ 1,   6,  17,  34, 57,  86, 121, 162, 209, 262, 321 /)
  real(sp), dimension(11) :: ssig, sout
  real(sp), dimension(degree+1) :: sa, sw
  real(sp), dimension(degree+1,degree+1) :: sv
  real(sp) :: schisq
  real(sp), dimension(2) :: sla
  real(sp) :: sa1, sb1

  Write(*,*) ''
  Write(*,*) 'Test mo_fit.f90'
  Write(*,*) ''
  Write(*,*) 'Output should read 1, 2, 3'

  ! Double precision
  Write(*,*) ''
  isgood = .true.
  da = polyfit(dx, dy, degree)
  write(*,*) da
  isgood = isgood .and. all((anint(1000._dp*da)-(/1000._dp,2000._dp,3000._dp/)) == 0._dp)
  dsig = 0.1
  call fitfun(dx, dy, dsig, da, fpoly_dp)
  write(*,*) da
  isgood = isgood .and. all((anint(1000._dp*da)-(/1000._dp,2000._dp,3000._dp/)) == 0._dp)
  call svdfit(dx, dy, dsig, da, dv, dw, dchisq, fpoly_dp)
  write(*,*) da
  isgood = isgood .and. all((anint(1000._dp*da)-(/1000._dp,2000._dp,3000._dp/)) == 0._dp)
  Write(*,*) ''
  Write(*,*) 'Output should read -44, 32'
  Write(*,*) ''
  dla = polyfit(dx, dy, 1)
  write(*,*) dla
  dout = linfit(dx, dy, a=da1, b=db1)
  write(*,*) da1, db1
  isgood = isgood .and. all(anint((dla-(/da1, db1/))*1000._dp) == 0._dp)
  

  if (isgood) then
     write(*,*) 'mo_fit double precision o.k.'
  else
     write(*,*) 'mo_fit double precision failed!'
  endif

  ! Single precision
  Write(*,*) ''
  isgood = .true.
  sa = polyfit(sx, sy, degree)
  write(*,*) sa
  isgood = isgood .and. all((anint(1000._sp*sa)-(/1000._sp,2000._sp,3000._sp/)) == 0._sp)
  ssig = 0.1
  call fitfun(sx, sy, ssig, sa, fpoly_sp)
  write(*,*) sa
  isgood = isgood .and. all((anint(1000._sp*sa)-(/1000._sp,2000._sp,3000._sp/)) == 0._sp)
  call svdfit(sx, sy, ssig, sa, sv, sw, schisq, fpoly_sp)
  write(*,*) sa
  isgood = isgood .and. all((anint(1000._sp*sa)-(/1000._sp,2000._sp,3000._sp/)) == 0._sp)
  Write(*,*) ''
  Write(*,*) 'Output should read -44, 32'
  Write(*,*) ''
  sla = polyfit(sx, sy, 1)
  write(*,*) sla
  sout = linfit(sx, sy, a=sa1, b=sb1)
  write(*,*) sa1, sb1
  isgood = isgood .and. all(anint((sla-(/sa1, sb1/))*1000._sp) == 0._sp)

  if (isgood) then
     write(*,*) 'mo_fit single precision o.k.'
  else
     write(*,*) 'mo_fit single precision failed!'
  endif

END PROGRAM main
