PROGRAM main
  
  USE mo_kind, ONLY: dp, sp
  USE mo_fit,  ONLY: fitfun, fpoly_dp, fpoly_sp, polyfit, svdfit
  USE mo_fit,  ONLY: polyval
  use mo_ansi_colors, only: color, c_red, c_green

  IMPLICIT NONE
  
  LOGICAL                                :: isgood
  
  integer, parameter                     :: degree = 2
  integer                                :: i

  real(dp), dimension(11)                :: dx = (/ (i,i=0,10) /)
  real(dp), dimension(11)                :: dy = (/ 1,   6,  17,  34, 57,  86, 121, 162, 209, 262, 321 /)
  real(dp), dimension(11)                :: dsig
  real(dp), dimension(degree+1)          :: da, dw
  real(dp), dimension(degree+1,degree+1) :: dv
  real(dp)                               :: dchisq
  real(dp), dimension(2)                 :: dla
  real(dp)                               :: xdp = 3
  real(dp), dimension(4)                 :: xdpvec = (/ 0, 1, 2, 3 /)
  real(dp), dimension(3)                 :: polydp = (/ 1, 2, 1 /)

  real(sp), dimension(11)                :: sx = (/ (i,i=0,10) /)
  real(sp), dimension(11)                :: sy = (/ 1,   6,  17,  34, 57,  86, 121, 162, 209, 262, 321 /)
  real(sp), dimension(11)                :: ssig
  real(sp), dimension(degree+1)          :: sa, sw
  real(sp), dimension(degree+1,degree+1) :: sv
  real(sp)                               :: schisq
  real(sp), dimension(2)                 :: sla
  real(sp)                               :: xsp = 3
  real(sp), dimension(4)                 :: xspvec = (/ 0, 1, 2, 3 /)
  real(sp), dimension(3)                 :: polysp = (/ 1, 2, 1 /)

  Write(*,*) ''
  Write(*,*) 'Test mo_fit.f90'
  Write(*,*) ''
  Write(*,*) 'Output should read 1, 2, 3'

  ! Double precision
  Write(*,*) ''
  isgood = .true.
  da = polyfit(dx, dy, degree)
  write(*,*) da
  isgood = isgood .and. all((nint(1000._dp*da)-(/1000, 2000, 3000/)) == 0)
  dsig = 0.1
  call fitfun(dx, dy, dsig, da, fpoly_dp)
  write(*,*) da
  isgood = isgood .and. all((nint(1000._dp*da)-(/1000, 2000, 3000/)) == 0)
  call svdfit(dx, dy, dsig, da, dv, dw, dchisq, fpoly_dp)
  write(*,*) da
  isgood = isgood .and. all((nint(1000._dp*da)-(/1000, 2000, 3000/)) == 0)
  Write(*,*) ''
  Write(*,*) 'Output should read -44, 32'
  Write(*,*) ''
  dla = polyfit(dx, dy, 1)
  write(*,*) dla
  Write(*,*) ''
  Write(*,*) 'Output should read 16 & 1 4 9 16'
  Write(*,*) ''
  write(*,*) polyval( polydp, xdp)
  write(*,*) polyval( polydp, xdpvec)
  isgood = isgood .and. all((nint(1000._dp*polyval( polydp, xdpvec))-(/1000, 4000, 9000, 16000/)) == 0)

  if (isgood) then
     write(*,*) 'mo_fit double precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_fit double precision ', color('failed!', c_red)
  endif

  ! Single precision
  Write(*,*) ''
  isgood = .true.
  sa = polyfit(sx, sy, degree)
  write(*,*) sa
  isgood = isgood .and. all((nint(1000._sp*sa)-(/1000, 2000, 3000/)) == 0)
  ssig = 0.1
  call fitfun(sx, sy, ssig, sa, fpoly_sp)
  write(*,*) sa
  isgood = isgood .and. all((nint(1000._sp*sa)-(/1000, 2000, 3000/)) == 0)
  call svdfit(sx, sy, ssig, sa, sv, sw, schisq, fpoly_sp)
  write(*,*) sa
  isgood = isgood .and. all((nint(1000._sp*sa)-(/1000, 2000, 3000/)) == 0)
  Write(*,*) ''
  Write(*,*) 'Output should read -44, 32'
  Write(*,*) ''
  sla = polyfit(sx, sy, 1)
  write(*,*) sla
  Write(*,*) ''
  Write(*,*) 'Output should read 16 & 1 4 9 16'
  Write(*,*) ''
  write(*,*) polyval( polysp, xsp)
  write(*,*) polyval( polysp, xspvec)
  isgood = isgood .and. all((nint(1000._sp*polyval( polysp, xspvec))-(/1000, 4000, 9000, 16000/)) == 0)

  if (isgood) then
     write(*,*) 'mo_fit single precision ', color('o.k.', c_green)
  else
     write(*,*) 'mo_fit single precision ', color('failed!', c_red)
  endif

END PROGRAM main
