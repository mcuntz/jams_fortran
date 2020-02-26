PROGRAM main

  USE mo_kind,          ONLY: i4, dp, sp
  USE mo_utils,         ONLY: ne
  USE mo_distributions, ONLY: laplace, laplace01, normal, normal01
  USE mo_distributions, ONLY: ep, ep01, sep, sep01
  USE mo_distributions, ONLY: st, st01, t, t01
  use mo_ansi_colors, only: color, c_red, c_green

  IMPLICIT NONE

  ! data
  INTEGER(i4), PARAMETER  :: nn = 10000

  REAL(dp)                :: out0, iloc, isig, ddat
  REAL(dp), DIMENSION(nn) :: dat, out
  REAL(dp)                :: loc, sca, sig, xi, beta, nu
  REAL(dp)                :: dmean, dvar, dstd

  REAL(sp)                :: isloc, issig, dsat
  REAL(sp), DIMENSION(nn) :: sat, sout
  REAL(sp)                :: sloc, ssca, ssig, sxi, sbeta, snu
  REAL(dp)                :: smean, svar, sstd

  INTEGER(i4) :: i
  LOGICAL     :: isgood, allgood

  Write(*,*) ''
  Write(*,*) 'Test mo_distributions.f90'

  allgood = .true.

  ! DP
  forall(i=1:nn) dat(i) = real(i,dp)/real(nn,dp)*300.0_dp - 150.0_dp
  ddat = dat(2)-dat(1)
  
  loc  = 1.1_dp
  sca  = 2.2_dp
  sig  = 2.2_dp
  xi   = 3.3_dp
  beta = 0.5_dp
  nu   = 4.4_dp

  ! integral
  do i=1, 12
     isgood = .true.
     select case(i)
     case(1) ! Exponential Power
        iloc = loc
        isig = sig
        out = ep(dat, loc, beta=beta, sig=sig)
     case(2)
        iloc = 0.0_dp
        isig = 1.0_dp
        out = ep01(dat, beta)
     case(3)! Laplace
        iloc = loc
        isig = sig
        out = laplace(dat, loc, sig=sig)
     case(4)
        iloc = 0.0_dp
        isig = sqrt(2.0_dp)
        out = laplace01(dat)
     case(5) ! Normal
        iloc = loc
        isig = sig
        out = normal(dat, loc, sig=sig)
     case(6)
        iloc = 0.0_dp
        isig = 1.0_dp
        out = normal01(dat)
     case(7) ! Skew Exponential Power
        iloc = loc
        isig = sig
        out = sep(dat, loc, xi=xi, beta=beta, sig=sig)
     case(8)
        iloc = 0.0_dp
        isig = 1.0_dp
        out = sep01(dat, xi, beta)
     case(9) ! Skewed Student t
        iloc = loc
        isig = sig
        out = st(dat, nu, loc, xi=xi, sig=sig)
     case(10)
        iloc = 0.0_dp
        isig = sqrt(nu/(nu-2.0_dp))
        out = st01(dat, nu, xi)
     case(11) ! Student t
        iloc = loc
        isig = sig
        out = t(dat, nu, loc, sig=sig)
     case(12)
        iloc = 0.0_dp
        isig = sqrt(nu/(nu-2.0_dp))
        out = t01(dat, nu)
     case default
        continue
     end select
     dmean = sum(dat * out * ddat)
     dvar  = sum((dat-dmean)**2 * out * ddat)
     dstd  = sqrt(dvar)
     if (abs(dmean-iloc) > 1.0e-2_dp) isgood = .false.
     if (abs(dstd-isig) > 1.0e-2_dp) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions integral DP failed: ', i, iloc, dmean, isig, dstd
  enddo

  ! normal vs normal01
  do i=1, 6
     isgood = .true.
     select case(i)
     case(1)
        out = ep(dat, beta=beta) - ep01(dat, beta)
     case(2)
        out = laplace(dat, 0.0_dp) - laplace01(dat)
     case(3)
        out = normal(dat, 0.0_dp, 1.0_dp) - normal01(dat)
     case(4)
        out = sep(dat, 0.0_dp, 1.0_dp, xi, beta) - sep01(dat, xi, beta)
     case(5)
        out = st(dat, nu, xi=xi) - st01(dat, nu, xi)
     case(6)
        out = t(dat, nu, sca=1.0_dp) - t01(dat, nu)
     case default
        continue
     end select
     if (any(ne(out, 0.0_dp))) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions comparison 01 versions failed: ', i, out
  enddo

  ! transformation
  do i=1, 6
     isgood = .true.
     select case(i)
     case(1)
        out = ep(dat, loc, sca, beta) - ep((dat-loc)/sca, beta=beta)/sca
     case(2)
        out = laplace(dat, loc, sca) - laplace((dat-loc)/sca)/sca
     case(3)
        out = normal(dat, loc, sca) - normal((dat-loc)/sca)/sca
     case(4)
        out = sep(dat, loc, sca, xi=xi, beta=beta) - sep((dat-loc)/sca, xi=xi, beta=beta)/sca
     case(5)
        out = st(dat, nu, loc, sca, xi=xi) - st((dat-loc)/sca, nu, xi=xi)/sca
     case(6)
        out = t(dat, nu, loc, sca) - t((dat-loc)/sca, nu)/sca
     case default
        continue
     end select
     ! if (any(ne(out, 0.0_dp))) isgood = .false.
     if (any(abs(out) > 1.0e-15_dp)) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions transformation 1D failed: ', i, maxval(abs(out))
  enddo

  ! transformation 0D
  do i=1, 6
     isgood = .true.
     select case(i)
     case(1)
        out0 = ep(0.5_dp, loc, sca, beta=beta) - ep((0.5_dp-loc)/sca, beta=beta)/sca
     case(2)
        out0 = laplace(0.5_dp, loc, sca) - laplace((0.5_dp-loc)/sca)/sca
     case(3)
        out0 = normal(0.5_dp, loc, sca) - normal((0.5_dp-loc)/sca)/sca
     case(4)
        out0 = sep(0.5_dp, loc, sca, xi=xi, beta=beta) - sep((0.5_dp-loc)/sca, xi=xi, beta=beta)/sca
     case(5)
        out0 = st(0.5_dp, nu, loc, sca, xi=xi) - st((0.5_dp-loc)/sca, nu, xi=xi)/sca
     case(6)
        out0 = t(0.5_dp, nu, loc, sca) - t((0.5_dp-loc)/sca, nu)/sca
     case default
        continue
     end select
     ! if (ne(out0, 0.0_dp)) isgood = .false.
     if (abs(out0) > 1.0e-15_dp) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions transformation 0D failed: ', i, out0
  enddo

  ! ep/sep vs normal/laplace
  do i=1, 9
     isgood = .true.
     select case(i)
     case(1)
        out = ep(dat, 0.0_dp, 1.0_dp, 0.0_dp) - normal(dat, 0.0_dp, 1.0_dp)
     case(2)
        out = ep01(dat, beta=0.0_dp) - normal01(dat)
     case(3)
        out = sep(dat, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp) - normal(dat, 0.0_dp, 1.0_dp)
     case(4)
        out = sep01(dat) - normal01(dat)
     case(5)
        out = ep(dat, 0.0_dp, 1.0_dp, 1.0_dp) - laplace(dat, 0.0_dp, sig=1.0_dp)
     case(6)
        out = ep(dat, 0.0_dp, 1.0_dp, 1.0_dp) - laplace(dat, 0.0_dp, 1.0_dp/sqrt(2.0_dp))
     case(7)
        out = ep(dat, 0.0_dp, sqrt(2.0_dp), 1.0_dp) - laplace(dat, 0.0_dp, 1.0_dp)
     case(8)
        out = sep(dat, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp) - laplace(dat, 0.0_dp, sig=1.0_dp)
     case(9)
        out = sep(dat, 0.0_dp, sqrt(2.0_dp), 1.0_dp, 1.0_dp) - laplace(dat, 0.0_dp, 1.0_dp)
     case default
        continue
     end select
     if (any(abs(out) > 1.0e-10)) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions ep/sep comparison to normal/laplace failed: ', i, out
  enddo

  ! skewed vs unskewed
  do i=1, 2
     isgood = .true.
     select case(i)
     case(1)
        out = sep(dat, loc, sca, 1.0_dp, beta) - ep(dat, loc, sca, beta)
     case(2)
        out = st(dat, nu, loc, sca, 1.0_dp) - t(dat, nu, loc, sca)
     case default
        continue
     end select
     if (any(abs(out) > 1.0e-10)) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions skewed vs non-skewed failed: ', i, out
  enddo


  ! SP
  forall(i=1:nn) sat(i) = real(i,sp)/real(nn,sp)*300.0_sp - 150.0_sp
  dsat = sat(2)-sat(1)
  
  sloc  = 1.1_sp
  ssca  = 2.2_sp
  ssig  = 2.2_sp
  sxi   = 3.3_sp
  sbeta = 0.5_sp
  snu   = 4.4_sp

  ! integral
  do i=1, 12
     isgood = .true.
     select case(i)
     case(1) ! Exponential Power
        isloc = sloc
        issig = ssig
        sout = ep(sat, sloc, beta=sbeta, sig=ssig)
     case(2)
        isloc = 0.0_sp
        issig = 1.0_sp
        sout = ep01(sat, sbeta)
     case(3)! Laplace
        isloc = sloc
        issig = ssig
        sout = laplace(sat, sloc, sig=ssig)
     case(4)
        isloc = 0.0_sp
        issig = sqrt(2.0_sp)
        sout = laplace01(sat)
     case(5) ! Normal
        isloc = sloc
        issig = ssig
        sout = normal(sat, sloc, ssca)
     case(6)
        isloc = 0.0_sp
        issig = 1.0_sp
        sout = normal01(sat)
     case(7) ! Skew Exponential Power
        isloc = sloc
        issig = ssig
        sout = sep(sat, sloc, xi=sxi, beta=sbeta, sig=ssig)
     case(8)
        isloc = 0.0_sp
        issig = 1.0_sp
        sout = sep01(sat, sxi, sbeta)
     case(9) ! Skewed Student t
        isloc = sloc
        issig = ssig
        sout = st(sat, snu, sloc, xi=sxi, sig=ssig)
     case(10)
        isloc = 0.0_sp
        issig = sqrt(snu/(snu-2.0_sp))
        sout = st01(sat, snu, sxi)
     case(11) ! Student t
        isloc = sloc
        issig = ssig
        sout = t(sat, snu, sloc, sig=ssig)
     case(12)
        isloc = 0.0_sp
        issig = sqrt(snu/(snu-2.0_sp))
        sout = t01(sat, snu)
     case default
        continue
     end select
     smean = sum(sat * sout * dsat)
     svar  = sum((sat-smean)**2 * sout * dsat)
     sstd  = sqrt(svar)
     if (abs(smean-isloc) > 1.0e-2_sp) isgood = .false.
     if (abs(sstd-issig) > 1.0e-2_sp) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions integral SP failed: ', i, isloc, dmean, issig, dstd
  enddo


  ! Finish
  allgood = allgood .and. isgood
  if (allgood) then
     write(*,*) 'mo_distributions ', color('o.k.', c_green)
  else
     write(*,*) 'mo_distributions ', color('failed!', c_red)
  endif

END PROGRAM main
