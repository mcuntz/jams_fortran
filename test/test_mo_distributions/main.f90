PROGRAM main

  USE mo_kind,          ONLY: i4, dp, sp
  USE mo_utils,         ONLY: ne
  USE mo_distributions, ONLY: laplace, laplace01, normal, normal01
  USE mo_distributions, ONLY: ep, ep01, sep, sep01
  USE mo_distributions, ONLY: st, st01, t, t01
  USE mo_distributions, ONLY: beta_den, gamma_pdf, weibull_pdf
  USE mo_functions,     ONLY: gamm
  use mo_ansi_colors,   only: color, c_red, c_green

  IMPLICIT NONE

  ! data
  INTEGER(i4), PARAMETER  :: nn = 10000

  REAL(dp)                :: out0, iloc, isig, ddat, ddat_beta, ddat_gamma, ddat_weibull
  REAL(dp), DIMENSION(nn) :: dat, dat_beta, dat_gamma, dat_weibull, out
  REAL(dp)                :: loc, sca, sig, xi, beta, nu
  REAL(dp)                :: dmean, dvar, dstd

  REAL(sp)                :: isloc, issig, dsat, dsat_beta, dsat_gamma, dsat_weibull
  REAL(sp), DIMENSION(nn) :: sat, sat_beta, sat_gamma, sat_weibull, sout
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
  forall(i=1:nn) dat_beta(i) = real(i,dp)/real(nn,dp)
  ddat_beta = dat_beta(2)-dat_beta(1)
  forall(i=1:nn) dat_gamma(i) = real(i,dp)/real(nn,dp) * 20._dp
  ddat_gamma = dat_gamma(2)-dat_gamma(1)
  forall(i=1:nn) dat_weibull(i) = real(i,dp)/real(nn,dp) * 5._dp
  ddat_weibull = dat_weibull(2)-dat_weibull(1)

  loc  = 1.1_dp
  sca  = 2.2_dp
  sig  = 2.2_dp
  xi   = 3.3_dp
  beta = 0.5_dp
  nu   = 4.4_dp

  ! integral
  do i=1, 15
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
     case(13) ! beta distribution
        iloc = 0.5_dp
        isig = 0.1507556722888821_dp
        out = beta_den(dat_beta, a=5._dp, b=5._dp, l=0._dp, s=1._dp)
     case(14) ! gamma distribution
        iloc = 0.75_dp * 2._dp
        isig = sqrt(0.75_dp * 2._dp**2._dp)
        out = gamma_pdf(dat_gamma, c=0.75_dp, l=0._dp, s=2._dp)
     case(15) ! weibull distribution
        iloc = 1._dp * gamm(1._dp + 1._dp / 5._dp)
        isig = 1._dp * sqrt(gamm(1._dp + 2._dp / 5._dp) - gamm(1._dp + 1._dp / 5._dp)**2._dp )
        out = weibull_pdf(dat_weibull, c=5._dp, l=0._dp, s=1._dp)
     case default
        out = 0.0_dp
     end select
     if (i .eq. 13) then
        dmean = sum(dat_beta * out * ddat_beta)
        dvar  = sum((dat_beta-dmean)**2 * out * ddat_beta)
     else if (i .eq. 14) then
        dmean = sum(dat_gamma * out * ddat_gamma)
        dvar  = sum((dat_gamma-dmean)**2 * out * ddat_gamma)
     else if (i .eq. 15) then
        dmean = sum(dat_weibull * out * ddat_weibull)
        dvar  = sum((dat_weibull-dmean)**2 * out * ddat_weibull)
     else 
        dmean = sum(dat * out * ddat)
        dvar  = sum((dat-dmean)**2 * out * ddat)
     end if
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
        out = 0.0_dp
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
        out = 0.0_dp
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
        out0 = 0.0_dp
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
        out = 0.0_dp
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
        out = 0.0_dp
     end select
     if (any(abs(out) > 1.0e-10)) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions skewed vs non-skewed failed: ', i, out
  enddo


  ! SP
  forall(i=1:nn) sat(i) = real(i,sp)/real(nn,sp)*300.0_sp - 150.0_sp
  dsat = sat(2)-sat(1)
  forall(i=1:nn) sat_beta(i) = real(i,sp)/real(nn,sp)
  dsat_beta = sat_beta(2)-sat_beta(1)
  forall(i=1:nn) sat_gamma(i) = real(i,sp)/real(nn,sp) * 20._sp
  dsat_gamma = sat_gamma(2)-sat_gamma(1)
  forall(i=1:nn) sat_weibull(i) = real(i,sp)/real(nn,sp) * 5._sp
  dsat_weibull = sat_weibull(2)-sat_weibull(1)
  
  sloc  = 1.1_sp
  ssca  = 2.2_sp
  ssig  = 2.2_sp
  sxi   = 3.3_sp
  sbeta = 0.5_sp
  snu   = 4.4_sp

  ! integral
  do i=1, 15
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
     case(13) ! beta distribution
        isloc = 0.5_sp
        issig = 0.1507556722888821_sp
        sout = beta_den(sat_beta, a=5._sp, b=5._sp, l=0._sp, s=1._sp)
     case(14) ! gamma distribution
        isloc = 0.75_sp * 2._sp
        issig = sqrt(0.75_sp * 2._sp**2._sp)
        sout = gamma_pdf(sat_gamma, c=0.75_sp, l=0._sp, s=2._sp)
     case(15) ! weibull distribution
        isloc = 1._sp * gamm(1._sp + 1._sp / 5._sp)
        issig = 1._sp * sqrt(gamm(1._sp + 2._sp / 5._sp) - gamm(1._sp + 1._sp / 5._sp)**2._sp )
        sout = weibull_pdf(sat_weibull, c=5._sp, l=0._sp, s=1._sp)
     case default
        sout = 0.0_sp
     end select
     if (i .eq. 15) then
        smean = sum(sat_weibull * sout * dsat_weibull)
        svar  = sum((sat_weibull-smean)**2 * sout * dsat_weibull)
     else if (i .eq. 14) then
        smean = sum(sat_gamma * sout * dsat_gamma)
        svar  = sum((sat_gamma-smean)**2 * sout * dsat_gamma)
     else if (i .eq. 13) then
        smean = sum(sat_beta * sout * dsat_beta)
        svar  = sum((sat_beta-smean)**2 * sout * dsat_beta)
     else
        smean = sum(sat * sout * dsat)
        svar  = sum((sat-smean)**2 * sout * dsat)
     end if
     sstd  = sqrt(svar)
     if (abs(smean-isloc) > 1.0e-2_sp) isgood = .false.
     if (abs(sstd-issig) > 1.0e-2_sp) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions integral SP failed: ', i, isloc, smean, issig, sstd
  enddo


  ! Finish
  allgood = allgood .and. isgood
  if (allgood) then
     write(*,*) 'mo_distributions ', color('o.k.', c_green)
  else
     write(*,*) 'mo_distributions ', color('failed!', c_red)
  endif

END PROGRAM main
