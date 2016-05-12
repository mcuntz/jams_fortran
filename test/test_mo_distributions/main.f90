PROGRAM main

  USE mo_kind,          ONLY: i4, dp, sp
  USE mo_utils,         ONLY: ne
  USE mo_distributions, ONLY: ep, ep01, laplace, laplace01, normal, normal01, &
       sep, sep01, sstudentt, sstudentt01, studentt, studentt01

  IMPLICIT NONE

  ! data
  INTEGER(i4), PARAMETER  :: nn = 200
  REAL(dp)                :: out0
  REAL(dp), DIMENSION(nn) :: dat1, out1
  REAL(dp)                :: mu, sig, skew, kurt, nu
  REAL(sp)                :: sout0
  REAL(sp), DIMENSION(nn) :: sat1, sout1
  REAL(sp)                :: smu, ssig, sskew, skurt, snu

  INTEGER(i4) :: i
  LOGICAL     :: isgood, allgood

  Write(*,*) ''
  Write(*,*) 'Test mo_distributions.f90'

  allgood = .true.

  ! DP
  forall(i=1:nn) dat1(i) = real(i,dp)/real(nn,dp)*20.0_dp - 10.0_dp
  mu   = 1.0_dp
  sig  = 2.0_dp
  skew = 2.0_dp
  kurt = 0.5_dp
  nu   = 10.0_dp

  ! integral
  do i=1, 12
     isgood = .true.
     select case(i)
     case(1) ! Exponential Power
        out1 = ep(dat1, mu, sig, kurt)
     case(2)
        out1 = ep01(dat1, kurt)
     case(3)! Laplace
        out1 = laplace(dat1, mu, sig)
     case(4)
        out1 = laplace01(dat1)
     case(5) ! Normal
        out1 = normal(dat1, mu, sig)
     case(6)
        out1 = normal01(dat1)
     case(7) ! Skew Exponential Power
        out1 = sep(dat1, mu, sig, skew, kurt)
     case(8)
        out1 = sep01(dat1, skew, kurt)
     case(9) ! Skewed Student t
        out1 = sstudentt(dat1, nu, mu, sig, skew)
     case(10)
        out1 = sstudentt01(dat1, nu, skew)
     case(11) ! Student t
        out1 = studentt(dat1, nu, mu, sig)
     case(12)
        out1 = studentt01(dat1, nu)
     case default
        continue
     end select
     out0 = sum((dat1(2:nn)-dat1(1:nn-1))*out1(1:nn-1))+(dat1(nn)-dat1(nn-1))*out1(nn)
     if (abs(1.0_dp-out0) > 1.0e-2_dp) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions integral failed: ', i, out0
  enddo

  ! normal vs normal01
  do i=1, 6
     isgood = .true.
     select case(i)
     case(1)
        out1 = ep(dat1, kurt=kurt) - ep01(dat1, kurt)
     case(2)
        out1 = laplace(dat1, 0.0_dp) - laplace01(dat1)
     case(3)
        out1 = normal(dat1, 0.0_dp) - normal01(dat1)
     case(4)
        out1 = sep(dat1, 0.0_dp, 1.0_dp, skew, kurt) - sep01(dat1, skew, kurt)
     case(5)
        out1 = sstudentt(dat1, nu, skew=skew) - sstudentt01(dat1, nu, skew)
     case(6)
        out1 = studentt(dat1, nu, sig=1.0_dp) - studentt01(dat1, nu)
     case default
        continue
     end select
     if (any(ne(out1, 0.0_dp))) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions comparison to 01 versions failed: ', i, out1
  enddo

  ! transformation
  do i=1, 6
     isgood = .true.
     select case(i)
     case(1)
        out1 = ep(dat1, mu, sig, kurt=kurt) - ep((dat1-mu)/sig, kurt=kurt)/sig
     case(2)
        out1 = laplace(dat1, mu, sig) - laplace((dat1-mu)/sig)/sig
     case(3)
        out1 = normal(dat1, mu, sig) - normal((dat1-mu)/sig)/sig
     case(4)
        out1 = sep(dat1, mu, sig, skew=skew, kurt=kurt) - sep((dat1-mu)/sig, skew=skew, kurt=kurt)/sig
     case(5)
        out1 = sstudentt(dat1, nu, mu, sig, skew=skew) - sstudentt((dat1-mu)/sig, nu, skew=skew)/sig
     case(6)
        out1 = studentt(dat1, nu, mu, sig) - studentt((dat1-mu)/sig, nu)/sig
     case default
        continue
     end select
     if (any(ne(out1, 0.0_dp))) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions transformation failed: ', i, out1
  enddo

  ! transformation 0D
  do i=1, 6
     isgood = .true.
     select case(i)
     case(1)
        out0 = ep(0.5_dp, mu, sig, kurt=kurt) - ep((0.5_dp-mu)/sig, kurt=kurt)/sig
     case(2)
        out0 = laplace(0.5_dp, mu, sig) - laplace((0.5_dp-mu)/sig)/sig
     case(3)
        out0 = normal(0.5_dp, mu, sig) - normal((0.5_dp-mu)/sig)/sig
     case(4)
        out0 = sep(0.5_dp, mu, sig, skew=skew, kurt=kurt) - sep((0.5_dp-mu)/sig, skew=skew, kurt=kurt)/sig
     case(5)
        out0 = sstudentt(0.5_dp, nu, mu, sig, skew=skew) - sstudentt((0.5_dp-mu)/sig, nu, skew=skew)/sig
     case(6)
        out0 = studentt(0.5_dp, nu, mu, sig) - studentt((0.5_dp-mu)/sig, nu)/sig
     case default
        continue
     end select
     if (ne(out0, 0.0_dp)) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions transformation failed: ', i, out0
  enddo

  ! ep/sep vs normal/laplace
  do i=1, 9
     isgood = .true.
     select case(i)
     case(1)
        out1 = ep(dat1, 0.0_dp, 1.0_dp, 0.0_dp) - normal(dat1, 0.0_dp, 1.0_dp)
     case(2)
        out1 = ep01(dat1, kurt=0.0_dp) - normal01(dat1)
     case(3)
        out1 = sep(dat1, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp) - normal(dat1, 0.0_dp, 1.0_dp)
     case(4)
        out1 = sep01(dat1) - normal01(dat1)
     case(5)
        out1 = ep(dat1, 0.0_dp, 1.0_dp, 1.0_dp) - laplace(dat1, 0.0_dp, 1.0_dp/sqrt(2.0_dp))
     case(6)
        out1 = ep(dat1, 0.0_dp, sqrt(2.0_dp), 1.0_dp) - laplace(dat1, 0.0_dp, 1.0_dp)
     case(7)
        out1 = ep(dat1, 0.0_dp, sqrt(2.0_dp), 1.0_dp) - laplace01(dat1)
     case(8)
        out1 = sep(dat1, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp) - laplace(dat1, 0.0_dp, 1.0_dp/sqrt(2.0_dp))
     case(9)
        out1 = sep(dat1, 0.0_dp, sqrt(2.0_dp), 1.0_dp, 1.0_dp) - laplace(dat1, 0.0_dp, 1.0_dp)
     case default
        continue
     end select
     if (any(abs(out1) > 1.0e-10)) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions ep/sep comparison to normal/laplace failed: ', i, out1
  enddo

  ! skewed vs unskewed
  do i=1, 2
     isgood = .true.
     select case(i)
     case(1)
        out1 = sep(dat1, mu, sig, 1.0_dp, kurt) - ep(dat1, mu, sig, kurt)
     case(2)
        out1 = sstudentt(dat1, nu, mu, sig, 1.0_dp) - studentt(dat1, nu, mu, sig)
     case default
        continue
     end select
     if (any(abs(out1) > 1.0e-10)) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions skewed vs non-skewed failed: ', i, out1
  enddo


  ! SP
  forall(i=1:nn) sat1(i) = real(i,sp)/real(nn,sp)*20.0_sp - 10.0_sp
  smu   = 1.0_sp
  ssig  = 2.0_sp
  sskew = 2.0_sp
  skurt = 0.5_sp
  snu   = 10.0_sp

  isgood = .true.
  ! integral
  do i=1, 12
     isgood = .true.
     select case(i)
     case(1)
        sout1 = ep(sat1, smu, ssig, skurt)
     case(2)
        sout1 = ep01(sat1, skurt)
     case(3)
        sout1 = laplace(sat1, smu, ssig)
     case(4)
        sout1 = laplace01(sat1)
     case(5)
        sout1 = normal(sat1, smu, ssig)
     case(6)
        sout1 = normal01(sat1)
     case(7)
        sout1 = sep(sat1, smu, ssig, sskew, skurt)
     case(8)
        sout1 = sep01(sat1, sskew, skurt)
     case(9)
        sout1 = sstudentt(sat1, snu, smu, ssig, sskew)
     case(10)
        sout1 = sstudentt01(sat1, snu, sskew)
     case(11)
        sout1 = studentt(sat1, snu, smu, ssig)
     case(12)
        sout1 = studentt01(sat1, snu)
     case default
        continue
     end select
     sout0 = sum((sat1(2:nn)-sat1(1:nn-1))*sout1(1:nn-1))+(sat1(nn)-sat1(nn-1))*sout1(nn)
     if (abs(1.0_sp-sout0) > 1.0e-2_sp) isgood = .false.
     allgood = allgood .and. isgood
     if (.not. isgood) write(*,*) 'mo_distributions integral-sp failed: ', i, sout0
  enddo


  ! Finish
  allgood = allgood .and. isgood
  if (allgood) then
     write(*,*) 'mo_distributions o.k.'
  else
     write(*,*) 'mo_distributions failed!'
  endif

END PROGRAM main
