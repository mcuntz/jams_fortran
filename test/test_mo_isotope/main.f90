program main

  use mo_kind,    only: dp, i4
  use mo_utils,   only: ne
  use mo_isotope, only: delta, delta1000, isoratio
  use mo_isotope, only: ratio_diff_air2vap, ratio_boundary_air2vap, &
       ratio_diff_vap2air, ratio_boundary_vap2air, &
       vpdbc13, Mc13, Mc13o2

  implicit none

  integer(i4), parameter :: nn = 10
  real(dp), dimension(nn) :: C, Ci, r, d, d1000, eps

  integer(i4) :: i

  logical :: isgood, allgood

  write(*,*) ''
  write(*,*) 'Test mo_isotope'
  
  allgood = .true.
  
  ! -----------------------
  ! Constants

  isgood = .true.

  if (ne(ratio_diff_air2vap,1._dp/ratio_diff_vap2air)) isgood = .false.
  if (ne(ratio_boundary_air2vap,1._dp/ratio_boundary_vap2air)) isgood = .false.
  if (ne(vpdbc13,0.0112372_dp)) isgood = .false.
  if (ne(Mc13*1.e3_dp,13._dp)) isgood = .false.
  if (ne(Mc13o2*1.e3_dp,45._dp)) isgood = .false.

  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope constants o.k.'
  else
     write(*,*) 'mo_isotope constants failed!'
  endif
  
  ! -----------------------
  ! Routines

  isgood = .true.

  C   = 1._dp
  eps = (/ (real(i,dp),i=1,10) /) / 10._dp
  Ci  = C * vpdbc13 * (1._dp+eps)

  ! scalar
  do i=1, nn
     r(i) = isoratio(Ci(i), C(i))
     d(i) = delta(Ci(i), C(i), vpdbc13)
     d1000(i) = delta1000(Ci(i), C(i), vpdbc13)
  end do
  if (any(abs(r/vpdbc13-(1._dp+eps)) > epsilon(1._dp))) isgood = .false.
  if (any(abs(d-eps) > epsilon(1._dp))) isgood = .false.
  if (any(abs(d1000/1.e3_dp-eps) > epsilon(1._dp))) isgood = .false.

  ! vector
  r = isoratio(Ci, C)
  d = delta(Ci, C, vpdbc13)
  d1000 = delta1000(Ci, C, vpdbc13)
  if (any(abs(r/vpdbc13-(1._dp+eps)) > epsilon(1._dp))) isgood = .false.
  if (any(abs(d-eps) > epsilon(1._dp))) isgood = .false.
  if (any(abs(d1000/1.e3_dp-eps) > epsilon(1._dp))) isgood = .false.

  ! default
  C(1)   = 0._dp
  eps(1) = 1._dp
  r = isoratio(Ci, C, 2._dp*vpdbc13)
  d = delta(Ci, C, vpdbc13, 1._dp)
  d1000 = delta1000(Ci, C, vpdbc13, 1000._dp)
  if (any(abs(r/vpdbc13-(1._dp+eps)) > epsilon(1._dp))) isgood = .false.
  if (any(abs(d-eps) > epsilon(1._dp))) isgood = .false.
  if (any(abs(d1000/1.e3_dp-eps) > epsilon(1._dp))) isgood = .false.

  ! precision
  C(2)   = 0.01_dp
  C(3)   = 0.02_dp
  Ci     = C * vpdbc13 * (1._dp+eps)
  eps(2) = 1._dp
  r = isoratio(Ci, C, 2._dp*vpdbc13, 0.01_dp)
  d = delta(Ci, C, vpdbc13, 1._dp, 0.01_dp)
  d1000 = delta1000(Ci, C, vpdbc13, 1000._dp, 0.01_dp)
  if (any(abs(r/vpdbc13-(1._dp+eps)) > epsilon(1._dp))) isgood = .false.
  if (any(abs(d-eps) > epsilon(1._dp))) isgood = .false.
  if (any(abs(d1000/1.e3_dp-eps) > epsilon(1._dp))) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model delta o.k.'
  else
     write(*,*) 'mo_isotope_pool_model delta failed!'
  endif

  ! ----------------------
  
  if (allgood) then
     write(*,*) 'mo_isotope o.k.'
  else
     write(*,*) 'mo_isotope failed!'
  endif

END PROGRAM main
