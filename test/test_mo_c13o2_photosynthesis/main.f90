program main

  use mo_kind,                 only: dp, i4
  use mo_utils,                only: ne
  use mo_isotope,              only: vpdbc13
  use mo_c13o2_photosynthesis, only: Vstarch, Rsucrose, Rphoto, Rstarch
  use mo_c13o2_photosynthesis, only: c13o2_discrimination_simple, init_ratio_leaf_pools!, c13o2_discrimination

  implicit none

  integer(i4), parameter :: nn = 10
  logical, dimension(nn) :: isc3
  real(dp), dimension(nn) :: gs, Ca, Ci, A, Rd, GPP, Tl, Ra, disc, A13, ntmp
  real(dp) :: Rinitc3, Rinitc4

  logical :: isgood, allgood

  write(*,*) ''
  write(*,*) 'Test mo_c13o2'
  
  allgood = .true.
  
  ! -----------------------
  ! Init

  isgood = .true.

  isc3(1:nn:2) = .true.
  isc3(2:nn:2) = .false.
  Rinitc3 = (1._dp-29.e-3_dp)*vpdbc13
  Rinitc4 = (1._dp-12e-3_dp)*vpdbc13
  call init_ratio_leaf_pools(Vstarch, Rsucrose, Rphoto, Rstarch, Rinitc3, Rinitc4, isc3)

  if (any(ne(Vstarch,0._dp))) isgood = .false.
  if (any(ne(Rsucrose(1::2),Rinitc3))) isgood = .false.
  if (any(ne(Rphoto(2::2),Rinitc4))) isgood = .false.
  if (any(ne(Rstarch,merge(Rinitc3,Rinitc4,isc3)))) isgood = .false.

  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_c13o2 init o.k.'
  else
     write(*,*) 'mo_c13o2 init failed!'
  endif
  
  ! -----------------------
  ! Simple discrimination

  isgood = .true.

  Ca = 353._dp  
  Rd = 1._dp
  Ra = vpdbc13
  gs = 0.2_dp

  ! C3
  isc3 = .true.
  Ci   = 0.7_dp*Ca
  A    = gs*(Ca-Ci)
  GPP  = A + Rd
  Tl   = 290._dp
  call c13o2_discrimination_simple(isc3, GPP, Rd, Ca, Ci, Tl, Ra, disc, A13)

  ntmp = 4.4_dp + (27.0_dp-4.4_dp)*0.7_dp
  if (any(abs(disc*1000._dp-ntmp) > epsilon(1._dp))) isgood = .false.
  ntmp = 1._dp - A13/A/Ra
  if (any(abs(disc-ntmp) > epsilon(1._dp))) isgood = .false.

  ! C4
  ! from Cousins et al (2006), Plant Physiology, 141(1), 232-242 - Fig. 3 high kCa -> Disc ~ 3 permil
  isc3 = .false.
  Ci   = 0.4_dp*Ca
  A    = gs*(Ca-Ci)
  GPP  = A + Rd
  Tl   = 303.15_dp
  call c13o2_discrimination_simple(isc3, GPP, Rd, Ca, Ci, Tl, Ra, disc, A13) ! Disc = 2.9342514824082233

  if (any(abs(disc*1000._dp-2.9342514824082233_dp) > epsilon(1._dp))) isgood = .false.

  ! Night - C3 and C4
  isc3(::2)  = .true.
  isc3(2::2) = .false.
  GPP  = 0._dp
  A    = -Rd
  Ci   = Ca - A/gs
  Tl   = 290._dp
  call c13o2_discrimination_simple(isc3, GPP, Rd, Ca, Ci, Tl, Ra, disc, A13)

  ntmp = 4.4_dp
  if (any(abs(disc*1000._dp-ntmp) > epsilon(1._dp))) isgood = .false.
  ntmp = 1._dp - A13/A/Ra
  if (any(abs(disc-ntmp) > epsilon(1._dp))) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model delta o.k.'
  else
     write(*,*) 'mo_isotope_pool_model delta failed!'
  endif

  
  ! -----------------------
  ! Discrimination

  ! ToDo
  
  
  ! ----------------------
  
  if (allgood) then
     write(*,*) 'mo_c13o2 o.k.'
  else
     write(*,*) 'mo_c13o2 failed!'
  endif

end program main
