program main
  !
  use mo_kind,   only: i4, dp, sp
  use mo_ncread, only: Get_NcVar
  use mo_MVN,    only: MVN
  !
  implicit none
  !
  integer(i4)                 :: i
  real(dp), dimension(36,36)  :: B
  real(sp), dimension(36,36)  :: B_sp
  real(dp), dimension(36,100) :: rand
  real(sp), dimension(36,100) :: rand_sp
  real(dp), dimension(36,100) :: res
  real(sp), dimension(36,100) :: res_sp
  real(dp), parameter         :: epsdp = epsilon(1._dp)
  real(sp), parameter         :: epssp = epsilon(1._sp)
  logical                     :: isgood
  !
  isgood = .true.
  !
  ! read lower cholesky factor
  call Get_NcVar( '../FORTRAN_chs_lib/test/test_mo_MVN/CholeskyFactor.nc', 'B', B)
  B_sp = real( B, sp )
  !
  ! do not initialize with time seed to compare results
  do i = 1, size(rand,2)
     rand(:,i) = MVN(B)
  end do
  do i = 1, size(rand_sp,2)
     rand_sp(:,i) = MVN(B_sp)
  end do
  !
  ! read result of above loop from txt file
  open( unit = 250, file = '../FORTRAN_chs_lib/test/test_mo_MVN/Result.txt', action= 'read')
  do i = 1, size(rand,2)
     read(250, *) res(:,i)
  end do
  close(250)
  open( unit = 250, file = '../FORTRAN_chs_lib/test/test_mo_MVN/Result_sp.txt', action= 'read')
  do i = 1, size(rand_sp,2)
     read(250, *) res_sp(:,i)
  end do
  close(250)
  !
  isgood = isgood .and. ( abs(sum(rand) - sum(res)) < epsdp )
  isgood = isgood .and. ( abs(sum(rand_sp)- sum(res_sp))< epssp )
  !
  ! initialize random number generator with timeseed
  rand(:,1) = MVN(B,.true.)
  !
  ! result has to be different then the first time
  do i = 1, size(rand,2)
     rand(:,i) = MVN(B)
  end do
  do i = 1, size(rand_sp,2)
     rand_sp(:,i) = MVN(B_sp)
  end do
  !
  ! result has to be different from before
  isgood = isgood .and. maxval(abs(rand - res)) > epsdp
  isgood = isgood .and. maxval(abs(rand_sp - res_sp)) > epssp
  !
  if ( isgood ) then
     print *, 'mo_MVN is o.k.'
  else
     print *, 'mo_MVN failed!'
  end if
  !
end program main
