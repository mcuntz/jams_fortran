program test_sce

  use mo_kind,             only: i4, i8, dp
  use mo_sce,              only: sce
  use mo_opt_functions,    only: ackley

  implicit none

  real(dp), dimension(:),   allocatable :: pini
  real(dp), dimension(:,:), allocatable :: prange
  real(dp), dimension(:),   allocatable :: opt

  real(dp)                              :: bestf      ! best function value = history(-1) = ackley(opt)
  integer(i4)                           :: neval      ! number of function evaluations = size(history,1)
  real(dp), dimension(:),   allocatable :: history    ! best function value found after ith iteration

  integer(i4)   :: n, iPar
  logical       :: isgood

  isgood = .true.

  ! Dimension of test function = number of parameters
  n = 30

  ! Allocation
  allocate(pini(n))
  allocate(prange(n,2))
  allocate(opt(n))
     
  ! Initialization
  pini        =  0.5_dp
  opt         =  pini
  do iPar = 1, n
     prange(iPar,1) = -10.0_dp
     prange(iPar,2) =  10.0_dp
  end do

  opt  = sce(ackley, pini, prange,                                                                  & ! Mandatory IN 
       mymaxn=30000_i8, mymaxit=.false., mykstop=10_i4, mypcento=0.0001_dp, myseed=10987_i8,        & ! Optional  IN
       myngs=2_i4, mynpg=2*n+1, mynps=n+1, mynspl=2*n+1, mymings=2_i4, myiniflg=1_i4, myprint=0_i4, & ! Optional  IN
       myalpha=0.8_dp, mybeta=0.45_dp,                                                              & ! Optional  IN
       bestf=bestf, neval=neval, history=history                                                    & ! Optional  OUT
       )

  write(*,*) 'number of function evaluations: neval = ', neval
  write(*,*) 'best function value found:      bestf = ', ackley(opt)
  write(*,*) 'global minimal function value:  optf  = ', 0.0_dp
  write(*,*) ''

  isgood = isgood .and. (neval .eq. 4455_i4)
  isgood = isgood .and. (nint(1.0443851441259699E-02*10000000_dp) .eq. 104439)

  if (isgood) then
     write(*,*) 'mo_sce o.k.'
  else
     write(*,*) 'mo_sce failed!'
  endif

end program test_sce
