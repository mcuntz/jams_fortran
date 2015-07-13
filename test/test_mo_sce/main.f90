program test_sce

  use mo_kind,             only: i4, i8, dp
  use mo_sce,              only: sce
  use mo_opt_functions,    only: ackley

  implicit none

  real(dp), dimension(:),   allocatable :: pini
  real(dp), dimension(:,:), allocatable :: prange
  real(dp), dimension(:),   allocatable :: opt
  logical,  dimension(:),   allocatable :: mask

  real(dp)                              :: bestf      ! best function value = history(-1) = ackley(opt)
  integer(i8)                           :: neval      ! number of function evaluations = size(history,1)
  real(dp), dimension(:),   allocatable :: history    ! best function value found after ith iteration

  integer(i4)   :: n, iPar
  logical       :: isgood
  logical       :: parallel

  ! Dimension of test function = number of parameters
  n = 30

  ! Allocation
  allocate(pini(n))
  allocate(prange(n,2))
  allocate(opt(n))
  allocate(mask(n))

  ! Initialization
  pini        =  0.5_dp
  opt         =  pini
  do iPar = 1, n
     prange(iPar,1) = -10.0_dp
     prange(iPar,2) =  10.0_dp
  end do
  mask(:) = .true.
  pini(:) = 0.5_dp
  parallel=.false.

  opt  = sce(ackley,       & ! Mandatory IN: Objective function
       pini,               & ! Mandatory IN: initial guess
       prange,             & ! Mandatory IN: range for each parameter (min, max)
       mymaxn=30000_i8,    & ! Optional  IN: maximal number of function evaluations
       mymaxit=.false.,    & ! Optional  IN: maximize (true) or minimize (false) objective
       mykstop=10_i4,      & ! Optional  IN: number of shuffling loops in which the criterion value must
       !                                     change by given percentage before optimiz. is terminated
       mypcento=0.0001_dp, & ! Optional  IN: percentage by which the criterion value must change in
       !                                     <kstop> number of shuffling loops
       mypeps=0.001_dp,    & ! Optional  IN: optimization is terminated if volume of complex has
       !                                     converged to given percentage of feasible space
       myseed=10987_i8,    & ! Optional  IN: seed of random number generator
       myngs=2_i4,         & ! Optional  IN: number of complexes in the initial population
       mynpg=2*n+1,        & ! Optional  IN: number of points in each complex
       mynps=n+1,          & ! Optional  IN: number of points in a sub-complex
       mynspl=2*n+1,       & ! Optional  IN: number of evolution steps allowed for each complex before
       !                                     complex shuffling
       mymings=2_i4,       & ! Optional  IN: minimum number of complexes required, if the number of
       !                                     complexes is allowed to reduce as the
       !                                     optimization proceeds
       myiniflg=1_i4,      & ! Optional  IN: flag on whether to include the initial point
       !                                     in population (1) or not (0)
       myprint=0_i4,       & ! Optional  IN: flag for controlling print-out
       !                                        0, print information on the best point of the population
       !                                        1, print information on every point of the population
       !                                        2, no printing
       mymask=mask,        & ! OPTIONAL  IN:
       myalpha=0.8_dp,     & ! Optional  IN: parameter for reflection  of points in complex
       mybeta=0.45_dp,     & ! Optional  IN: parameter for contraction of points in complex
       tmp_file='best_objective_make_check_test_file', & ! Optional IN: writes best objective after each evolution loop
       popul_file='population_make_check_test_file',   & ! Optional IN: writes whole population and respective objectives
       !                                !                                  for each evolution loop
       restart_file='restart_make_check_test_file',    & ! Optional IN: restart file name
       parallel=parallel,             & ! Optional IN: if .false. sce will not run in parallel (although compiled with OpenMP)
       !                                               if .true.  sce will run in parallel when compiled with OpenMP and
       !                                                          in single thread when not compiled with OpenMP
       bestf=bestf,        & ! Optional OUT: objective function value of best parameter set
       neval=neval,        & ! Optional OUT: number of function evaluations needed
       history=history     & ! Optional OUT: history(k) = best objective after k function evaluations
       )

  write(*,*) ''
  write(*,*) 'number of function evaluations: neval = ', neval
  write(*,*) 'best function value found:      bestf = ', ackley(opt)
  write(*,*) 'global minimal function value:  optf  = ', 0.0_dp
  write(*,*) ''

  if (parallel) then
     ! compiled without OpenMP
     isgood = .true.
     isgood = isgood .and. (neval .eq. 4455_i4)
     isgood = isgood .and. (nint(bestf*10000000_dp) .eq. 104439)

     ! compiled with OpenMP
     !$ isgood = .true.
     !$ write(*,*) 'mo_sce: It is not possible to check if sce runs properly with OpenMP.'

  else
     ! no matter if compiler with or without OpenMP
     isgood = .true.
     isgood = isgood .and. (neval .eq. 4455_i4)
     isgood = isgood .and. (nint(bestf*10000000_dp) .eq. 104439)
  end if

  if (isgood) then
     write(*,*) 'mo_sce o.k.'
  else
     write(*,*) 'mo_sce failed!'
  endif


  ! Check restart
  bestf   = 0.
  neval   = 0
  history = 0.

  opt = sce(ackley,       & ! Mandatory IN  : Objective function
       pini,               & ! Mandatory IN  : initial guess
       prange,             & ! Mandatory IN  : range for each parameter (min, max)
       restart=.true.,     & ! Do the restart
       restart_file='restart_make_check_test_file', & ! Optional IN: restart file name
       bestf=bestf,        & ! Optional  OUT : objective function value of best parameter set
       neval=neval,        & ! Optional  OUT : number of function evaluations needed
       history=history     & ! Optional  OUT : history(k) = best objective after k function evaluations
       )

  write(*,*) ''
  write(*,*) 'number of function evaluations: neval = ', neval
  write(*,*) 'best function value found:      bestf = ', ackley(opt)
  write(*,*) 'global minimal function value:  optf  = ', 0.0_dp
  write(*,*) ''

  if (parallel) then
     ! compiled without OpenMP
     isgood = .true.
     isgood = isgood .and. (neval .eq. 4455_i4)
     isgood = isgood .and. (nint(bestf*10000000_dp) .eq. 104439)
     ! compiled with OpenMP
     !$ isgood = .true.
     !$ write(*,*) 'mo_sce: It is not possible to check if sce runs properly with OpenMP.'
  else
     ! no matter if compiler with or without OpenMP
     isgood = .true.
     isgood = isgood .and. (neval .eq. 4455_i4)
     isgood = isgood .and. (nint(bestf*10000000_dp) .eq. 104439)
  end if

  if (isgood) then
     write(*,*) 'mo_sce restart o.k.'
  else
     write(*,*) 'mo_sce restart failed!'
  endif

end program test_sce
