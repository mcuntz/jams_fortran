PROGRAM anneal_test

  use mo_kind,    only: dp, i4, i8
  use mo_anneal,  only: anneal_valid, GetTemperature_valid
  use mo_cost,    only: cost_dp, range_dp
  use mo_xor4096, only: get_timeseed

  IMPLICIT NONE

  REAL(DP), DIMENSION(4)            :: para
  REAL(DP)                          :: temperature
  REAL(DP)                          :: costbest
  REAL(DP), DIMENSION(size(para,1)) :: parabest

  REAL(DP)                          :: costbestAll
  REAL(DP), DIMENSION(size(para,1)) :: parabestAll

  INTEGER(I4)                       :: i, runs
  INTEGER(I8), DIMENSION(3)         :: seeds
  REAL(DP)                          :: Tstart, Tend

  ! time dependent seeds
  call get_timeseed(seeds)
  print*, 'time dependent seeds would be: ', seeds

  runs = 4_i4

  para(1) = 2.0_dp  
  para(2) = 10.0_dp
  para(3) = 0.5_dp
  para(4) = 0.4_dp

  ! Initialization
  print*, '-----------------------------------'
  print*, '   INITIALIZATION                  '
  print*, '-----------------------------------'
  print*, 'Initial parameter set:       ',para
  print*, 'Initial cost function value: ',cost_dp(para(:))
  costbestAll = cost_dp(para(:))
  parabestAll = para(:)
  print*, '-----------------------------------'
  print*, '   INITIAL TEMPERATURE             '
  print*, '-----------------------------------'
  print*, 'Estimation of Initial Temperature: '
  !temperature = 10.0_dp
  seeds(1) = 854_i8
  seeds(2) = seeds(1) + 1000_i8
  print*, 'seeds used:                        ', seeds(1:2)
  temperature = GetTemperature_valid( para, cost_dp, range_dp, 0.95_dp, samplesize_in=500_i4, &
       seeds_in=seeds(1:2), printflag_in=.true.)

  print*, '-----------------------------------'
  print*, '   SIMULATED ANNEALING             '
  print*, '-----------------------------------'

  ! Run Simulated Annealing <runs> times
  open(unit=1, file='Anneal.out',   status='unknown')
  write (1, *) 'cost           para(1)        para(2)        para(3)        para(4)        '
  do i=1,runs
     ! Setting the seeds: Only used to generate comparable results for test case
     seeds(1) = int(i,i8)*259_i8
     seeds(2) = seeds(1) + 1000_i8
     seeds(3) = seeds(2) + 1000_i8
     print*, 'seeds used: ', seeds(1:3)
     !
     call cpu_time(Tstart)
     parabest = anneal_valid(cost_dp, para, range_dp, maxit_in=.false., &
          temp_in=temperature, seeds_in=seeds,&
          LEN_in=250_i4,nITERmax_in=150000_i4,eps_in=0.00001_dp,&
          printflag_in=.false., &
          funcbest=costbest)
     call cpu_time(Tend)
     if (costbestAll .gt. costbest) then
        costbestAll = costbest
        parabestAll = parabest
     end if
     print*,'Run ',i,':   cost = ',costbest,'  (CPU time = ',Tend-Tstart,' sec)'
     write (1, '(F15.7,4(F15.7))') costbest, parabest
  end do
  close(unit=1,  status='keep')

  ! Print best overall solutions
  print*, '-----------------------------------'
  print*, 'Best solution found in ',runs,' Runs'
  print*, '   costbest = ', costbestAll
  print*, '   parabest = ', parabestAll

  ! Is program running properly?   costbestAll = 8.4057185000218151E-02
  print*, '-----------------------------------'
  if ( anint(costbestAll*100000) .eq. 8406._dp ) then
     print*, 'mo_anneal: o.k.'
  else
     print*, 'mo_anneal: failed '
  end if

END PROGRAM anneal_test
