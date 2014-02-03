PROGRAM anneal_test

  use mo_kind,    only: dp, i4, i8
  use mo_anneal,  only: anneal            !, anneal_valid
  use mo_anneal,  only: GetTemperature    !, GetTemperature_valid
  use mo_cost,    only: cost_dp, range_dp !, cost_valid_dp
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
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: history

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
  seeds(1) = 854_i8
  seeds(2) = seeds(1) + 1000_i8
  print*, 'seeds used:                        ', seeds(1:2)
  temperature = GetTemperature( para, cost_dp, 0.95_dp, &
       ! optionals
       prange_func=range_dp, &
       samplesize=500_i4, &
       seeds=seeds(1:2), &
       printflag=.true.)

  print*, '-----------------------------------'
  print*, '   SIMULATED ANNEALING             '
  print*, '-----------------------------------'

  ! Run Simulated Annealing <runs> times
  do i=1,runs
     ! Setting the seeds: Only used to generate comparable results for test case
     seeds(1) = int(i,i8)*259_i8
     seeds(2) = seeds(1) + 1000_i8
     seeds(3) = seeds(2) + 1000_i8
     print*, 'seeds used: ', seeds(1:3)
     !
     call cpu_time(Tstart)
     parabest = anneal(cost_dp, para, &
          ! optionals
          prange_func=range_dp, &
          maxit=.false., &
          temp=temperature, &
          seeds=seeds,&
          LEN=250_i4, &
          nST=5_i4, &
          nITERmax=150000_i4, &
          eps=0.00001_dp,&
          acc=0.1_dp,&
          DT=0.90_dp, &
          printflag=.false., &
          funcbest=costbest, history=history, &
          weight=(/ 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp  /), &
          maskpara=(/ .true., .true., .true., .true. /), &
          undef_funcval=9999.0_dp, &
          ! settings for standard simulated annealing
          reflectionFlag=.False., &
          pertubFlexFlag=.True., &
          changeParaMode=1_i4)
          ! settings for dds generation of new parameter sets
          ! reflectionFlag=.True., &
          ! pertubFlexFlag=.False., &
          ! changeParaMode=3_i4)
     call cpu_time(Tend)
     if (costbestAll .gt. costbest) then
        costbestAll = costbest
        parabestAll = parabest
     end if
     print*,'Run ',i,':   cost = ',costbest,'  (CPU time = ',Tend-Tstart,' sec)'
  end do

  ! Print best overall solutions
  print*, '-----------------------------------'
  print*, 'Best solution found in ',runs,' Runs'
  print*, '   costbest = ', costbestAll
  print*, '   parabest = ', parabestAll

  ! Is program running properly?   costbestAll = 3.1142480812726376E-02
  print*, '-----------------------------------'
  if ( nint(costbestAll*100000) .eq. 3114 ) then
     print*, 'mo_anneal: o.k.'
  else
     print*, 'mo_anneal: failed '
  end if

END PROGRAM anneal_test
