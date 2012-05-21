PROGRAM anneal_test

  use mo_kind,   only: dp, i4, i8
  use mo_anneal, only: anneal
  use mo_cost,   only: cost_dp, range_dp

  IMPLICIT NONE

  REAL(DP), DIMENSION(4)            :: para
  REAL(DP)                          :: temperature
  REAL(DP)                          :: costbest
  REAL(DP), DIMENSION(size(para,1)) :: parabest

  REAL(DP)                          :: costbestAll
  REAL(DP), DIMENSION(size(para,1)) :: parabestAll

  INTEGER(I4)                       :: i, runs
  INTEGER(I8), DIMENSION(3)         :: seeds

  runs = 4_i4

  para(1) = 2.0_dp  
  para(2) = 10.0_dp
  para(3) = 0.5_dp
  para(4) = 0.4_dp

  ! Initialization
  temperature = 10.0_dp
  print*, 'Initial parameter set:       ',para
  print*, 'Initial cost function value: ',cost_dp(para(:))
  print*, '-----------------------------------'
  costbestAll = cost_dp(para(:))
  parabestAll = para(:)

  ! Run Simulated Annealing 10 times
  open(unit=1, file='Anneal.out',   status='unknown')
  write (1, *) 'cost           para(1)        para(2)        para(3)        para(4)        '
  do i=1,runs
     ! Setting the seeds: Only used to generate comparable results for test case
     seeds(1) = int(i,i8)*914_i8
     seeds(2) = seeds(1) + 1000_i8
     seeds(3) = seeds(2) + 1000_i8
     !
     call anneal(cost_dp, para, range_dp, temperature, costbest, parabest, seeds_in=seeds,&
                 nITERmax_in=150000_i4,eps_in=0.00001_dp,printflag_in=.false.)
     if (costbestAll .gt. costbest) then
        costbestAll = costbest
        parabestAll = parabest
     end if
     print*,'Run ',i,':   cost = ',costbest
     write (1, '(F15.7,4(F15.7))') costbest, parabest
  end do
  close(unit=1,  status='keep')

  ! Print best overall solutions
  print*, '-----------------------------------'
  print*, 'Best solution found in ',runs,' Runs'
  print*, '   costbest = ', costbestAll
  print*, '   parabest = ', parabestAll

  ! Is program running properly?   costbestAll = 7.8976644995166545E-02
  print*, '-----------------------------------'
  if ( anint(costbestAll*100000) .eq. 7898._dp ) then
     print*, 'mo_anneal: o.k.'
  else
     print*, 'failed '
  end if

END PROGRAM anneal_test
