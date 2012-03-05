PROGRAM anneal_test

  use mo_kind,   only: sp, dp,i4
  use mo_anneal, only: anneal
  use mo_cost,   only: cost_dp, cost_sp

  IMPLICIT NONE

  REAL(DP), DIMENSION(4, 3)         :: para
  REAL(DP)                          :: temperature
  REAL(DP)                          :: costbest
  REAL(DP), DIMENSION(size(para,1)) :: parabest

  REAL(DP)                          :: costbestAll
  REAL(DP), DIMENSION(size(para,1)) :: parabestAll

  INTEGER(I4)                       :: i, runs

  runs = 1_i4

  para(:,1) = 1.0_dp
  para(:,2) = 0.0_dp    ! min-value of parameters
  para(:,3) = 30.0_dp    ! max-value of parameters

  ! Initialization
  temperature = 10.0_dp
  print*, cost_sp(real(para(:,1),sp))
  pause
  costbestAll = cost_dp(para(:,1))
  parabestAll = para(:,1)

  ! Run Simulated Annealing 10 times
  open(unit=1, file='Anneal.out',   status='unknown')
  write (1, *) 'cost           para(1)        para(2)        para(3)        para(4)        '
  do i=1,runs
     call anneal(cost_dp, para, temperature, costbest, parabest, &
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

END PROGRAM anneal_test
