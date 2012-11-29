program test_mcmc

  use mo_kind,       only: dp,i4, i8
  use mo_likelihood, only: setmeas, loglikelihood_dp, stdev_data_dp
  use mo_mcmc,       only: mcmc
  use mo_moment,     only: mean, stddev

  ! for running mcmc
  real(dp)                              :: p
  real(dp)                              :: sigma           ! standarddeviation at each data point
  real(dp), dimension(3)                :: parabest
  real(dp), dimension(3,2)              :: rangePar
  logical,  dimension(3)                :: maskpara
  real(dp), dimension(:,:), allocatable :: burnin_paras
  real(dp), dimension(:,:), allocatable :: mcmc_paras
  integer(i8), dimension(5,3)           :: seeds
  ! for result handling
  integer(i4)                           :: i, samples

  ! loading the data points observed
  call setmeas()

  ! best parameter set found for dataset using e.g. Simlated Annealing or DDS
  parabest = (/ 0.9999815795293934_dp, 2.0035191695594499_dp, 2.8675479872139675_dp /)
  
  ! your best parameter set should have the largest likelihood
  sigma = stdev_data_dp(parabest)
  p = loglikelihood_dp(parabest,sigma)
  print*, 'best log-likelihood = ',p

  ! initializing the ranges of parameters
  rangePar(:,1) = -10.0_dp
  rangePar(:,2) =  10.0_dp

  ! which parameter will be sampled
  maskpara = (/ .true., .true., .true. /)

  ! fixing seeds for test case:
  seeds(1,:) = (/ 62001519_i8, 62002519_i8, 62003519_i8 /)
  seeds(2,:) = (/ 62004519_i8, 62005519_i8, 62006519_i8 /)
  seeds(3,:) = (/ 62007519_i8, 62008519_i8, 62009519_i8 /)
  seeds(4,:) = (/ 62010519_i8, 62011519_i8, 62012519_i8 /)
  seeds(5,:) = (/ 62013519_i8, 62014519_i8, 62015519_i8 /)

  ! starting MCMC: 
  !     (1) Burn-in will be performed to optimize settings for MCMC
  !     (2) posterior distribution of the parameters at the minimum (best parameterset) 
  !         will be sampled by MCMC
  call mcmc(loglikelihood_dp, stdev_data_dp, parabest, rangePar, mcmc_paras, burnin_paras, &
            ParaSelectMode_in=2_i4,tmp_file='tmp_parasets.nc', &
            seeds_in=seeds, loglike_in=.true., maskpara_in=maskpara, printflag_in=.true.)

  print*,''
  print*, '--------------------------------------------------'
  print*, 'Results '
  print*, '--------------------------------------------------'
  print*, ''
  print*, 'Number of parameter sets sampled'
  ! number of samples per chain x number of chains
  samples = size(mcmc_paras,1) 
  print*, samples

  print*, ''
  print*, 'Mean parameter value +/- standard deviation'
  do i=1,size(parabest)
     if (maskpara(i)) then
        print*, 'para #',i, ' = ',mean(mcmc_paras(:,i)), &
                            '+/-',stddev(mcmc_paras(:,i))
     end if
  end do

  ! Is program running properly?  
  !      stddev(para_1) = 6.6512652152650692E-05
  !      stddev(para_2) = 6.9419315068863776E-03
  !      stddev(para_3) = 0.1550643250843179
  !
  print*, '-----------------------------------'
  if ( (nint(stddev(mcmc_paras(:,1))*10000000,i4) .eq. 632_i4)   .and. &
       (nint(stddev(mcmc_paras(:,2))*10000000,i4) .eq. 66281_i4) .and. &
       (nint(stddev(mcmc_paras(:,3))*10000000,i4) .eq. 1498534_i4) ) then
     print*, 'mo_mcmc: o.k.'
  else
     print*, 'mo_mcmc: failed '
  end if

end program test_mcmc
