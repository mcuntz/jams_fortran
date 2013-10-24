program test_mcmc

  use mo_kind,       only: dp,i4, i8
  use mo_likelihood, only: setmeas, loglikelihood_dp
  use mo_mcmc,       only: mcmc
  use mo_moment,     only: mean, stddev

  ! for running mcmc
  real(dp)                              :: p
  real(dp)                              :: likeli_new
  real(dp)                              :: stddev_new
  real(dp), dimension(3)                :: parabest
  real(dp), dimension(3,2)              :: rangePar
  logical,  dimension(3)                :: maskpara
  real(dp), dimension(:,:), allocatable :: burnin_paras
  real(dp), dimension(:,:), allocatable :: mcmc_paras
  integer(i8)                           :: seed
  ! for result handling
  integer(i4)                           :: i, samples

  ! loading the data points observed
  call setmeas()

  ! best parameter set found for dataset using e.g. Simlated Annealing or DDS
  parabest = (/ 0.9999815795293934_dp, 2.0035191695594499_dp, 2.8675479872139675_dp /)

  ! your best parameter set should have the largest likelihood
  ! sigma of errors is unknown --> initial guess e.g. 1.0_dp
  !    stddev_new = standard deviation of errors using paraset
  !    likeli_new = likelihood using stddev_new
  p = loglikelihood_dp(parabest,1.0_dp,likeli_new=likeli_new,stddev_new=stddev_new)
  write(*,*) 'guessed log-likelihood = ',p
  write(*,*) 'best log-likelihood    = ',likeli_new

  ! initializing the ranges of parameters
  rangePar(:,1) = -10.0_dp
  rangePar(:,2) =  10.0_dp

  ! which parameter will be sampled
  maskpara = (/ .true., .true., .true. /)

  ! fixing seeds for test case:
  seed = 62001519_i8

  ! starting MCMC: 
  !     (1) Burn-in will be performed to optimize settings for MCMC
  !     (2) posterior distribution of the parameters at the minimum (best parameterset) 
  !         will be sampled by MCMC
  call mcmc(loglikelihood_dp, parabest, rangePar, mcmc_paras, burnin_paras, &
       ParaSelectMode_in=2_i4,tmp_file='make_check_test_file',              &
       seed_in=seed, loglike_in=.true., maskpara_in=maskpara, printflag_in=.true.)

  write(*,*)''
  write(*,*) '--------------------------------------------------'
  write(*,*) 'Results '
  write(*,*) '--------------------------------------------------'
  write(*,*) ''
  write(*,*) 'Number of parameter sets sampled'
  ! number of samples per chain x number of chains
  samples = size(mcmc_paras,1) 
  write(*,*) samples

  write(*,*) ''
  write(*,*) 'Mean parameter value +/- standard deviation'
  do i=1,size(parabest)
     if (maskpara(i)) then
        write(*,*) 'para #',i, ' = ',mean(mcmc_paras(:,i)), '+/-',stddev(mcmc_paras(:,i))
     end if
  end do

  ! Does it run properly?
  ! para # 1  =    0.9999781778559107 +/-   6.3209436919767900E-05
  ! para # 2  =    2.0038078713380076 +/-   6.6280930370445947E-03
  ! para # 3  =    2.8654065450088502 +/-   0.1498533843248984
  write(*,*) '-----------------------------------'
  if ( (nint(stddev(mcmc_paras(:,1))*10000000,i4) .eq. 632_i4)   .and. &
       (nint(stddev(mcmc_paras(:,2))*10000000,i4) .eq. 66281_i4) .and. &
       (nint(stddev(mcmc_paras(:,3))*10000000,i4) .eq. 1498534_i4) ) then
     write(*,*) 'mo_mcmc: o.k.'
  else
     write(*,*) 'mo_mcmc: failed '
  end if

end program test_mcmc
