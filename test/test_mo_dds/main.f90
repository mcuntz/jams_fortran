program main

  use mo_kind, only: i4, i8, dp
  use mo_dds,  only: dds
  use mo_opt_functions, only: griewank

  implicit none

  integer(i4), parameter :: np = 10

  ! Input to DDS
  real(dp),    dimension(np)   :: dv_ini   ! inital value of decision variables 
  real(dp),    dimension(np,2) :: dv_range ! Min/max values of decision variables 

  ! Optional Input to DDS
  real(dp)                     :: r_val    ! DDS perturbation parameter (-> 0.2 by default)
  integer(i8) :: seed                      ! User seed to initialise the random number generator
  integer(i8) :: nIterMax                  ! Maximum number of iteration or function evaluation
  logical     :: to_max                    ! Maximization or minimization of function

  ! Output of DDS
  real(dp), dimension(np) :: dv_opt        ! Best value of decision variables

  ! Test
  logical :: isgood

  write(*,*) '--------------------------------------------------'
  write(*,*) '  Dynamically Dimensioned Search (DDS) Algorithm  '
  write(*,*) '         version 1.1 - by Bryan Tolson            '
  write(*,*) '--------------------------------------------------'

  dv_ini        = (/ -.226265E+01, -.130187E+01, -.151219E+01, 0.133983E+00, 0.988159E+00, &
                     -.495074E+01, -.126574E+02, 0.572684E+00, 0.303864E+01, 0.343031E+01 /)
  dv_range(:,1) = (/ -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0 /)
  dv_range(:,2) = (/ 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0 /)
  r_val    = 0.20_dp
  nIterMax = 100000_i8
  to_max   = .False.
  seed     = 123456789_i8
  dv_opt   = DDS(griewank, dv_ini, dv_range, r=r_val, seed=seed, maxiter=nIterMax, maxit=to_max)

  write(*,*) ''
  write(*,'(A,10F9.4)') 'Should be    3.0825   0.0629   0.0362  -0.0063   0.0026  -7.6358   0.0718  -0.0528  -0.1008   0.0037'
  write(*,'(A,10F9.4)') 'Output is ', dv_opt

  isgood = .true.
  isgood = isgood .and. (anint(10000._dp*dv_opt(1))  ==  30825._dp)
  isgood = isgood .and. (anint(10000._dp*dv_opt(2))  ==    629._dp)
  isgood = isgood .and. (anint(10000._dp*dv_opt(3))  ==    362._dp)
  isgood = isgood .and. (anint(10000._dp*dv_opt(4))  ==    -63._dp)
  isgood = isgood .and. (anint(10000._dp*dv_opt(5))  ==     26._dp)
  isgood = isgood .and. (anint(10000._dp*dv_opt(6))  == -76358._dp)
  isgood = isgood .and. (anint(10000._dp*dv_opt(7))  ==    718._dp)
  isgood = isgood .and. (anint(10000._dp*dv_opt(8))  ==   -528._dp)
  isgood = isgood .and. (anint(10000._dp*dv_opt(9))  ==  -1008._dp)
  isgood = isgood .and. (anint(10000._dp*dv_opt(10)) ==     37._dp)

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_dds double precision o.k.'
  else
     write(*,*) 'mo_dds double precision failed!'
  endif

end program main
