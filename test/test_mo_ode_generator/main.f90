!! =============================================================================
!! Name			: main.f90
!! Author      	: Giovanni Dalmasso
!! Contacts	  	: giovanni.dalmasso@ufz.de - http://www.ufz.de/index.php?en=30665
!!                Department of Computational Hydrosystems CHS - UZF Leipzig
!! Date	      	: Mar 20, 2014
!! Version     	: 1.0
!! Copyright   	: 2012-2014 Giovanni Dalmasso
!! Description 	: main for testing ODE generator
!! Note         : the following module from the fortran are needed:
!!                  - mo_constant.f90
!!                  - mo_kind.f90
!!                  - mo_nrutil.f90
!!                  - mo_ode_solver.f90
!! =============================================================================


program main

    use mo_kind,                    only    :   i4, dp
    use mo_ansi_colors, only: color, c_red, c_green
    use mo_ode_generator,           only    :   ode_generator

    implicit none

    integer(i4)                                 ::  n_reac                      ! number of reactants involved in the network
    real(dp)                                    ::  t_start, t_end              ! initial and final time step
    real(dp)                                    ::  dt                          ! initial size time step
    real(dp),   dimension(:),       allocatable ::  x_initial                   ! initial conditions
    real(dp),   dimension(:,:),     allocatable ::  para                        ! array of parameters
    real(dp),   dimension(:,:,:),   allocatable ::  ode_results                 ! results
    logical                                     ::  isgood

    !! ============================================================================

    !! number of reactants involved in the network
    n_reac = 2_i4
    allocate( x_initial(n_reac) )
    allocate( para(n_reac,n_reac) )

    !! initial conditions
    x_initial(:) = 0.5_dp

    !! initial and final time step
    t_start = 0.0_dp
    t_end = 5.0_dp

    !! initial size time step
    dt = 0.001_dp

    !! array of parameters
    para(:,:) = 1.0_dp

    !! calling mo_ode_generator to generate all the possible equations and solve it
    call ode_generator( n_reac, dt, t_start, t_end, x_initial, para, ode_results )

    isgood = .true.
    isgood = isgood .and. ( nint(ode_results(1,1,2)*1000.0_dp) .eq. 500 )
    isgood = isgood .and. ( all(nint(ode_results(:,5,2)*1000.0_dp) .eq. &
                            (/500, 500, 449, 453, 551, 578, 500, 502, 557, 557, 500, 500, 601, 604, 546, 557/)) )

    if (isgood) then
        write(*,*) 'mo_ode_generator ', color('o.k.', c_green)
    else
        write(*,*) 'mo_ode_generator ', color('failed!', c_red)
    endif


end program main
