! ============================================================================
! Name        : main.f90
! Author      : Giovanni Dalmasso
! Date        : Mar 12, 2013
! Version     : 1.1
! Copyright   : Department of Computational Hydrosystems CHS - UZF Leipzig
! Description : main test mo_ode_solver
! ============================================================================

program main_test_mo_ode_solver

    use mo_kind,        only    : i4, sp
    use mo_eq_lv,       only    : LV_eqn_sp
    use mo_ode_solver,  only    : EULER, RK4, RK4as

    implicit none

    integer(i4)                             :: nEq      ! number of equation
    real(sp)                                :: x1       ! initial time point
    real(sp)                                :: x2       ! final time point
    real(sp)                                :: h        ! time step
    real(sp)                                :: hmin     ! minimum allowed stepsize
    real(sp)                                :: eps      ! accurancy
    real(sp)                                :: test     ! test variable
    real(sp), dimension(:),     allocatable :: yi       ! initial conditions
    real(sp), dimension(:),     allocatable :: xout     ! output time
    real(sp), dimension(:,:),   allocatable :: yout     ! output solutions

    logical                                 :: isgood

    isgood = .true.     ! test variable
    test   = 1e4_sp     ! test variable

    nEq    = 2_i4       ! number of equations
    x1     = 0._sp      ! intial time
    x2     = 10._sp     ! final time
    h      = 0.1_sp     ! incremental time step
    hmin   = h*1e-11_sp ! minimum allowed stepsize
    eps    = 1e-8_sp    ! accuracy

    allocate( yi(nEq) )
    yi = (/ 1._sp, 1._sp /)         ! initial conditions

    ! Euler method
    call Euler( yi, x1, x2, h, LV_eqn_sp, xout, yout )
    if ( nint(xout(2)*test) .ne. 1000_i4 .or. any(nint(yout(2,:)*test) .ne. (/ 10000_i4, 13000_i4 /)) ) isgood = .false.    ! test

    ! 4th order RUNGE-KUTTA method
    call RK4( yi, x1, x2, h, LV_eqn_sp, xout, yout )
    if ( nint(xout(2)*test) .ne. 1000_i4 .or. any(nint(yout(2,:)*test) .ne. (/ 9850_i4, 12835_i4 /)) ) isgood = .false.     ! test

    ! 4th order RUNGE-KUTTA Adaptive Step-size
    call RK4as( yi, x1, x2, h, LV_eqn_sp, xout, yout, hmin=hmin, eps=eps )
    if ( nint(xout(2)*test) .ne. 580_i4 .or. any(nint(yout(2,:)*test) .ne. (/ 9950_i4, 11687_i4 /)) ) isgood = .false.     ! test

    ! Is the program running properly?
    write(*,*) ''
    write(*,*) '------------------------------------------'
    if (isgood) then
        write(*,*) 'mo_ode_solver o.k.'
    else
        write(*,*) 'mo_ode_solver failed!'
    endif
    write(*,*) '-------------------END---------------------'

end program main_test_mo_ode_solver
