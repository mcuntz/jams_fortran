! ============================================================================
! Name        : main.f90
! Author      : Giovanni Dalmasso
! Version     : 1.0
! Copyright   : Department of Computational Hydrosystems CHS - UZF Leipzig
! Description : main test mo_ode_solver
! ============================================================================

program main_test_mo_ode_solver

    use mo_kind,        only : i4, sp, lgt
    use mo_eq_lv,       only : LV_eqn_sp
    use mo_ode_solver,   only : EULERdumb, RKdumb, RKas

    implicit none

    integer(i4)                             :: nEq
    real(sp)                                :: x1, x2, h, Tout, eps, hmin, test
    real(sp), dimension(:),     allocatable :: xi, xout
    real(sp), dimension(:,:),   allocatable :: yout

    logical(lgt)                            :: isgood

    isgood = .true.                 ! test variable
    test = 10._sp**4._sp            ! test variable

    nEq = 2_i4                      ! number of equations
    x1 = 0._sp                      ! intial time
    x2 = 10._sp                     ! final time
    h = 0.1_sp                      ! incremental time step
    hmin = h*0.00000000001_sp       ! minimum allowed stepsize
    eps = 10._sp**(-6._sp)          ! accuracy

    allocate( xi(nEq) )
    xi = (/ 1._sp, 1._sp /)         ! initial conditions

    call EULERdumb( xi, x1, x2, h, LV_eqn_sp, Tout, xout, yout )           ! EULER method
    if ( nint(xout(2)*test) .ne. 1000_i4 .or. any(nint(yout(:,2)*test) .ne. (/ 10000_i4, 13000_i4 /)) ) isgood = .false.    ! test

    call RKdumb( xi, x1, x2, h, LV_eqn_sp, Tout, xout, yout )                 ! RUNGE-KUTTA method
    if ( nint(xout(2)*test) .ne. 1000_i4 .or. any(nint(yout(:,2)*test) .ne. (/ 9850_i4, 12835_i4 /)) ) isgood = .false.     ! test

    call RKas( xi, x1, x2, h, LV_eqn_sp, Tout, xout, yout, hminIN=hmin, epsIN=eps )  ! RUNGE-KUTTA with adaptive stepsize control
    if ( nint(xout(2)*test) .ne. 1000_i4 .or. any(nint(yout(:,2)*test) .ne. (/ 9850_i4, 12835_i4 /)) ) isgood = .false.     ! test

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
