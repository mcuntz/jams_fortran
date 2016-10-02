! ============================================================================
! Name        : main.f90
! Author      : Giovanni Dalmasso
! Date        : Mar 12, 2013
! Version     : 1.1
! Copyright   : Department of Computational Hydrosystems CHS - UZF Leipzig
! Description : main test mo_ode_solver
! ============================================================================

program main_test_mo_ode_solver

    use mo_kind,        only    : i4, sp, dp
    use mo_eq_lv,       only    : LV_eqn_sp, LV_eqn_dp, LV_eqn_para_sp, LV_eqn_para_dp,&
                                  deriv2_testRB_sp, jacobn2_testRB_sp,&
                                  deriv2_testRB_dp, jacobn2_testRB_dp,&
                                  deriv2_testRB_para_sp, jacobn2_testRB_para_sp,&
                                  deriv2_testRB_para_dp, jacobn2_testRB_para_dp
    use mo_ode_solver,  only    : EULER, RK4, RK4as, RBstiff

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
    real(sp), dimension(3)                  :: para     ! parameter for derivs
    
    real(dp)                                :: dx1       ! initial time point
    real(dp)                                :: dx2       ! final time point
    real(dp)                                :: dh        ! time step
    real(dp)                                :: dhmin     ! minimum allowed stepsize
    real(dp)                                :: deps      ! accurancy
    real(dp)                                :: dtest     ! test variable
    real(dp), dimension(:),     allocatable :: dyi       ! initial conditions
    real(dp), dimension(:),     allocatable :: dxout     ! output time
    real(dp), dimension(:,:),   allocatable :: dyout     ! output solutions
    real(dp), dimension(3)                  :: dpara     ! parameter for derivs

    logical                                 :: isgood


    isgood = .true.       ! test variable

    ! sp
    test   = 1.0e4_sp     ! test variable

    nEq    = 2_i4         ! number of equations
    x1     = 0._sp        ! intial time
    x2     = 10._sp       ! final time
    h      = 0.1_sp       ! incremental time step
    hmin   = h*1.0e-11_sp ! minimum allowed stepsize
    eps    = 1e-8_sp      ! accuracy

    para   = (/0.5_sp, 2.0_sp, -1.0_sp/)

    allocate( yi(nEq) )
    yi = (/ 1._sp, 1._sp /)         ! initial conditions

    ! Euler method
    call Euler( yi, x1, x2, h, LV_eqn_sp, xout, yout )
    if ( nint(xout(2)*test) .ne. 1000_i4 .or. any(nint(yout(2,:)*test) .ne. (/ 10000_i4, 13000_i4 /)) ) isgood = .false.    ! test
    if (.not. isgood) write(*,*) "euler sp"

    ! 4th order RUNGE-KUTTA method
    call RK4( yi, x1, x2, h, LV_eqn_sp, xout, yout )
    if ( nint(xout(2)*test) .ne. 1000_i4 .or. any(nint(yout(2,:)*test) .ne. (/ 9850_i4, 12835_i4 /)) ) isgood = .false.     ! test
    if (.not. isgood) write(*,*) "rk sp"

    ! 4th order RUNGE-KUTTA Adaptive Step-size
    call RK4as( yi, x1, x2, h, LV_eqn_sp, xout, yout, hmin=hmin, eps=eps )
	
#if defined (pgiFortran) || defined (INTEL)
    ! PGI compiler other rounding errors in single precision, Intel in release mode
    if ( nint(xout(2)*test*1.0e-1_sp) .ne. 58_i4 .or. &
         any(nint(yout(2,:)*test*1.0e-1_sp) .ne. (/ 995_i4, 1169_i4 /)) ) isgood = .false.
#else
    if ( nint(xout(2)*test) .ne. 580_i4 .or. any(nint(yout(2,:)*test) .ne. (/ 9950_i4, 11687_i4 /)) ) isgood = .false.
#endif
    ! test
    if (.not. isgood) write(*,*) "rk-as sp"

    ! Rosenbrock Method
    call RBstiff( (/1.0_sp,1.0_sp,0.0_sp/), 0.0_sp, 50.0_sp, 2.9e-4_sp, deriv2_testRB_sp,&
                    jacobn2_testRB_sp, xout, yout, hmin=hmin, eps=1e-4_sp )
    if (    nint(xout(2)*test) .ne. 3_i4 .or.&
            any(nint(yout(2,:)*test) .ne. (/ 10000_i4, 10000_i4, 0_i4 /)) )&
       isgood = .false.    ! test
    if (.not. isgood) write(*,*) "rb-stiff sp"


    !with parameters

    ! Euler method with parameter
    call Euler( yi, x1, x2, h, LV_eqn_para_sp, para, xout, yout )
    if ( nint(xout(2)*test) .ne. 1000_i4 .or. any(nint(yout(2,:)*test) .ne. (/ 10000_i4, 13000_i4 /)) ) isgood = .false.    ! test
    if (.not. isgood) write(*,*) "euler sp para"

    ! 4th order RUNGE-KUTTA method with parameter
    call RK4( yi, x1, x2, h, LV_eqn_para_sp, para, xout, yout )
#if defined (pgiFortran) || defined (INTEL)
    if ( nint(xout(2)*test*1e-1) .ne. 100_i4 .or. any(nint(yout(2,:)*test*1e-1) .ne. (/ 985_i4, 1284_i4 /)) ) isgood = .false.     ! test
#else
    if ( nint(xout(2)*test) .ne. 1000_i4 .or. any(nint(yout(2,:)*test) .ne. (/ 9850_i4, 12835_i4 /)) ) isgood = .false.     ! test
#endif
    if (.not. isgood) write(*,*) "rk sp para"

    ! 4th order RUNGE-KUTTA Adaptive Step-size with parameter
    call RK4as( yi, x1, x2, h, LV_eqn_para_sp, para, xout, yout, hmin=hmin, eps=eps )
#if defined (pgiFortran) || defined (INTEL)
    if ( nint(xout(2)*test*1e-1) .ne. 58_i4 .or. any(nint(yout(2,:)*test*1e-1) .ne. (/ 995_i4, 1169_i4 /)) ) isgood = .false.     ! test
#else
    if ( nint(xout(2)*test) .ne. 580_i4 .or. any(nint(yout(2,:)*test) .ne. (/ 9950_i4, 11687_i4 /)) ) isgood = .false.     ! test
#endif
    if (.not. isgood) write(*,*) "rk-as sp para"

    ! Rosenbrock Method with parameter
    call RBstiff( (/1.0_sp,1.0_sp,0.0_sp/), 0.0_sp, 50.0_sp, 2.9e-4_sp, deriv2_testRB_para_sp,&
                    jacobn2_testRB_para_sp, (/1000.0_sp/), xout, yout, hmin=hmin, eps=1e-4_sp )
#if defined (pgiFortran) || defined (INTEL)
    if (    nint(xout(2)*test*1e-1) .ne. 0_i4 .or.&
            any(nint(yout(2,:)*test*1e-1) .ne. (/ 1000_i4, 1000_i4, 0_i4 /)) )&
       isgood = .false.    ! test
#else
    if (    nint(xout(2)*test) .ne. 3_i4 .or.&
            any(nint(yout(2,:)*test) .ne. (/ 10000_i4, 10000_i4, 0_i4 /)) )&
       isgood = .false.    ! test
#endif
    if (.not. isgood) write(*,*) "rb-stiff sp"


!dp
    dtest   = 1e8_dp     ! test variable

    nEq    = 2_i4       ! number of equations
    dx1     = 0._dp      ! intial time
    dx2     = 10._dp     ! final time
    dh      = 0.1_dp     ! incremental time step
    dhmin   = h*1e-11_dp ! minimum allowed stepsize
    deps    = 1e-8_dp    ! accuracy

    dpara   = (/0.5_dp, 2.0_dp, -1.0_dp/)

    allocate( dyi(nEq) )
    dyi = (/ 1._dp, 1._dp /)         ! initial conditions

    ! Euler method
    call Euler( dyi, dx1, dx2, dh, LV_eqn_dp, dxout, dyout )
    if ( nint(dxout(2)*dtest) .ne. 10000000_i4 .or. &
         any(nint(dyout(2,:)*dtest) .ne. (/ 100000000_i4, 130000000_i4 /)) ) isgood = .false.  ! test
    if (.not. isgood) write(*,*) "euler dp"

    ! 4th order RUNGE-KUTTA method
    call RK4( dyi, dx1, dx2, dh, LV_eqn_dp, dxout, dyout )
    if ( nint(dxout(2)*dtest) .ne. 10000000_i4 .or. &
         any(nint(dyout(2,:)*dtest) .ne. (/ 98503750_i4, 128353750_i4 /)) ) isgood = .false.   ! test
    if (.not. isgood) write(*,*) "rk dp"

    ! 4th order RUNGE-KUTTA Adaptive Step-size
    call RK4as( dyi, dx1, dx2, dh, LV_eqn_dp, dxout, dyout, hmin=dhmin, eps=deps )
    if ( nint(dxout(2)*dtest) .ne. 6008771_i4 .or. &
         any(nint(dyout(2,:)*dtest) .ne. (/ 99458909_i4, 117452697_i4 /)) ) isgood = .false.    ! test
    if (.not. isgood) write(*,*) "rk-as dp"

    ! Rosenbrock Method
    call RBstiff( (/1.0_dp,1.0_dp,0.0_dp/), 0.0_dp, 50.0_dp, 2.9e-4_dp, deriv2_testRB_dp,&
                    jacobn2_testRB_dp, dxout, dyout, hmin=dhmin, eps=1e-4_dp )
    if (    nint(dxout(2)*dtest) .ne. 29000_i4 .or.&
            any(nint(dyout(2,:)*dtest) .ne. (/ 99999663_i4, 100000100_i4, -237_i4 /)) )&
       isgood = .false.    ! test
    if (.not. isgood) write(*,*) "rb-stiff dp"

    !with parameters

    ! Euler method with parameters
    call Euler( dyi, dx1, dx2, dh, LV_eqn_para_dp, dpara, dxout, dyout )
    if ( nint(dxout(2)*dtest) .ne. 10000000_i4 .or. &
         any(nint(dyout(2,:)*dtest) .ne. (/ 100000000_i4, 130000000_i4 /)) ) isgood = .false.  ! test
    if (.not. isgood) write(*,*) "euler dp para"

!    write(*,*) dxout(2), dyout(2,:)
!    write(*,*) nint(dxout(2)*dtest), nint(dyout(2,:)*dtest)

    ! 4th order RUNGE-KUTTA method with parameters
    call RK4( dyi, dx1, dx2, dh, LV_eqn_para_dp, dpara, dxout, dyout )
    if ( nint(dxout(2)*dtest) .ne. 10000000_i4 .or. &
         any(nint(dyout(2,:)*dtest) .ne. (/ 98503750_i4, 128353750_i4 /)) ) isgood = .false.   ! test
    if (.not. isgood) write(*,*) "rk dp para"

!    write(*,*) dxout(2), dyout(2,:)
!    write(*,*) nint(dxout(2)*dtest), nint(dyout(2,:)*dtest)

    ! 4th order RUNGE-KUTTA Adaptive Step-size with parameters
    call RK4as( dyi, dx1, dx2, dh, LV_eqn_para_dp, dpara, dxout, dyout, hmin=dhmin, eps=deps )
    if ( nint(dxout(2)*dtest) .ne. 6008771_i4 .or. &
         any(nint(dyout(2,:)*dtest) .ne. (/ 99458909_i4, 117452697_i4 /)) ) isgood = .false.    ! test
    if (.not. isgood) write(*,*) "rk-as dp para"

!    write(*,*) dxout(2), dyout(2,:)
!    write(*,*) nint(dxout(2)*dtest), nint(dyout(2,:)*dtest)

    ! Rosenbrock Method with parameter
    call RBstiff( (/1.0_dp,1.0_dp,0.0_dp/), 0.0_dp, 50.0_dp, 2.9e-4_dp, deriv2_testRB_para_dp,&
                    jacobn2_testRB_para_dp, (/1000.0_dp/), dxout, dyout, hmin=dhmin, eps=1e-4_dp )
    if (    nint(dxout(2)*dtest) .ne. 29000_i4 .or.&
            any(nint(dyout(2,:)*dtest) .ne. (/ 99999663_i4, 100000100_i4, -237_i4 /)) )&
       isgood = .false.    ! test
    if (.not. isgood) write(*,*) "rb-stiff dp para"

!    write(*,*) dxout(2), dyout(2,:)
!    write(*,*) nint(dxout(2)*dtest), nint(dyout(2,:)*dtest)


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
