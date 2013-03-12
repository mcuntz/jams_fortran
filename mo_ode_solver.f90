module mo_ode_solver

    ! This module provides a set of iterative methods for the approximation of solutions
    ! of Ordinary Differential Equations (ODE).
    !
    ! WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
    !    Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
    !    Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996
    !
    ! Written Numerical Recipes, 1996
    ! Modified Giovanni Dalmasso, Jul 2012 - adapted to FORTRAN_chs_lib structure
    !                                      - collected together different methodos
    !
    ! License
    ! -------
    ! This file is part of the UFZ Fortran library.
    !
    ! The UFZ Fortran library is free software: you can redistribute it and/or modify
    ! it under the terms of the GNU Lesser General Public License as published by
    ! the Free Software Foundation, either version 3 of the License, or
    ! (at your option) any later version.
    !
    ! The UFZ Fortran library is distributed in the hope that it will be useful,
    ! but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    ! GNU Lesser General Public License for more details.
    !
    ! You should have received a copy of the GNU Lesser General Public License
    ! along with the UFZ Fortran library. If not, see <http://www.gnu.org/licenses/>.
    !
    ! Copyright 2012 Giovanni Dalmasso
    !
    !
    ! Note on Numerical Recipes License
    ! ---------------------------------
    ! Be aware that some code is under the Numerical Recipes License 3rd
    ! edition <http://www.nr.com/aboutNR3license.html>
    !
    ! The Numerical Recipes Personal Single-User License lets you personally
    ! use Numerical Recipes code ("the code") on any number of computers,
    ! but only one computer at a time. You are not permitted to allow anyone
    ! else to access or use the code. You may, under this license, transfer
    ! precompiled, executable applications incorporating the code to other,
    ! unlicensed, persons, providing that (i) the application is
    ! noncommercial (i.e., does not involve the selling or licensing of the
    ! application for a fee), and (ii) the application was first developed,
    ! compiled, and successfully run by you, and (iii) the code is bound
    ! into the application in such a manner that it cannot be accessed as
    ! individual routines and cannot practicably be unbound and used in
    ! other programs. That is, under this license, your application user
    ! must not be able to use Numerical Recipes code as part of a program
    ! library or "mix and match" workbench.
    !
    ! Businesses and organizations that purchase the disk or code download,
    ! and that thus acquire one or more Numerical Recipes Personal
    ! Single-User Licenses, may permanently assign those licenses, in the
    ! number acquired, to individual employees. Such an assignment must be
    ! made before the code is first used and, once made, it is irrevocable
    ! and can not be transferred.
    !
    ! If you do not hold a Numerical Recipes License, this code is only for
    ! informational and educational purposes but cannot be used.

    implicit none

    private

    public :: Euler     ! Euler method with equal time-steps increments
    public :: RK4       ! Fourth-order Runge-Kutta method with equal time-steps increments
    public :: RK4as     ! Fourth-order Runge-Kutta method with adaptive stepsize control

    ! Interfaces for single and double precision routines
    interface Euler
        module procedure Euler_sp, Euler_dp
    end interface Euler

    ! ------------------------------------------------------------------

    !     NAME
    !         Euler, RK4

    !     PURPOSE
    !         Integration of Ordinary Differential Equations using Euler or fourth-order Runge-Kutta methods
    !         with fixed time-steps increments.
    !         Starting from N initial values vstart known at x1, use Euler or Runge-Kutta to advance nstep equal increments to x2.
    !         Results are stored in the variables xout and yout.

    !     CALLING SEQUENCE
    !         call Euler( vstart, x1, x2, h, derivs, xout, yout )
    !         call RK4( vstart, x1, x2, h, xout, yout )
    !
    !         --> see also example in test directory --> test/test_mo_ode_solver

    !     INTENT(IN)
    !         real(sp/dp),  dimension(:),   ::  vstart          ... initial conditions
    !                                                               N inital values known at the time x1 (given N ODEs)
    !
    !         real(sp/dp)                   ::  x1              ... initial time
    !
    !         real(sp/dp)                   ::  x2              ... final time
    !
    !         real(sp/dp)                   ::  h               ... size of the incremental time step (fixed)
    !
    !         interface                     ::  derivs_sp/dp    ... returns derivatives dydx of y at x

    !     INTENT(INOUT)
    !         none

    !     INTENT(OUT)
    !         real(sp/dp),  dimension(:),   allocatable    ::  xout     ... storage for outputs (time)
    !
    !         real(sp/dp),  dimension(:,:), allocatable    ::  yout     ... storage for outputs (incremented variables)
    !                                                                       dim 1 = function evaluation at any time point
    !                                                                       dim 2 = number of equations

    !     INTENT(IN), OPTIONAL
    !         none
    !
    !     INTENT(INOUT), OPTIONAL
    !         none

    !     INTENT(OUT), OPTIONAL
    !         none

    !     RESTRICTIONS
    !         The user has to supply the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.

    !     EXAMPLE
    !         --> see example in test directory --> test/test_mo_ode_solver

    !     LITERATURE
    !        1) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 77 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 1 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1992
    !        2) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1996

    !     HISTORY
    !         Written,  Giovanni Dalmasso, Jul 2012
    !         Modified, Giovanni Dalmasso, Mar 2013

    ! Interfaces for single and double precision routines
    interface RK4
        module procedure RK4_sp, RK4_dp
    end interface RK4

    ! ------------------------------------------------------------------

    !     NAME
    !         RK4as

    !     PURPOSE
    !         Integration of Ordinary Differential Equations using Runge-Kutta with adaptive stepsize control.
    !         Integrate the array of starting values ystart from x1 to x2 with accuracy eps, storing intermediate results
    !         in the module variables. h1 should be set as a guessed first stepsize,
    !         hmin as the minimum allowed stepsize (can be zero).
    !         On output ystart is replaced by values at the end of the integration interval.

    !     CALLING SEQUENCE
    !         call RK4as( ystart, x1, x2, h, derivs, xout, yout, hmin, eps )
    !
    !         --> see also example in test directory --> test/test_mo_ode_solver

    !     INTENT(IN)
    !         real(sp/dp),  dimension(:),   ::  ystart          ... initial conditions
    !                                                               N inital values known at the time x1 (given N ODEs)
    !
    !         real(sp/dp)                   ::  x1              ... initial time
    !
    !         real(sp/dp)                   ::  x2              ... final time
    !
    !         real(sp/dp)                   ::  h               ... guessed first stepsize
    !
    !         interface                     ::  derivs_sp/dp    ... returns derivatives dydx of y at x

    !     INTENT(INOUT)
    !         none

    !     INTENT(OUT)
    !         real(sp/dp),  dimension(:),   allocatable    ::  xout     ... storage for outputs (time)
    !
    !         real(sp/dp),  dimension(:,:), allocatable    ::  yout     ... storage for outputs (incremented variables)
    !                                                                       dim 1 = function evaluation at any time point
    !                                                                       dim 2 = number of equations

    !     INTENT(IN), OPTIONAL
    !         real(sp/dp)   ::  hmin        ... minimum allowed stepsize (can be zero)
    !                                           DEFAULT: 0.0
    !
    !         real(sp/dp)   ::  eps         ... accuracy (overall tolerance level)
    !                                           DEFAULT: 10.0**(-6.0)
    !
    !     INTENT(INOUT), OPTIONAL
    !         none

    !     INTENT(OUT), OPTIONAL
    !         none

    !     RESTRICTIONS
    !         The user has to supply the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.

    !     EXAMPLE
    !         --> see example in test directory --> test/test_mo_ode_solver

    !     LITERATURE
    !        1) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 77 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 1 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1992
    !        2) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1996

    !     HISTORY
    !         Written,  Giovanni Dalmasso, Jul 2012
    !         Modified, Giovanni Dalmasso, Mar 2013

    ! Interfaces for single and double precision routines
    interface RK4as
        module procedure RK4as_sp, RK4as_dp
    end interface RK4as

    interface CashKarpRK
        module procedure CashKarpRK_sp, CashKarpRK_dp
    end interface CashKarpRK

    ! ------------------------------------------------------------------

contains

    ! ------------------------------------------------------------------

    ! SINGLE PRECISION Euler
    subroutine Euler_sp( ystart, x1, x2, h, derivs, xout, yout )    ! all obligatory

        use mo_kind,    only : i4, sp
        use mo_nrutil,  only : nrerror

        implicit none

        ! Intent IN
        real(sp),   dimension(:),                  intent(in)  :: ystart   ! initial conditions
        real(sp),                                  intent(in)  :: x1, x2   ! initial and final time
        real(sp),                                  intent(in)  :: h        ! step size

        ! Intent OUT
        real(sp),   dimension(:),   allocatable,   intent(out) :: xout
        real(sp),   dimension(:,:), allocatable,   intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind,    only    : sp
                implicit none
                real(sp),                   intent(in)  :: x        ! time
                real(sp),   dimension(:),   intent(in)  :: y        ! unknowns of the equations
                real(sp),   dimension(:),   intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: j        ! counter
        integer(i4)                             :: nstep    ! nuber of steps
        real(sp)                                :: x
        real(sp),   dimension( size(ystart) )   :: dy, y

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout)       ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )                  ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        do j=1, nstep                   ! take nstep steps
            call derivs( x, y, dy )
            y = y + h*dy                ! step EULER
            if (x+h .eq. x) call nrerror( 'stepsize not significant in Euler_sp' )
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine Euler_sp


    ! DOUBLE PRECISION Euler
    subroutine Euler_dp( ystart, x1, x2, h, derivs, xout, yout )    ! all obligatory

        use mo_kind,    only : i4, dp
        use mo_nrutil,  only : nrerror

        implicit none

        ! Intent IN
        real(dp),   dimension(:),                  intent(in)  :: ystart   ! initial conditions
        real(dp),                                  intent(in)  :: x1, x2   ! initial and final time
        real(dp),                                  intent(in)  :: h        ! step size

        ! Intent OUT
        real(dp),   dimension(:),   allocatable,   intent(out) :: xout
        real(dp),   dimension(:,:), allocatable,   intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind,    only    : dp
                implicit none
                real(dp),                   intent(in)  :: x        ! time
                real(dp),   dimension(:),   intent(in)  :: y        ! unknowns of the equations
                real(dp),   dimension(:),   intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: j        ! counter
        integer(i4)                             :: nstep    ! nuber of steps
        real(dp)                                :: x
        real(dp),   dimension( size(ystart) )   :: dy, y

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout)       ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )                  ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        do j=1, nstep                   ! take nstep steps
            call derivs( x, y, dy )
            y = y + h*dy                ! step EULER
            if (x+h .eq. x) call nrerror( 'stepsize not significant in Euler_dp' )
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine Euler_dp


      ! SINGLE PRECISION 4th order RUNGE-KUTTA
    subroutine RK4_sp( ystart, x1, x2, h, derivs, xout, yout )    ! all obligatory

        use mo_kind,    only : i4, sp
        use mo_nrutil,  only : nrerror

        implicit none

        ! Intent IN
        real(sp),   dimension(:),                  intent(in)  :: ystart   ! initial conditions
        real(sp),                                  intent(in)  :: x1, x2   ! initial and final time
        real(sp),                                  intent(in)  :: h        ! step size

        ! Intent OUT
        real(sp),   dimension(:),   allocatable,   intent(out) :: xout
        real(sp),   dimension(:,:), allocatable,   intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind,    only    : sp
                implicit none
                real(sp),                   intent(in)  :: x        ! time
                real(sp),   dimension(:),   intent(in)  :: y        ! unknowns of the equations
                real(sp),   dimension(:),   intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: j            ! counter
        integer(i4)                             :: nstep        ! nuber of steps
        real(sp)                                :: x
        real(sp)                                :: hh, h6, xh
        real(sp),   dimension( size(ystart) )   :: dy, dyt, dym, y, yt

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout)       ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )                  ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        hh = h*0.5_sp
        h6 = h/6.0_sp
        xh = x+hh

        do j=1, nstep                          ! take nstep steps
            call derivs( x, y, dy )
            yt = y + hh*dy                     ! first step
            call derivs( xh, yt, dyt )         ! second step
            yt = y + hh*dyt
            call derivs( xh, yt, dym )         ! third step
            yt = y + h*dym
            dym = dyt + dym
            call derivs( x+h, yt, dyt )        ! fourth step
            y = y + h6*( dy+dyt+2.0_sp*dym )   ! accumulate increments with proper weights
            if (x+h .eq. x) call nrerror( 'stepsize not significant in RK4_sp' )
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine RK4_sp


    ! DOUBLE PRECISION 4th order RUNGE-KUTTA
    subroutine RK4_dp( ystart, x1, x2, h, derivs, xout, yout )    ! all obligatory

        use mo_kind,    only : i4, dp
        use mo_nrutil,  only : nrerror

        implicit none

        ! Intent IN
        real(dp),   dimension(:),                  intent(in)  :: ystart   ! initial conditions
        real(dp),                                  intent(in)  :: x1, x2   ! initial and final time
        real(dp),                                  intent(in)  :: h        ! step size

        ! Intent OUT
        real(dp),   dimension(:),   allocatable,   intent(out) :: xout
        real(dp),   dimension(:,:), allocatable,   intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind,    only    : dp
                implicit none
                real(dp),                   intent(in)  :: x        ! time
                real(dp),   dimension(:),   intent(in)  :: y        ! unknowns of the equations
                real(dp),   dimension(:),   intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: j            ! counter
        integer(i4)                             :: nstep        ! nuber of steps
        real(dp)                                :: x
        real(dp)                                :: hh, h6, xh
        real(dp),   dimension( size(ystart) )   :: dy, dyt, dym, y, yt

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout)       ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )                  ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        hh = h*0.5_dp
        h6 = h/6.0_dp
        xh = x+hh

        do j=1, nstep                          ! take nstep steps
            call derivs( x, y, dy )
            yt = y + hh*dy                     ! first step
            call derivs( xh, yt, dyt )         ! second step
            yt = y + hh*dyt
            call derivs( xh, yt, dym )         ! third step
            yt = y + h*dym
            dym = dyt + dym
            call derivs( x+h, yt, dyt )        ! fourth step
            y = y + h6*( dy+dyt+2.0_dp*dym )   ! accumulate increments with proper weights
            if (x+h .eq. x) call nrerror( 'stepsize not significant in RK4_dp' )
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine RK4_dp


    ! ------------------------------------------------------------------


    ! SINGLE PRECISION 4th order RUNGE-KUTTA Adaptive Step-size
    subroutine RK4as_sp( ystart, x1, x2, h, derivs, xout, yout, &   ! obligatory
        hmin, eps )                                                 ! optional

        use mo_kind,    only    : i4, sp
        use mo_nrutil,  only    : nrerror, reallocate

        implicit none

        ! Intent IN
        real(sp),                                   intent(in) :: x1, x2    ! initial and final time
        real(sp),                                   intent(in) :: h         ! guessed step size
        real(sp),   dimension(:),                   intent(in) :: ystart    ! initial conditions
        real(sp),                   optional,       intent(in) :: hmin
        real(sp),                   optional,       intent(in) :: eps

        ! Intent OUT
        real(sp),   dimension(:),   allocatable,    intent(out) :: xout
        real(sp),   dimension(:,:), allocatable,    intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind,    only    : sp
                implicit none
                real(sp),                   intent(in)  :: x        ! time
                real(sp),   dimension(:),   intent(in)  :: y        ! unknowns of the equations
                real(sp),   dimension(:),   intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                         :: nstep        ! nuber of steps
        integer(i4)                         :: nok, nbad    ! number of good and bad (but retried and fixed) steps taken
        integer(i4)                         :: kount        ! total number of saved steps
        real(sp)                            :: hIN, hdid, hnext, x, hminIN, epsIN, xsav, dxsav, errmax, htemp, xnew
        real(sp),   dimension(size(ystart)) :: dydx, y, yscal, yerr, ytemp
        ! Pointers
        real(sp),   dimension(:),   pointer :: xp
        real(sp),   dimension(:,:), pointer :: yp

        ! parameters
        integer(i4),    parameter   :: MAXstp = 1000000_i4  ! max nuber of steps
        real(sp),       parameter   :: safety = 0.9_sp, pgrow = -0.2_sp, pshrnk = -0.25_sp, &
            errcon = exp( (1._sp/pgrow)*log(5._sp/safety) )

        if( present(hmin) ) then
            hminIN = hmin
        else
            hminIN = 0.0_sp
        end if

        if( present(eps) ) then
            epsIN = eps
        else
            epsIN = 10._sp**(-6_i4)
        end if

        x = x1
        hIN = sign( h, x2-x1 )
        nok = 0_i4
        nbad = 0_i4
        kount = 0_i4
        y(:) = ystart(:)
        dxsav = tiny(1._sp)

        xsav = x-2.0_sp*dxsav                    ! assures storage of first step
        nullify( xp, yp )
        allocate( xp(256) )
        allocate( yp(size(xp), size(ystart)) )

        call save_a_step                        ! save initial step

        do nstep=1, MAXstp      ! take at most MAXstp steps

            call derivs( x, y, dydx )
            yscal(:) = abs( y(:) ) + abs( hIN*dydx(:) ) + tiny(1._sp)   ! scaling used to monitor accuracy --> CAN BE MODIFIED...

            if ( abs(x-xsav) .gt. abs(dxsav) ) call save_a_step         ! store intermediate results

            if ( (x+hIN-x2)*(x+hIN-x1) .gt. 0.0_sp ) hIN = x2-x         ! if stepsize can overshoot, decrease

            do
                call CashKarpRK( y, dydx, x, hIN, ytemp, yerr, derivs ) ! take a step
                errmax = maxval( abs(yerr(:)/yscal(:)) )/epsIN          ! evaluate accuracy
                if ( errmax .le. 1.0_sp )   exit                        ! step succeeded
                htemp = safety*hIN*(errmax**pshrnk)                     ! truncation error too large, reduce stepsize
                hIN = sign( max( abs(htemp), 0.1_sp*abs(hIN) ), hIN )   ! no more than a factor of 10
                xnew = x+hIN
                if ( xnew .eq. x )  call nrerror( 'hey!!! stepsize underflow in RK4as_sp' )
            end do

            if ( errmax .gt. errcon )   then                            ! compute size of next step
                hnext = safety*hIN*(errmax**pgrow)
            else                                                        ! no more than a factor of 5 increase
                hnext = 5.0_sp*hIN
            end if

            hdid = hIN
            x = x+hIN
            y(:) = ytemp(:)

            if ( hdid .eq. hIN ) then
                nok = nok+1_i4
            else
                nbad = nbad+1_i4
            end if

            if ( (x-x2)*(x2-x1) .ge. 0.0_sp )   then        ! are we done?!?!
                call save_a_step                            ! save final step
                allocate( xout(kount) )                     ! allocate storage for outputs
                xout(:) = xp(1:kount)
                allocate( yout(kount, size(yp,2)) )         ! allocate storage for outputs
                yout(:,:) = yp(1:kount, :)
                deallocate( xp, yp )                        ! clear out old stored variables
                return                                      ! normal exit
            end if

            if ( abs(hnext) .lt. hminIN ) call nrerror( 'WTF! ...stepsize smaller than minimum in RKas_sp!!!' )
            hIN = hnext

        end do

        call nrerror( 'dude, too many steps in RK4as_sp!!!' )

    contains

        subroutine save_a_step      ! --> like a macro in C
            kount = kount+1
            if ( kount .gt. size(xp) ) then
                xp=>reallocate( xp, 2*size(xp) )
                yp=>reallocate( yp, size(xp), size(yp,2) )
            end if
            xp(kount) = x
            yp(kount, :) = y(:)
            xsav = x
        end subroutine save_a_step

    end subroutine RK4as_sp


    ! DOUBLE PRECISION 4th order RUNGE-KUTTA Adaptive Step-size
    subroutine RK4as_dp( ystart, x1, x2, h, derivs, xout, yout, &   ! obligatory
        hmin, eps )                                                 ! optional

        use mo_kind,    only    : i4, dp
        use mo_nrutil,  only    : nrerror, reallocate

        implicit none

        ! Intent IN
        real(dp),                                   intent(in) :: x1, x2    ! initial and final time
        real(dp),                                   intent(in) :: h         ! guessed step size
        real(dp),   dimension(:),                   intent(in) :: ystart    ! initial conditions
        real(dp),                   optional,       intent(in) :: hmin
        real(dp),                   optional,       intent(in) :: eps

        ! Intent OUT
        real(dp),   dimension(:),   allocatable,    intent(out) :: xout
        real(dp),   dimension(:,:), allocatable,    intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind,    only    : dp
                implicit none
                real(dp),                   intent(in)  :: x        ! time
                real(dp),   dimension(:),   intent(in)  :: y        ! unknowns of the equations
                real(dp),   dimension(:),   intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                         :: nstep        ! nuber of steps
        integer(i4)                         :: nok, nbad    ! number of good and bad (but retried and fixed) steps taken
        integer(i4)                         :: kount        ! total number of saved steps
        real(dp)                            :: hIN, hdid, hnext, x, hminIN, epsIN, xsav, dxsav, errmax, htemp, xnew
        real(dp),   dimension(size(ystart)) :: dydx, y, yscal, yerr, ytemp
        ! Pointers
        real(dp),   dimension(:),   pointer :: xp
        real(dp),   dimension(:,:), pointer :: yp

        ! parameters
        integer(i4),    parameter   :: MAXstp = 1000000_i4  ! max nuber of steps
        real(dp),       parameter   :: safety = 0.9_dp, pgrow = -0.2_dp, pshrnk = -0.25_dp, &
            errcon = exp( (1._dp/pgrow)*log(5._dp/safety) )

        if( present(hmin) ) then
            hminIN = hmin
        else
            hminIN = 0.0_dp
        end if

        if( present(eps) ) then
            epsIN = eps
        else
            epsIN = 10._dp**(-6_i4)
        end if

        x = x1
        hIN = sign( h, x2-x1 )
        nok = 0_i4
        nbad = 0_i4
        kount = 0_i4
        y(:) = ystart(:)
        dxsav = tiny(1._dp)

        xsav = x-2.0_dp*dxsav                    ! assures storage of first step
        nullify( xp, yp )
        allocate( xp(256) )
        allocate( yp(size(xp), size(ystart)) )

        call save_a_step                        ! save initial step

        do nstep=1, MAXstp      ! take at most MAXstp steps

            call derivs( x, y, dydx )
            yscal(:) = abs( y(:) ) + abs( hIN*dydx(:) ) + tiny(1._dp)   ! scaling used to monitor accuracy --> CAN BE MODIFIED...

            if ( abs(x-xsav) .gt. abs(dxsav) ) call save_a_step         ! store intermediate results

            if ( (x+hIN-x2)*(x+hIN-x1) .gt. 0.0_dp ) hIN = x2-x         ! if stepsize can overshoot, decrease

            do
                call CashKarpRK( y, dydx, x, hIN, ytemp, yerr, derivs ) ! take a step
                errmax = maxval( abs(yerr(:)/yscal(:)) )/epsIN          ! evaluate accuracy
                if ( errmax .le. 1.0_dp )   exit                        ! step succeeded
                htemp = safety*hIN*(errmax**pshrnk)                     ! truncation error too large, reduce stepsize
                hIN = sign( max( abs(htemp), 0.1_dp*abs(hIN) ), hIN )   ! no more than a factor of 10
                xnew = x+hIN
                if ( xnew .eq. x )  call nrerror( 'hey!!! stepsize underflow in RK4as_dp' )
            end do

            if ( errmax .gt. errcon )   then                            ! compute size of next step
                hnext = safety*hIN*(errmax**pgrow)
            else                                                        ! no more than a factor of 5 increase
                hnext = 5.0_dp*hIN
            end if

            hdid = hIN
            x = x+hIN
            y(:) = ytemp(:)

            if ( hdid .eq. hIN ) then
                nok = nok+1_i4
            else
                nbad = nbad+1_i4
            end if

            if ( (x-x2)*(x2-x1) .ge. 0.0_dp )   then        ! are we done?!?!
                call save_a_step                            ! save final step
                allocate( xout(kount) )                     ! allocate storage for outputs
                xout(:) = xp(1:kount)
                allocate( yout(kount, size(yp,2)) )         ! allocate storage for outputs
                yout(:,:) = yp(1:kount, :)
                deallocate( xp, yp )                        ! clear out old stored variables
                return                                      ! normal exit
            end if

            if ( abs(hnext) .lt. hminIN ) call nrerror( 'WTF! ...stepsize smaller than minimum in RKas_dp!!!' )
            hIN = hnext

        end do

        call nrerror( 'dude, too many steps in RK4as_dp!!!' )

    contains

        subroutine save_a_step      ! --> like a macro in C
            kount = kount+1
            if ( kount .gt. size(xp) ) then
                xp=>reallocate( xp, 2*size(xp) )
                yp=>reallocate( yp, size(xp), size(yp,2) )
            end if
            xp(kount) = x
            yp(kount, :) = y(:)
            xsav = x
        end subroutine save_a_step

    end subroutine RK4as_dp


    ! ============================================================================
    !   PRIVATE METHODS
    ! ============================================================================

    ! SINGLE PRECISION CASH-KARP RUNGE-KUTTA step
    subroutine CashKarpRK_sp( y, dydx, x, h, yout, yerr, derivs )
        ! Given values for N variables y and their derivatives dydx known at x, use the fifth-order
        ! Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
        ! the incremented variables as yout. Also return an estimate of the local truncation error
        ! in yout using the embedded fourth order method.

        use mo_kind,    only    : i4, sp
        use mo_nrutil,  only    : assert_eq

        implicit none

        ! Intent IN
        real(sp),   dimension(:),   intent(in)  :: y, dydx
        real(sp),                   intent(in)  :: x, h

        ! Intent OUT
        real(sp),   dimension(:),   intent(out) :: yout, yerr

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind,    only    : sp
                implicit none
                real(sp),                   intent(in)  :: x
                real(sp),   dimension(:),   intent(in)  :: y
                real(sp),   dimension(:),   intent(out) :: dydx
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: ndum
        real(sp),   dimension(:),   allocatable :: ak2, ak3, ak4, ak5, ak6, ytemp

        ! parameters
        real(sp),   parameter   :: A2=.2_sp, A3=.3_sp, A4=.6_sp, A5=1._sp, A6=.875_sp, B21=.2_sp, B31=3._sp/40._sp, &
            B32=9._sp/40._sp, B41=.3_sp, B42=-.9_sp, B43=1.2_sp, B51=-11._sp/54._sp, B52=2.5_sp, &
            B53=-70._sp/27._sp, B54=35._sp/27._sp, B61=1631._sp/55296._sp, B62=175._sp/512._sp, &
            B63=575._sp/13824._sp, B64=44275._sp/110592._sp, B65=253._sp/4096._sp, &
            C1=37._sp/378._sp, C3=250._sp/621._sp, C4=125._sp/594._sp, C6=512._sp/1771._sp, &
            DC1=C1-2825._sp/27648._sp, DC3=C3-18575._sp/48384._sp, &
            DC4=C4-13525._sp/55296._sp, DC5=-277._sp/14336._sp, DC6=C6-.25_sp

        ndum = assert_eq( size(y), size(dydx) ,size(yout) ,size(yerr) , 'CashKarpRK_sp' )

        allocate( ak2(ndum), ak3(ndum), ak4(ndum), ak5(ndum), ak6(ndum), ytemp(ndum) )

        ytemp = y + B21*h*dydx                                  ! first step
        call derivs( x+A2*h, ytemp, ak2 )                       ! second step
        ytemp = y + h*(B31*dydx+B32*ak2)
        call derivs( x+A3*h, ytemp, ak3 )                       ! third step
        ytemp = y + h*(B41*dydx+B42*ak2+B43*ak3)
        call derivs( x+A4*h, ytemp, ak4 )                       ! fourth step
        ytemp = y + h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
        call derivs( x+A5*h, ytemp, ak5 )                       ! fifth step
        ytemp = y + h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
        call derivs( x+A6*h, ytemp, ak6 )                       ! sixth step
        yout = y + h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)             ! accumulate increments with proper weights
        yerr = h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)

    end subroutine CashKarpRK_sp


    ! DOUBLE PRECISION CASH-KARP RUNGE-KUTTA step
    subroutine CashKarpRK_dp( y, dydx, x, h, yout, yerr, derivs )
        ! Given values for N variables y and their derivatives dydx known at x, use the fifth-order
        ! Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
        ! the incremented variables as yout. Also return an estimate of the local truncation error
        ! in yout using the embedded fourth order method.

        use mo_kind,    only    : i4, dp
        use mo_nrutil,  only    : assert_eq

        implicit none

        ! Intent IN
        real(dp),   dimension(:),   intent(in)  :: y, dydx
        real(dp),                   intent(in)  :: x, h

        ! Intent OUT
        real(dp),   dimension(:),   intent(out) :: yout, yerr

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind,    only    : dp
                implicit none
                real(dp),                   intent(in)  :: x
                real(dp),   dimension(:),   intent(in)  :: y
                real(dp),   dimension(:),   intent(out) :: dydx
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: ndum
        real(dp),   dimension(:),   allocatable :: ak2, ak3, ak4, ak5, ak6, ytemp

        ! parameters
        real(dp),   parameter   :: A2=.2_dp, A3=.3_dp, A4=.6_dp, A5=1._dp, A6=.875_dp, B21=.2_dp, B31=3._dp/40._dp, &
            B32=9._dp/40._dp, B41=.3_dp, B42=-.9_dp, B43=1.2_dp, B51=-11._dp/54._dp, B52=2.5_dp, &
            B53=-70._dp/27._dp, B54=35._dp/27._dp, B61=1631._dp/55296._dp, B62=175._dp/512._dp, &
            B63=575._dp/13824._dp, B64=44275._dp/110592._dp, B65=253._dp/4096._dp, &
            C1=37._dp/378._dp, C3=250._dp/621._dp, C4=125._dp/594._dp, C6=512._dp/1771._dp, &
            DC1=C1-2825._dp/27648._dp, DC3=C3-18575._dp/48384._dp, &
            DC4=C4-13525._dp/55296._dp, DC5=-277._dp/14336._dp, DC6=C6-.25_dp

        ndum = assert_eq( size(y), size(dydx) ,size(yout) ,size(yerr) , 'CashKarpRK_dp' )

        allocate( ak2(ndum), ak3(ndum), ak4(ndum), ak5(ndum), ak6(ndum), ytemp(ndum) )

        ytemp = y + B21*h*dydx                                  ! first step
        call derivs( x+A2*h, ytemp, ak2 )                       ! second step
        ytemp = y + h*(B31*dydx+B32*ak2)
        call derivs( x+A3*h, ytemp, ak3 )                       ! third step
        ytemp = y + h*(B41*dydx+B42*ak2+B43*ak3)
        call derivs( x+A4*h, ytemp, ak4 )                       ! fourth step
        ytemp = y + h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
        call derivs( x+A5*h, ytemp, ak5 )                       ! fifth step
        ytemp = y + h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
        call derivs( x+A6*h, ytemp, ak6 )                       ! sixth step
        yout = y + h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)             ! accumulate increments with proper weights
        yerr = h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)

    end subroutine CashKarpRK_dp

end module mo_ode_solver
