!> \file mo_ode_solver.f90

!> \brief This module provides a set of iterative methods for the approximation of solutions
!> of Ordinary Differential Equations (ODE).

!> \details
!> It includes the possibilities to integrate a system of Ordinary Differential Equations using Euler,
!> a fourth-order Runge-Kutta with fixed time-steps increments or a fourth-order Runge-Kutta with
!> adaptive stepsize control.

!> \authors Giovanni Dalmasso, Sebastian Mueller
!> \date Mar 2015
module mo_ode_solver

    ! This module provides a set of iterative methods for the approximation of solutions
    ! of Ordinary Differential Equations (ODE).

    ! Written  Giovanni Dalmasso, Jul 2012
    ! Modified Giovanni Dalmasso, Mar 2013 - adapted to FORTRAN_chs_lib structure
    !                                      - collected together different methods
    !                                      - speeded up
    !                             Apr 2013 - added documentation
    ! Modified Sebastian Mueller, Jan 2015 - added a parameter input for the derivatives to use the solver dynamically
    !                             Mar 2015 - added RBstiff as a solver for stiff ODEs


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
    ! along with the UFZ Fortran library (cf. gpl.txt and lgpl.txt).
    ! If not, see <http://www.gnu.org/licenses/>.
    !
    ! Copyright 2012-2013 Giovanni Dalmasso


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

    public :: Euler     ! Euler method with equal time-steps increments
    public :: RK4       ! Fourth-order Runge-Kutta method with equal time-steps increments
    public :: RK4as     ! Fourth-order Runge-Kutta method with adaptive stepsize control
    public :: RBstiff   ! Rosenbrock Method for stiff ODEs with adaptive stepsize control

    ! ------------------------------------------------------------------

    !     NAME
    !         Euler

    !     PURPOSE
    !         Integration of Ordinary Differential Equations using Euler method.
    !         Starting from N initial values vstart known at x1, use Euler to advance nstep equal increments to x2.
    !         Results are stored in the variables xout and yout.

    !>        \brief Euler.

    !>        \details Starting from $N$ initial values $v_{start}$ known at $x_1$,
    !>                 use Euler to advance $n$-steps equal increments to $x_2$.
    !>                 Results are stored in the variables $x_{out}$ and $y_{out}$.
    !>                 If you use "para" as parameters for derivs, the interface for derivs has to look like:
    !>
    !>                            interface
    !>                                subroutine derivs( x, y, para, dydx )
    !>                                    use mo_kind, only: sp/dp
    !>                                    implicit none
    !>                                    real(sp/dp),                 intent(in)  :: x        ! time
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: para     ! parameter for the derivatives
    !>                                    real(sp/dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
    !>                                end subroutine derivs
    !>                            end interface
    !>
    !>                 Elsewise "para" must be left out:
    !>
    !>                            interface
    !>                                subroutine derivs( x, y, dydx )
    !>                                    use mo_kind, only: sp/dp
    !>                                    implicit none
    !>                                    real(sp/dp),                 intent(in)  :: x        ! time
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
    !>                                    real(sp/dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
    !>                                end subroutine derivs
    !>                            end interface

    !     INTENT(IN)
    !>        \param[in] "real(sp/dp),  dimension(:)    ::  vstart"          initial conditions.
    !>                                                                       $N$ inital values known at the time $x_1$
    !>                                                                       (given $N$ ODEs)

    !>        \param[in] "real(sp/dp)                   ::  x1"              initial time.
    !>        \param[in] "real(sp/dp)                   ::  x2"              final time.
    !>        \param[in] "real(sp/dp)                   ::  h"               size of the incremental time step (fixed).
    !>        \param[in] "interface                     ::  derivs_sp/dp"    returns derivatives $dydx$ of $y$ at $x$.

    !     INTENT(INOUT)
    !         None

    !     INTENT(OUT)
    !>        \param[out] "real(sp/dp),  dimension(:), allocatable :: xout"  storage for outputs (time).
    !>        \param[out] "real(sp/dp),  dimension(:), allocatable :: yout"  storage for outputs
    !>                                                                       (incremented variables).
    !>                                                                       dim 1 = function evaluation
    !>                                                                               at any time point.
    !>                                                                       dim 2 = number of equations.

    !     INTENT(IN), OPTIONAL
    !>        \param[in] "real(sp/dp),  dimension(:)    ::  para"           parameter for the derivatives $dydx$

    !     INTENT(INOUT), OPTIONAL
    !         None

    !     INTENT(OUT), OPTIONAL
    !         None

    !     RESTRICTIONS
    !>       \note The user has to supply the subroutine derivs(x,y,dydx), which returns derivatives $dydx$ at $x$.

    !     EXAMPLE
    !         call Euler( vstart, x1, x2, h, derivs, xout, yout )
    !         --> see example in test directory --> test/test_mo_ode_solver

    !     LITERATURE
    !         http://en.wikipedia.org/wiki/Euler_method

    !     HISTORY
    !>        \author Giovanni Dalmasso
    !>        \date Jul 2012
    !         Modified, Giovanni Dalmasso, Mar 2013
    !         Modified, Sebastian Mueller, Jan 2015
        interface Euler
        module procedure Euler_sp, Euler_dp, Euler_para_sp, Euler_para_dp
    end interface Euler

    ! ------------------------------------------------------------------

    !     NAME
    !         RK4

    !     PURPOSE
    !         Integration of Ordinary Differential Equations using a fourth-order Runge-Kutta method
    !         with fixed time-steps increments.
    !         Starting from N initial values vstart known at x1, use Euler to advance nstep equal increments to x2.
    !         Results are stored in the variables xout and yout.

    !>        \brief fourth-order Runge-Kutta.

    !>        \details Starting from $N$ initial values $v_{start}$ known at $x_1$,
    !>                 use a fourth-order Runge-Kutta with fixed time-steps increments
    !>                 to advance $n$-steps equal increments to $x_2$.
    !>                 Results are stored in the variables $x_{out}$ and $y_{out}$.
    !>                 If you use "para" as parameters for derivs, the interface for derivs has to look like:
    !>
    !>                            interface
    !>                                subroutine derivs( x, y, para, dydx )
    !>                                    use mo_kind, only: sp/dp
    !>                                    implicit none
    !>                                    real(sp/dp),                 intent(in)  :: x        ! time
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: para     ! parameter for the derivatives
    !>                                    real(sp/dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
    !>                                end subroutine derivs
    !>                            end interface
    !>
    !>                 Elsewise "para" must be left out:
    !>
    !>                            interface
    !>                                subroutine derivs( x, y, dydx )
    !>                                    use mo_kind, only: sp/dp
    !>                                    implicit none
    !>                                    real(sp/dp),                 intent(in)  :: x        ! time
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
    !>                                    real(sp/dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
    !>                                end subroutine derivs
    !>                            end interface

    !     INTENT(IN)
    !>        \param[in] "real(sp/dp),  dimension(:)    ::  vstart"          initial conditions.
    !>                                                                       $N$ inital values known at the time $x_1$
    !>                                                                       (given $N$ ODEs)

    !>        \param[in] "real(sp/dp)                   ::  x1"              initial time.
    !>        \param[in] "real(sp/dp)                   ::  x2"              final time.
    !>        \param[in] "real(sp/dp)                   ::  h"               size of the incremental time step (fixed).
    !>        \param[in] "interface                     ::  derivs_sp/dp"    returns derivatives $dydx$ of $y$ at $x$.

    !     INTENT(INOUT)
    !         None

    !     INTENT(OUT)
    !>        \param[out] "real(sp/dp),  dimension(:), allocatable :: xout"  storage for outputs (time).
    !>        \param[out] "real(sp/dp),  dimension(:), allocatable :: yout"  storage for outputs
    !>                                                                       (incremented variables).
    !>                                                                       dim 1 = function evaluation
    !>                                                                               at any time point.
    !>                                                                       dim 2 = number of equations.

    !     INTENT(IN), OPTIONAL
    !>        \param[in] "real(sp/dp),  dimension(:)    ::  para"           parameter for the derivatives $dydx$

    !     INTENT(INOUT), OPTIONAL
    !         None

    !     INTENT(OUT), OPTIONAL
    !         None

    !     RESTRICTIONS
    !>       \note The user has to supply the subroutine derivs(x,y,dydx), which returns derivatives $dydx$ at $x$.

    !     EXAMPLE
    !           call RK4( vstart, x1, x2, h, derivs, xout, yout )
    !           --> see example in test directory --> test/test_mo_ode_solver

    !     LITERATURE
    !        1) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 77 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 1 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1992
    !        2) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1996

    !     HISTORY
    !>        \author Giovanni Dalmasso
    !>        \date Jul 2012
    !         Modified, Giovanni Dalmasso, Mar 2013
    !         Modified, Sebastian Mueller, Jan 2015
    interface RK4
        module procedure RK4_sp, RK4_dp, RK4_para_sp, RK4_para_dp
    end interface RK4

    ! ------------------------------------------------------------------

    !     NAME
    !         RK4as

    !     PURPOSE
    !         Integration of Ordinary Differential Equations using a fourth-order Runge-Kutta method
    !         with adaptive stepsize control.
    !         Integrate the array of starting values ystart from x1 to x2 with accuracy eps, storing intermediate results
    !         in the module variables. h1 should be set as a guessed first stepsize,
    !         hmin as the minimum allowed stepsize (can be zero).
    !         On output ystart is replaced by values at the end of the integration interval.

    !>        \brief fourth-order Runge-Kutta with adaptive stepsize control.

    !>        \details Integrate the array of starting values $y_{start} from $x_1$ to $x_2$ with accuracy $\varepsilon$,
    !>                 storing intermediate results in the module variables.
    !>                 $h_1$ should be set as a guessed first stepsize, $h_{min} as the minimum allowed stepsize (can be zero).
    !>                 On output $y_{start}$ is replaced by values at the end of the integration interval.
    !>                 If you use "para" as parameters for derivs, the interface for derivs has to look like:
    !>
    !>                            interface
    !>                                subroutine derivs( x, y, para, dydx )
    !>                                    use mo_kind, only: sp/dp
    !>                                    implicit none
    !>                                    real(sp/dp),                 intent(in)  :: x        ! time
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: para     ! parameter for the derivatives
    !>                                    real(sp/dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
    !>                                end subroutine derivs
    !>                            end interface
    !>
    !>                 Elsewise "para" must be left out:
    !>
    !>                            interface
    !>                                subroutine derivs( x, y, dydx )
    !>                                    use mo_kind, only: sp/dp
    !>                                    implicit none
    !>                                    real(sp/dp),                 intent(in)  :: x        ! time
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
    !>                                    real(sp/dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
    !>                                end subroutine derivs
    !>                            end interface

    !     INTENT(IN)
    !>        \param[in] "real(sp/dp),  dimension(:)    ::  vstart"          initial conditions.
    !>                                                                       $N$ inital values known at the time $x_1$
    !>                                                                       (given $N$ ODEs)
    !>        \param[in] "real(sp/dp)                   ::  x1"              initial time.
    !>        \param[in] "real(sp/dp)                   ::  x2"              final time.
    !>        \param[in] "real(sp/dp)                   ::  h"               guessed first stepsize.
    !>        \param[in] "interface                     ::  derivs_sp/dp"     derivatives $dydx$ of $y$ at $x$.

    !     INTENT(INOUT)
    !         None

    !     INTENT(OUT)
    !>        \param[out] "real(sp/dp),  dimension(:), allocatable :: xout"  storage for outputs (time).
    !>        \param[out] "real(sp/dp),  dimension(:), allocatable :: yout"  storage for outputs
    !>                                                                       (incremented variables).
    !>                                                                       dim 1 = function evaluations
    !>                                                                               at any time point.
    !>                                                                       dim 2 = number of equations.

    !     INTENT(IN), OPTIONAL
    !>        \param[in] "real(sp/dp),  optional         ::  hmin"                  minimum allowed stepsize (can be 0.)
    !>                                                                              DEFAULT: 0.0
    !>        \param[in] "real(sp/dp),  optional         ::  eps"                   accuracy (overall tolerance level)
    !>                                                                              DEFAULT: 1E-6
    !>        \param[in] "real(sp/dp),  dimension(:)     ::  para"                  parameter for the derivatives $dydx$

    !     INTENT(INOUT), OPTIONAL
    !         None

    !     INTENT(OUT), OPTIONAL
    !         None

    !     RESTRICTIONS
    !>       \note The user has to supply the subroutine derivs(x,y,dydx), which returns derivatives $dydx$ at $x$.

    !     EXAMPLE
    !           call RK4as( ystart, x1, x2, h, derivs, xout, yout, hmin, eps )
    !           --> see example in test directory --> test/test_mo_ode_solver

    !     LITERATURE
    !        1) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 77 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 1 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1992
    !        2) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1996

    !     HISTORY
    !>        \author Giovanni Dalmasso
    !>        \date Jul 2012
    !         Modified, Giovanni Dalmasso, Mar 2013
    !         Modified, Sebastian Mueller, Jan 2015
    interface RK4as
        module procedure RK4as_sp, RK4as_dp, RK4as_para_sp, RK4as_para_dp
    end interface RK4as


    ! ------------------------------------------------------------------

    !     NAME
    !         RBstiff

    !     PURPOSE
    !         Integration of Ordinary Differential Equations using a Rosenbrock Method
    !         with adaptive stepsize control.
    !         Integrate the array of starting values ystart from x1 to x2 with accuracy eps, storing intermediate results
    !         in the module variables. h1 should be set as a guessed first stepsize,
    !         hmin as the minimum allowed stepsize (can be zero).
    !         On output ystart is replaced by values at the end of the integration interval.

    !>        \brief Rosenbrock Method with adaptive stepsize control.

    !>        \details Integrate the array of starting values $y_{start} from $x_1$ to $x_2$ with accuracy $\varepsilon$,
    !>                 storing intermediate results in the module variables.
    !>                 $h_1$ should be set as a guessed first stepsize, $h_{min} as the minimum allowed stepsize (can be zero).
    !>                 On output $y_{start}$ is replaced by values at the end of the integration interval.
    !>                 If you use "para" as parameters for derivs/jacobn, the interface for derivs has to look like:
    !>
    !>                            interface
    !>                                subroutine derivs( x, y, para, dydx )
    !>                                    use mo_kind, only: sp/dp
    !>                                    implicit none
    !>                                    real(sp/dp),                 intent(in)  :: x        ! time
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: para     ! parameter for the derivatives
    !>                                    real(sp/dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
    !>                                end subroutine derivs
    !>                                subroutine jacobn(x, y, para, dfdx, dfdy)
    !>                                    use mo_kind, only: sp
    !>                                    implicit none
    !>                                    real(sp/dp),                intent(in)   :: x        ! time
    !>                                    real(sp/dp),dimension(:),   intent(in)   :: y        ! unknowns of the equations
    !>                                    real(sp/dp),dimension(:),   intent(in)   :: para     ! parameters for derivs
    !>                                    real(sp/dp),dimension(:),   intent(out)  :: dfdx     ! derivatives of f for x
    !>                                    real(sp/dp),dimension(:,:), intent(out)  :: dfdy     ! jacobi-matrix df/dy
    !>                                end subroutine jacobn
    !>                            end interface
    !>
    !>                 Elsewise "para" must be left out:
    !>
    !>                            interface
    !>                                subroutine derivs( x, y, dydx )
    !>                                    use mo_kind, only: sp/dp
    !>                                    implicit none
    !>                                    real(sp/dp),                 intent(in)  :: x        ! time
    !>                                    real(sp/dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
    !>                                    real(sp/dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
    !>                                end subroutine derivs
    !>                                subroutine jacobn(x, y, dfdx, dfdy)
    !>                                    use mo_kind, only: sp
    !>                                    implicit none
    !>                                    real(sp/dp),                intent(in)   :: x        ! time
    !>                                    real(sp/dp),dimension(:),   intent(in)   :: y        ! unknowns of the equations
    !>                                    real(sp/dp),dimension(:),   intent(out)  :: dfdx     ! derivatives of f for x
    !>                                    real(sp/dp),dimension(:,:), intent(out)  :: dfdy     ! jacobi-matrix df/dy
    !>                                end subroutine jacobn
    !>                            end interface

    !     INTENT(IN)
    !>        \param[in] "real(sp/dp),  dimension(:)    ::  vstart"          initial conditions.
    !>                                                                       $N$ inital values known at the time $x_1$
    !>                                                                       (given $N$ ODEs)
    !>        \param[in] "real(sp/dp)                   ::  x1"              initial time.
    !>        \param[in] "real(sp/dp)                   ::  x2"              final time.
    !>        \param[in] "real(sp/dp)                   ::  h"               guessed first stepsize.
    !>        \param[in] "interface                     ::  derivs_sp/dp"    derivatives $dydx$ of $y$ at $x$.
    !>        \param[in] "interface                     ::  jacobn_sp/dp"    derivatives of the RHS: $dfdx$ and $dfdy$

    !     INTENT(INOUT)
    !         None

    !     INTENT(OUT)
    !>        \param[out] "real(sp/dp),  dimension(:), allocatable :: xout"  storage for outputs (time).
    !>        \param[out] "real(sp/dp),  dimension(:), allocatable :: yout"  storage for outputs
    !>                                                                       (incremented variables).
    !>                                                                       dim 1 = function evaluations
    !>                                                                               at any time point.
    !>                                                                       dim 2 = number of equations.

    !     INTENT(IN), OPTIONAL
    !>        \param[in] "real(sp/dp),  optional         ::  hmin"                  minimum allowed stepsize (can be 0.)
    !>                                                                              DEFAULT: 0.0
    !>        \param[in] "real(sp/dp),  optional         ::  eps"                   accuracy (overall tolerance level)
    !>                                                                              DEFAULT: 1E-6
    !>        \param[in] "real(sp/dp),  dimension(:)     ::  para"                  parameter for the derivatives

    !     INTENT(INOUT), OPTIONAL
    !         None

    !     INTENT(OUT), OPTIONAL
    !         None

    !     RESTRICTIONS
    !>       \note The user has to supply the subroutine derivs(x,y,dydx), which returns derivatives $dydx$ at $x$.

    !     EXAMPLE
    !           call RK4as( ystart, x1, x2, h, derivs, jacobn, xout, yout, hmin, eps )
    !           --> see example in test directory --> test/test_mo_ode_solver

    !     LITERATURE
    !        1) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 77 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 1 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1992
    !        2) Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
    !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
    !             Cambridge University Press, UK, 1996

    !     HISTORY
    !>        \author Sebastian Mueller
    !>        \date Mar 2015

    interface RBstiff
        module procedure RBstiff_sp, RBstiff_dp, RBstiff_para_sp, RBstiff_para_dp
    end interface RBstiff

    private

    ! Private method
    interface CashKarpRK
        module procedure CashKarpRK_sp, CashKarpRK_dp, CashKarpRK_para_sp, CashKarpRK_para_dp
    end interface CashKarpRK

    interface RBstep
        module procedure RBstep_sp, RBstep_dp, RBstep_para_sp, RBstep_para_dp
    end interface RBstep

    ! ------------------------------------------------------------------

contains

    ! ------------------------------------------------------------------

    ! SINGLE PRECISION Euler
    subroutine Euler_sp( ystart, x1, x2, h, derivs, xout, yout )    ! all obligatory

        use mo_kind, only: i4, sp

        implicit none

        ! Intent IN
        real(sp),   dimension(:),                intent(in)  :: ystart   ! initial conditions
        real(sp),                                intent(in)  :: x1, x2   ! initial and final time
        real(sp),                                intent(in)  :: h        ! step size

        ! Intent OUT
        real(sp),   dimension(:),   allocatable, intent(out) :: xout
        real(sp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)  :: x        ! time
                real(sp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(sp),   dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                           :: j        ! counter
        integer(i4)                           :: nstep    ! number of steps
        real(sp)                              :: x
        real(sp),   dimension( size(ystart) ) :: dy, y

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout) ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )            ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        do j=1, nstep                   ! take nstep steps
            call derivs( x, y, dy )
            y = y + h*dy                ! step EULER
            if ( h .lt. tiny(1._sp) )    stop 'Euler_sp --> stepsize not significant!'
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine Euler_sp


    ! DOUBLE PRECISION Euler
    subroutine Euler_dp( ystart, x1, x2, h, derivs, xout, yout )    ! all obligatory

        use mo_kind, only: i4, dp

        implicit none

        ! Intent IN
        real(dp),   dimension(:),                intent(in)  :: ystart   ! initial conditions
        real(dp),                                intent(in)  :: x1, x2   ! initial and final time
        real(dp),                                intent(in)  :: h        ! step size

        ! Intent OUT
        real(dp),   dimension(:),   allocatable, intent(out) :: xout
        real(dp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)  :: x        ! time
                real(dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                           :: j        ! counter
        integer(i4)                           :: nstep    ! nuber of steps
        real(dp)                              :: x
        real(dp),   dimension( size(ystart) ) :: dy, y

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout) ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )            ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        do j=1, nstep                   ! take nstep steps
            call derivs( x, y, dy )
            y = y + h*dy                ! step EULER
            if ( h .lt. tiny(1._dp) )    stop 'Euler_dp --> stepsize not significant!'
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine Euler_dp

    ! ------------------------------------------------------------------

    ! SINGLE PRECISION 4th order RUNGE-KUTTA
    subroutine RK4_sp( ystart, x1, x2, h, derivs, xout, yout )    ! all obligatory

        use mo_kind, only: i4, sp

        implicit none

        ! Intent IN
        real(sp),   dimension(:),                intent(in)  :: ystart   ! initial conditions
        real(sp),                                intent(in)  :: x1, x2   ! initial and final time
        real(sp),                                intent(in)  :: h        ! step size

        ! Intent OUT
        real(sp),   dimension(:),   allocatable, intent(out) :: xout
        real(sp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)  :: x        ! time
                real(sp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(sp),   dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                           :: j            ! counter
        integer(i4)                           :: nstep        ! nuber of steps
        real(sp)                              :: x
        real(sp)                              :: hh, h6, xh
        real(sp),   dimension( size(ystart) ) :: dy, dyt, dym, y, yt

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout) ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )            ! allocate storage for saved values
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
            if ( h .lt. tiny(1._sp) )    stop 'RK4_sp --> stepsize not significant!'
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine RK4_sp


    ! DOUBLE PRECISION 4th order RUNGE-KUTTA
    subroutine RK4_dp( ystart, x1, x2, h, derivs, xout, yout )    ! all obligatory

        use mo_kind, only: i4, dp

        implicit none

        ! Intent IN
        real(dp),   dimension(:),                intent(in)  :: ystart   ! initial conditions
        real(dp),                                intent(in)  :: x1, x2   ! initial and final time
        real(dp),                                intent(in)  :: h        ! step size

        ! Intent OUT
        real(dp),   dimension(:),   allocatable, intent(out) :: xout
        real(dp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)  :: x        ! time
                real(dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                           :: j            ! counter
        integer(i4)                           :: nstep        ! nuber of steps
        real(dp)                              :: x
        real(dp)                              :: hh, h6, xh
        real(dp),   dimension( size(ystart) ) :: dy, dyt, dym, y, yt

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout) ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )            ! allocate storage for saved values
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
            if ( h .lt. tiny(1._dp) )    stop 'RK4_sp --> stepsize not significant!'
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine RK4_dp

    ! ------------------------------------------------------------------

    ! SINGLE PRECISION 4th order RUNGE-KUTTA Adaptive Step-size
    subroutine RK4as_sp( ystart, x1, x2, h, derivs, xout, yout, &   ! obligatory
        hmin, eps )                                                 ! optional

        use mo_kind,   only: i4, sp
        use mo_nrutil, only: reallocate

        implicit none

        ! Intent IN
        real(sp),                                intent(in) :: x1, x2    ! initial and final time
        real(sp),                                intent(in) :: h         ! guessed step size
        real(sp),   dimension(:),                intent(in) :: ystart    ! initial conditions
        real(sp),                   optional,    intent(in) :: hmin
        real(sp),                   optional,    intent(in) :: eps

        ! Intent OUT
        real(sp),   dimension(:),   allocatable, intent(out) :: xout
        real(sp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)  :: x        ! time
                real(sp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(sp),   dimension(:), intent(out) :: dydx     ! derivatives of y
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
        integer(i4),    parameter   :: MAXstp = 1000000_i4  ! max number of steps
        real(sp),       parameter   :: safety = 0.9_sp, pgrow = -0.2_sp, pshrnk = -0.25_sp
        real(sp) :: errcon

        errcon = exp( (1._sp/pgrow)*log(5._sp/safety) )

        if( present(hmin) ) then
            hminIN = hmin
        else
            hminIN = 0.0_sp
        end if

        if( present(eps) ) then
            epsIN = eps
        else
            epsIN = 1e-6_sp
        end if

        x = x1
        hIN = sign( h, x2-x1 )
        nok = 0_i4
        nbad = 0_i4
        kount = 0_i4
        y(:) = ystart(:)
        dxsav = tiny(1._sp)

        xsav = x-2.0_sp*dxsav                                                   ! assures storage of first step
        nullify( xp, yp )
        allocate( xp(256) )
        allocate( yp(size(xp), size(ystart)) )

        call save_a_step                                                        ! save initial step

        do nstep=1, MAXstp                                                      ! take at most MAXstp steps

            call derivs( x, y, dydx )
            yscal(:) = abs( y(:) ) + abs( hIN*dydx(:) ) + tiny(1._sp)           ! scaling used to monitor accuracy
                                                                                ! --> CAN BE MODIFIED...

            if ( abs(x-xsav) .gt. abs(dxsav) ) call save_a_step                 ! store intermediate results

            if ( (x+hIN-x2)*(x+hIN-x1) .gt. 0.0_sp ) hIN = x2-x                 ! if stepsize can overshoot, decrease

            do
                call CashKarpRK( y, dydx, x, hIN, ytemp, yerr, derivs )         ! take a step
                errmax = maxval( abs(yerr(:)/yscal(:)) )/epsIN                  ! evaluate accuracy
                if ( errmax .lt. 1.0_sp )   exit                                ! step succeeded
                htemp = safety*hIN*(errmax**pshrnk)                             ! truncation error too large, reduce stepsize
                hIN = sign( max( abs(htemp), 0.1_sp*abs(hIN) ), hIN )           ! no more than a factor of 10
                xnew = x+hIN
                if ( abs(xnew-x) .lt. epsilon(1.0_sp) )  stop 'RK4as_sp --> hey !!! stepsize underflow!'
            end do

            if ( errmax .gt. errcon )   then                                    ! compute size of next step
                hnext = safety*hIN*(errmax**pgrow)
            else                                                                ! no more than a factor of 5 increase
                hnext = 5.0_sp*hIN
            end if

            hdid = hIN
            x = x+hIN
            y(:) = ytemp(:)

            if ( abs(hdid-hIN) .gt. epsilon(1.0_sp) )    then
                nbad = nbad+1_i4
            else
                nok = nok+1_i4
            end if

            if ( (x-x2)*(x2-x1) .ge. 0.0_sp )   then            ! are we done?!?!
                call save_a_step                                ! save final step
                allocate( xout(kount) )                         ! allocate storage for outputs
                xout(:) = xp(1:kount)
                allocate( yout(kount, size(yp,2)) )             ! allocate storage for outputs
                yout(:,:) = yp(1:kount, :)
                deallocate( xp, yp )                            ! clear out old stored variables
                return                                          ! normal exit
            end if

            if ( abs(hnext-hminIN) .lt. epsilon(1.0_sp) )   stop 'RK4as_sp --> WTF! ...stepsize smaller than minimum!!!'
            hIN = hnext

        end do

        stop 'RK4as_sp --> dude, too many steps!!!'

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

        use mo_kind,   only: i4, dp
        use mo_nrutil, only: reallocate

        implicit none

        ! Intent IN
        real(dp),                                intent(in) :: x1, x2    ! initial and final time
        real(dp),                                intent(in) :: h         ! guessed step size
        real(dp),   dimension(:),                intent(in) :: ystart    ! initial conditions
        real(dp),                   optional,    intent(in) :: hmin
        real(dp),                   optional,    intent(in) :: eps

        ! Intent OUT
        real(dp),   dimension(:),   allocatable, intent(out) :: xout
        real(dp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)  :: x        ! time
                real(dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
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
        real(dp),       parameter   :: safety = 0.9_dp, pgrow = -0.2_dp, pshrnk = -0.25_dp
        real(dp) :: errcon

        errcon = exp( (1._dp/pgrow)*log(5._dp/safety) )

        if( present(hmin) ) then
            hminIN = hmin
        else
            hminIN = 0.0_dp
        end if

        if( present(eps) ) then
            epsIN = eps
        else
            epsIN = 1e-6_dp
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
                if ( errmax .lt. 1.0_dp )   exit                        ! step succeeded
                htemp = safety*hIN*(errmax**pshrnk)                     ! truncation error too large, reduce stepsize
                hIN = sign( max( abs(htemp), 0.1_dp*abs(hIN) ), hIN )   ! no more than a factor of 10
                xnew = x+hIN
                if ( abs(xnew-x) .lt. epsilon(1._dp) )  stop 'RK4as_dp --> hey!!! stepsize underflow!'
            end do

            if ( errmax .gt. errcon )   then                            ! compute size of next step
                hnext = safety*hIN*(errmax**pgrow)
            else                                                        ! no more than a factor of 5 increase
                hnext = 5.0_dp*hIN
            end if

            hdid = hIN
            x = x+hIN
            y(:) = ytemp(:)

            if ( abs(hdid-hIN) .gt. epsilon(1._dp) ) then
                nbad = nbad+1_i4
            else
                nok = nok+1_i4
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

            if ( abs(hnext-hminIN) .lt. epsilon(1.0_dp) )   stop 'RK4as_dp --> WTF! ...stepsize smaller than minimum!!!'
            hIN = hnext

        end do

        stop 'RK4as_dp --> dude, too many steps!!!'

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

    ! ------------------------------------------------------------------

    ! SINGLE PRECISION Rosenbrock Method for stiff ode's
    subroutine RBstiff_sp( ystart, x1, x2, h, derivs, jacobn, xout, yout, &                 ! obligatory
        hmin, eps )                                                                         ! optional

        use mo_kind,            only : i4, sp
        use mo_utils,           only : ge
        use mo_nrutil,          only : reallocate

        implicit none

        ! Intent IN
        real(sp),                                intent(in)     :: x1, x2    ! initial and final time
        real(sp),                                intent(in)     :: h         ! guessed step size
        real(sp),   dimension(:),                intent(in)     :: ystart    ! initial conditions
        real(sp),                   optional,    intent(in)     :: hmin
        real(sp),                   optional,    intent(in)     :: eps

        ! Intent OUT
        real(sp),   dimension(:),   allocatable, intent(out)    :: xout
        real(sp),   dimension(:,:), allocatable, intent(out)    :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)            :: x        ! time
                real(sp),   dimension(:), intent(in)            :: y        ! unknowns of the equations
                real(sp),   dimension(:), intent(out)           :: dydx     ! derivatives of y
            end subroutine derivs

            subroutine jacobn(x, y, dfdx, dfdy)
                use mo_kind, only: sp
                implicit none
                real(sp),                   intent(in)          :: x
                real(sp),   dimension(:),   intent(in)          :: y
                real(sp),   dimension(:),   intent(out)         :: dfdx
                real(sp),   dimension(:,:), intent(out)         :: dfdy
            end subroutine jacobn
        end interface

        ! Internal variables
        integer(i4)                         :: n    
        integer(i4)                         :: nstep                        ! nuber of steps
        integer(i4)                         :: kount                        ! total number of saved steps
        real(sp)                            :: hIN, hdid, hnext, x, hminIN, epsIN, xsav, dxsav
        real(sp),   dimension(size(ystart)) :: dydx, y, yscal, C

        ! Pointers
        real(sp),   dimension(:),   pointer :: xp
        real(sp),   dimension(:,:), pointer :: yp

        ! parameters
        integer(i4),    parameter           :: MAXstp = 1000000_i4          ! max number of steps


        if( present(hmin) ) then
            hminIN = abs(hmin)
        else
            hminIN = 0.0_sp
        end if

        if( present(eps) ) then
            epsIN = abs(eps)
        else
            epsIN = 1e-6_sp
        end if
        
        C       = 1.0_sp                                                !lower boundary for scaling

        hIN     = sign( h, x2-x1 )
        hnext   = 0.0_sp
        hdid    = 0.0_sp

        kount   = 0_i4

        x       = x1
        y(:)    = ystart(:)
        dxsav   = tiny(1._sp)
        xsav    = x-2.0_sp*dxsav                                        ! assures storage of first step

        nullify( xp, yp )
        allocate( xp(256) )
        allocate( yp(size(xp), size(ystart)) )

        call save_a_step                                                ! save initial step

        do nstep=1, MAXstp                                              ! take at most MAXstp steps

            call derivs( x, y, dydx )
            forall( n = 1:size(yscal) ) yscal(n) = max (abs(y(n)),C(n)) ! scaling used to monitor accuracy --> CAN BE MODIFIED...

            if ( abs(x-xsav) .gt. abs(dxsav) ) call save_a_step         ! store intermediate results

            if ( (x+hIN-x2)*(x+hIN-x1) .gt. 0.0_sp ) hIN = x2-x         ! if stepsize can overshoot, decrease

            call RBstep( y, dydx, x, hIN, epsIN, yscal, hdid, hnext, derivs, jacobn) !do one Rosenbrock-Step

            if ( ge((x-x2)*(x2-x1), 0.0_sp) )   then                    ! are we done?!?!
                call save_a_step                                        ! save final step
                allocate( xout(kount) )                                 ! allocate storage for outputs
                xout(:) = xp(1:kount)
                allocate( yout(kount, size(yp,2)) )                     ! allocate storage for outputs
                yout(:,:) = yp(1:kount, :)
                deallocate( xp, yp )                                    ! clear out old stored variables
                return                                                  ! normal exit
            end if

            if ( abs(hnext-hminIN) .lt. epsilon(1.0_sp) )   stop 'RBstiff_sp --> WTF! ...stepsize smaller than minimum!!!'
            
            hIN     = hnext

        end do

        stop 'RBstiff_sp --> dude, too many steps!!!'

    contains

        subroutine save_a_step      ! --> like a macro in C
            kount = kount+1_i4
            if ( kount .gt. size(xp) ) then
                xp=>reallocate( xp, 2*size(xp) )
                yp=>reallocate( yp, size(xp), size(yp,2) )
            end if
            xp(kount) = x
            yp(kount, :) = y(:)
            xsav = x
        end subroutine save_a_step

    end subroutine RBstiff_sp

    ! DOUBLE PRECISION Rosenbrock Method for stiff ode's
    subroutine RBstiff_dp( ystart, x1, x2, h, derivs, jacobn, xout, yout, &                 ! obligatory
        hmin, eps )                                                                         ! optional

        use mo_kind,            only : i4, dp
        use mo_utils,           only : ge
        use mo_nrutil,          only : reallocate

        implicit none

        ! Intent IN
        real(dp),                                intent(in)     :: x1, x2    ! initial and final time
        real(dp),                                intent(in)     :: h         ! guessed step size
        real(dp),   dimension(:),                intent(in)     :: ystart    ! initial conditions
        real(dp),                   optional,    intent(in)     :: hmin
        real(dp),                   optional,    intent(in)     :: eps

        ! Intent OUT
        real(dp),   dimension(:),   allocatable, intent(out)    :: xout
        real(dp),   dimension(:,:), allocatable, intent(out)    :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)            :: x        ! time
                real(dp),   dimension(:), intent(in)            :: y        ! unknowns of the equations
                real(dp),   dimension(:), intent(out)           :: dydx     ! derivatives of y
            end subroutine derivs

            subroutine jacobn(x, y, dfdx, dfdy)
                use mo_kind, only: dp
                implicit none
                real(dp),                   intent(in)          :: x
                real(dp),   dimension(:),   intent(in)          :: y
                real(dp),   dimension(:),   intent(out)         :: dfdx
                real(dp),   dimension(:,:), intent(out)         :: dfdy
            end subroutine jacobn
        end interface

        ! Internal variables
        integer(i4)                         :: n    
        integer(i4)                         :: nstep                        ! nuber of steps
        integer(i4)                         :: kount                        ! total number of saved steps
        real(dp)                            :: hIN, hdid, hnext, x, hminIN, epsIN, xsav, dxsav
        real(dp),   dimension(size(ystart)) :: dydx, y, yscal, C

        ! Pointers
        real(dp),   dimension(:),   pointer :: xp
        real(dp),   dimension(:,:), pointer :: yp

        ! parameters
        integer(i4),    parameter           :: MAXstp = 1000000_i4          ! max number of steps


        if( present(hmin) ) then
            hminIN = abs(hmin)
        else
            hminIN = 0.0_dp
        end if

        if( present(eps) ) then
            epsIN = abs(eps)
        else
            epsIN = 1e-6_dp
        end if
        
        C       = 1.0_dp                                                !lower boundary for scaling

        hIN     = sign( h, x2-x1 )
        hnext   = 0.0_dp
        hdid    = 0.0_dp

        kount   = 0_i4

        x       = x1
        y(:)    = ystart(:)
        dxsav   = tiny(1._dp)
        xsav    = x-2.0_dp*dxsav                                        ! assures storage of first step

        nullify( xp, yp )
        allocate( xp(256) )
        allocate( yp(size(xp), size(ystart)) )

        call save_a_step                                                ! save initial step

        do nstep=1, MAXstp                                              ! take at most MAXstp steps

            call derivs( x, y, dydx )
            forall( n = 1:size(yscal) ) yscal(n) = max (abs(y(n)),C(n)) ! scaling used to monitor accuracy --> CAN BE MODIFIED...

            if ( abs(x-xsav) .gt. abs(dxsav) ) call save_a_step         ! store intermediate results

            if ( (x+hIN-x2)*(x+hIN-x1) .gt. 0.0_dp ) hIN = x2-x         ! if stepsize can overshoot, decrease

            call RBstep( y, dydx, x, hIN, epsIN, yscal, hdid, hnext, derivs, jacobn) !do one Rosenbrock-Step

            if ( ge((x-x2)*(x2-x1), 0.0_dp) )   then                    ! are we done?!?!
                call save_a_step                                        ! save final step
                allocate( xout(kount) )                                 ! allocate storage for outputs
                xout(:) = xp(1:kount)
                allocate( yout(kount, size(yp,2)) )                     ! allocate storage for outputs
                yout(:,:) = yp(1:kount, :)
                deallocate( xp, yp )                                    ! clear out old stored variables
                return                                                  ! normal exit
            end if

            if ( abs(hnext-hminIN) .lt. epsilon(1.0_dp) )   stop 'RBstiff_dp --> WTF! ...stepsize smaller than minimum!!!'
            
            hIN     = hnext

        end do

        stop 'RBstiff_dp --> dude, too many steps!!!'

    contains

        subroutine save_a_step      ! --> like a macro in C
            kount = kount+1_i4
            if ( kount .gt. size(xp) ) then
                xp=>reallocate( xp, 2*size(xp) )
                yp=>reallocate( yp, size(xp), size(yp,2) )
            end if
            xp(kount) = x
            yp(kount, :) = y(:)
            xsav = x
        end subroutine save_a_step

    end subroutine RBstiff_dp

    ! ------------------------------------------------------------------
    !   PRIVATE METHODS
    ! ------------------------------------------------------------------

    ! SINGLE PRECISION CASH-KARP RUNGE-KUTTA step
    subroutine CashKarpRK_sp( y, dydx, x, h, yout, yerr, derivs )
        ! Given values for N variables y and their derivatives dydx known at x, use the fifth-order
        ! Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
        ! the incremented variables as yout. Also return an estimate of the local truncation error
        ! in yout using the embedded fourth order method.

        use mo_kind,   only : i4, sp
        use mo_nrutil, only : assert_eq

        implicit none

        ! Intent IN
        real(sp),   dimension(:), intent(in)  :: y, dydx
        real(sp),                 intent(in)  :: x, h

        ! Intent OUT
        real(sp),   dimension(:), intent(out) :: yout, yerr

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)  :: x
                real(sp),   dimension(:), intent(in)  :: y
                real(sp),   dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: ndum
        real(sp),   dimension(:),   allocatable :: ak2, ak3, ak4, ak5, ak6, ytemp

        ! parameters
        real(sp),   parameter   :: A2=0.2_sp, A3=0.3_sp, A4=0.6_sp, A5=1._sp, A6=0.875_sp, B21=0.2_sp, B31=3._sp/40._sp, &
            B32=9._sp/40._sp, B41=0.3_sp, B42=-0.9_sp, B43=1.2_sp, B51=-11._sp/54._sp, B52=2.5_sp, &
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

        use mo_kind,   only: i4, dp
        use mo_nrutil, only: assert_eq

        implicit none

        ! Intent IN
        real(dp),   dimension(:), intent(in)  :: y, dydx
        real(dp),                 intent(in)  :: x, h

        ! Intent OUT
        real(dp),   dimension(:), intent(out) :: yout, yerr

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)  :: x
                real(dp),   dimension(:), intent(in)  :: y
                real(dp),   dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: ndum
        real(dp),   dimension(:),   allocatable :: ak2, ak3, ak4, ak5, ak6, ytemp

        ! parameters
        real(dp),   parameter   :: A2=0.2_dp, A3=0.3_dp, A4=0.6_dp, A5=1._dp, A6=0.875_dp, B21=0.2_dp, B31=3._dp/40._dp, &
            B32=9._dp/40._dp, B41=0.3_dp, B42=-0.9_dp, B43=1.2_dp, B51=-11._dp/54._dp, B52=2.5_dp, &
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

    ! ------------------------------------------------------------------

    ! SINGLE PRECISION ROSENBROCK step
    subroutine RBstep_sp(y, dydx, x, htry, eps, yscal, hdid, hnext, derivs, jacobn)

        use mo_kind,            only : i4, sp
        use mo_utils,           only : eq, le
        use mo_nrutil,          only : assert_eq, diagadd
        use mo_linear_algebra,  only : solve_linear_equations

        implicit none

        real(sp),       dimension(:),   intent(inout)   :: y
        real(sp),       dimension(:),   intent(in)      :: dydx, yscal
        real(sp),                       intent(inout)   :: x

        real(sp),                       intent(in)      :: htry, eps
        real(sp),                       intent(out)     :: hdid, hnext

        interface
            subroutine derivs(x, y, dydx)
                use mo_kind, only: sp
                implicit none
                real(sp),                   intent(in)  :: x
                real(sp),   dimension(:),   intent(in)  :: y
                real(sp),   dimension(:),   intent(out) :: dydx
            end subroutine derivs

            subroutine jacobn(x, y, dfdx, dfdy)
                use mo_kind, only: sp
                implicit none
                real(sp),                   intent(in)  :: x
                real(sp),   dimension(:),   intent(in)  :: y
                real(sp),   dimension(:),   intent(out) :: dfdx
                real(sp),   dimension(:,:), intent(out) :: dfdy
            end subroutine jacobn
        end interface

        integer(i4)                                     :: jtry,ndum
!        integer(i4),    dimension(size(y))              :: indx
        real(sp),       dimension(size(y))              :: dfdx,dytmp,err,g1,g2,g3,g4,ysav
        real(sp),       dimension(size(y),size(y))      :: a,dfdy
        real(sp)                                        :: errmax,h,xsav !,d

        integer(i4), parameter  ::  maxtry  = 40_i4
        real(sp), parameter     ::  safety  = 0.9_sp,               grow    = 1.5_sp,           pgrow   = -0.25_sp,&
                                    shrnk   = 0.5_sp,               pshrnk  = -1.0_sp/3.0_sp,   errcon  = 0.1296_sp,&
                                    gam     = 1.0_sp/2.0_sp,&
                                    a21     = 2.0_sp,               a31     = 48.0_sp/25.0_sp,  a32     = 6.0_sp/25.0_sp,&
                                    c21     = -8.0_sp,              c31     = 372.0_sp/25.0_sp, c32     = 12.0_sp/5.0_sp,&
                                    c41     = -112.0_sp/125.0_sp,   c42     = -54._sp/125.0_sp, c43     = -2.0_sp/5.0_sp,&
                                    b1      = 19.0_sp/9.0_sp,       b2      = 1.0_sp/2.0_sp,    b3      = 25.0_sp/108.0_sp,&
                                    b4      = 125.0_sp/108.0_sp,&
                                    e1      = 17.0_sp/54.0_sp,      e2      = 7.0_sp/36.0_sp,   e3      = 0.0_sp,&
                                    e4      = 125.0_sp/108.0_sp,&
                                    c1x     = 1.0_sp/2.0_sp,        c2x     = -3.0_sp/2.0_sp,   c3x     = 121.0_sp/50.0_sp,&
                                    c4x     = 29.0_sp/250.0_sp,&
                                    a2x     = 1.0_sp,               a3x     = 3.0_sp/5.0_sp


        ndum        =   assert_eq(size(y),size(dydx),size(yscal),'stiff')
        xsav        =   x
        ysav(:)     =   y(:)

        call jacobn(xsav,ysav,dfdx,dfdy)

        h           =   htry

        do jtry=1,maxtry

            a(:,:)  =   -dfdy(:,:)

            call diagadd(a,1.0_sp/(gam*h))
!            call ludcmp(a,indx,d)

            g1      =   dydx+h*c1x*dfdx
            g1      =   solve_linear_equations(a,g1)

 !           call lubksb(a,indx,g1)

            y       =   ysav+a21*g1
            x       =   xsav+a2x*h

            call derivs(x,y,dytmp)

            g2      =   dytmp+h*c2x*dfdx+c21*g1/h
            g2      =   solve_linear_equations(a,g2)

  !          call lubksb(a,indx,g2)

            y       =   ysav+a31*g1+a32*g2
            x       =   xsav+a3x*h

            call derivs(x,y,dytmp)

            g3      =   dytmp+h*c3x*dfdx+(c31*g1+c32*g2)/h
            g3      =   solve_linear_equations(a,g3)

   !         call lubksb(a,indx,g3)

            g4      =   dytmp+h*c4x*dfdx+(c41*g1+c42*g2+c43*g3)/h
            g4      =   solve_linear_equations(a,g4)

    !        call lubksb(a,indx,g4)

            y       =   ysav+b1*g1+b2*g2+b3*g3+b4*g4
            err     =   e1*g1+e2*g2+e3*g3+e4*g4
            x       =   xsav+h

            if (eq(x,xsav)) stop 'stepsize not significant in RBstep_sp'

            errmax  =   maxval(abs(err/yscal))/eps

            if (le(errmax,1.0_sp)) then
                hdid    =   h
                hnext   =  merge(safety*h*errmax**pgrow, grow*h, errmax > errcon)
                return
            else
                hnext   =   safety*h*errmax**pshrnk
                h       =   sign(max(abs(hnext),shrnk*abs(h)),h)
            end if

        end do

        stop 'exceeded maxtry in RBstep_sp'

    end subroutine RBstep_sp


    ! DOUBLE PRECISION ROSENBROCK step
    subroutine RBstep_dp(y, dydx, x, htry, eps, yscal, hdid, hnext, derivs, jacobn)

        use mo_kind,            only : i4, dp
        use mo_utils,           only : eq, le
        use mo_nrutil,          only : assert_eq, diagadd
        use mo_linear_algebra,  only : solve_linear_equations

        implicit none

        real(dp),       dimension(:),   intent(inout)   :: y
        real(dp),       dimension(:),   intent(in)      :: dydx, yscal
        real(dp),                       intent(inout)   :: x

        real(dp),                       intent(in)      :: htry, eps
        real(dp),                       intent(out)     :: hdid, hnext

        interface
            subroutine derivs(x, y, dydx)
                use mo_kind, only: dp
                implicit none
                real(dp),                   intent(in)  :: x
                real(dp),   dimension(:),   intent(in)  :: y
                real(dp),   dimension(:),   intent(out) :: dydx
            end subroutine derivs

            subroutine jacobn(x, y, dfdx, dfdy)
                use mo_kind, only: dp
                implicit none
                real(dp),                   intent(in)  :: x
                real(dp),   dimension(:),   intent(in)  :: y
                real(dp),   dimension(:),   intent(out) :: dfdx
                real(dp),   dimension(:,:), intent(out) :: dfdy
            end subroutine jacobn
        end interface

        integer(i4)                                     :: jtry,ndum
!        integer(i4),    dimension(size(y))              :: indx
        real(dp),       dimension(size(y))              :: dfdx,dytmp,err,g1,g2,g3,g4,ysav
        real(dp),       dimension(size(y),size(y))      :: a,dfdy
        real(dp)                                        :: errmax,h,xsav !,d

        integer(i4), parameter  ::  maxtry  = 40_i4
        real(dp), parameter     ::  safety  = 0.9_dp,               grow    = 1.5_dp,           pgrow   = -0.25_dp,&
                                    shrnk   = 0.5_dp,               pshrnk  = -1.0_dp/3.0_dp,   errcon  = 0.1296_dp,&
                                    gam     = 1.0_dp/2.0_dp,&
                                    a21     = 2.0_dp,               a31     = 48.0_dp/25.0_dp,  a32     = 6.0_dp/25.0_dp,&
                                    c21     = -8.0_dp,              c31     = 372.0_dp/25.0_dp, c32     = 12.0_dp/5.0_dp,&
                                    c41     = -112.0_dp/125.0_dp,   c42     = -54._dp/125.0_dp, c43     = -2.0_dp/5.0_dp,&
                                    b1      = 19.0_dp/9.0_dp,       b2      = 1.0_dp/2.0_dp,    b3      = 25.0_dp/108.0_dp,&
                                    b4      = 125.0_dp/108.0_dp,&
                                    e1      = 17.0_dp/54.0_dp,      e2      = 7.0_dp/36.0_dp,   e3      = 0.0_dp,&
                                    e4      = 125.0_dp/108.0_dp,&
                                    c1x     = 1.0_dp/2.0_dp,        c2x     = -3.0_dp/2.0_dp,   c3x     = 121.0_dp/50.0_dp,&
                                    c4x     = 29.0_dp/250.0_dp,&
                                    a2x     = 1.0_dp,               a3x     = 3.0_dp/5.0_dp


        ndum        =   assert_eq(size(y),size(dydx),size(yscal),'stiff')
        xsav        =   x
        ysav(:)     =   y(:)

        call jacobn(xsav,ysav,dfdx,dfdy)

        h           =   htry

        do jtry=1,maxtry

            a(:,:)  =   -dfdy(:,:)

            call diagadd(a,1.0_dp/(gam*h))
!            call ludcmp(a,indx,d)

            g1      =   dydx+h*c1x*dfdx
            g1      =   solve_linear_equations(a,g1)

 !           call lubksb(a,indx,g1)

            y       =   ysav+a21*g1
            x       =   xsav+a2x*h

            call derivs(x,y,dytmp)

            g2      =   dytmp+h*c2x*dfdx+c21*g1/h
            g2      =   solve_linear_equations(a,g2)

  !          call lubksb(a,indx,g2)

            y       =   ysav+a31*g1+a32*g2
            x       =   xsav+a3x*h

            call derivs(x,y,dytmp)

            g3      =   dytmp+h*c3x*dfdx+(c31*g1+c32*g2)/h
            g3      =   solve_linear_equations(a,g3)

   !         call lubksb(a,indx,g3)

            g4      =   dytmp+h*c4x*dfdx+(c41*g1+c42*g2+c43*g3)/h
            g4      =   solve_linear_equations(a,g4)

    !        call lubksb(a,indx,g4)

            y       =   ysav+b1*g1+b2*g2+b3*g3+b4*g4
            err     =   e1*g1+e2*g2+e3*g3+e4*g4
            x       =   xsav+h

            if (eq(x,xsav)) stop 'stepsize not significant in RBstep_dp'

            errmax  =   maxval(abs(err/yscal))/eps

            if (le(errmax,1.0_dp)) then
                hdid    =   h
                hnext   =  merge(safety*h*errmax**pgrow, grow*h, errmax > errcon)
                return
            else
                hnext   =   safety*h*errmax**pshrnk
                h       =   sign(max(abs(hnext),shrnk*abs(h)),h)
            end if

        end do

        stop 'exceeded maxtry in RBstep_dp'

    end subroutine RBstep_dp

    ! ------------------------------------------------------------------


!-----------------------------------------------------------|
!the second version with parameter-input for the derivatives|
!-----------------------------------------------------------|


    ! ------------------------------------------------------------------

    ! SINGLE PRECISION Euler
    subroutine Euler_para_sp( ystart, x1, x2, h, derivs, para, xout, yout )    ! all obligatory

        use mo_kind, only: i4, sp

        implicit none

        ! Intent IN
        real(sp),   dimension(:),                intent(in)  :: ystart   ! initial conditions
        real(sp),                                intent(in)  :: x1, x2   ! initial and final time
        real(sp),                                intent(in)  :: h        ! step size
        real(sp),   dimension(:),                intent(in)  :: para     ! parameters for derivs

        ! Intent OUT
        real(sp),   dimension(:),   allocatable, intent(out) :: xout
        real(sp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)  :: x        ! time
                real(sp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(sp),   dimension(:), intent(in)  :: para     ! parameters for derivs
                real(sp),   dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                           :: j        ! counter
        integer(i4)                           :: nstep    ! number of steps
        real(sp)                              :: x
        real(sp),   dimension( size(ystart) ) :: dy, y

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout) ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )            ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        do j=1, nstep                   ! take nstep steps
            call derivs( x, y, para, dy )
            y = y + h*dy                ! step EULER
            if ( h .lt. tiny(1._sp) )    stop 'Euler_para_sp --> stepsize not significant!'
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine Euler_para_sp

    ! DOUBLE PRECISION Euler
    subroutine Euler_para_dp( ystart, x1, x2, h, derivs, para, xout, yout )    ! all obligatory

        use mo_kind, only: i4, dp

        implicit none

        ! Intent IN
        real(dp),   dimension(:),                intent(in)  :: ystart   ! initial conditions
        real(dp),                                intent(in)  :: x1, x2   ! initial and final time
        real(dp),                                intent(in)  :: h        ! step size
        real(dp),   dimension(:),                intent(in)  :: para     ! parameters for derivs

        ! Intent OUT
        real(dp),   dimension(:),   allocatable, intent(out) :: xout
        real(dp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)  :: x        ! time
                real(dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(dp),   dimension(:), intent(in)  :: para     ! parameters for derivs
                real(dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                           :: j        ! counter
        integer(i4)                           :: nstep    ! nuber of steps
        real(dp)                              :: x
        real(dp),   dimension( size(ystart) ) :: dy, y

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout) ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )            ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        do j=1, nstep                   ! take nstep steps
            call derivs( x, y, para, dy )
            y = y + h*dy                ! step EULER
            if ( h .lt. tiny(1._dp) )    stop 'Euler_para_dp --> stepsize not significant!'
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine Euler_para_dp

    ! ------------------------------------------------------------------

    ! SINGLE PRECISION 4th order RUNGE-KUTTA
    subroutine RK4_para_sp( ystart, x1, x2, h, derivs, para, xout, yout )    ! all obligatory

        use mo_kind, only: i4, sp

        implicit none

        ! Intent IN
        real(sp),   dimension(:),                intent(in)  :: ystart   ! initial conditions
        real(sp),                                intent(in)  :: x1, x2   ! initial and final time
        real(sp),                                intent(in)  :: h        ! step size
        real(sp),   dimension(:),                intent(in)  :: para     ! parameters for derivs

        ! Intent OUT
        real(sp),   dimension(:),   allocatable, intent(out) :: xout
        real(sp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)  :: x        ! time
                real(sp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(sp),   dimension(:), intent(in)  :: para     ! parameters for derivs
                real(sp),   dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                           :: j            ! counter
        integer(i4)                           :: nstep        ! nuber of steps
        real(sp)                              :: x
        real(sp)                              :: hh, h6, xh
        real(sp),   dimension( size(ystart) ) :: dy, dyt, dym, y, yt

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout) ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )            ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        hh = h*0.5_sp
        h6 = h/6.0_sp
        xh = x+hh

        do j=1, nstep                          ! take nstep steps
            call derivs( x, y, para, dy )            ! first step
            yt = y + hh*dy
            call derivs( xh, yt, para, dyt )         ! second step
            yt = y + hh*dyt
            call derivs( xh, yt, para, dym )         ! third step
            yt = y + h*dym
            dym = dyt + dym
            call derivs( x+h, yt, para, dyt )        ! fourth step
            y = y + h6*( dy+dyt+2.0_sp*dym )   ! accumulate increments with proper weights
            if ( h .lt. tiny(1._sp) )    stop 'RK4_para_sp --> stepsize not significant!'
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine RK4_para_sp


    ! DOUBLE PRECISION 4th order RUNGE-KUTTA
    subroutine RK4_para_dp( ystart, x1, x2, h, derivs, para, xout, yout )    ! all obligatory

        use mo_kind, only: i4, dp

        implicit none

        ! Intent IN
        real(dp),   dimension(:),                intent(in)  :: ystart   ! initial conditions
        real(dp),                                intent(in)  :: x1, x2   ! initial and final time
        real(dp),                                intent(in)  :: h        ! step size
        real(dp),   dimension(:),                intent(in)  :: para     ! parameters for derivs

        ! Intent OUT
        real(dp),   dimension(:),   allocatable, intent(out) :: xout
        real(dp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)  :: x        ! time
                real(dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(dp),   dimension(:), intent(in)  :: para     ! parameters for derivs
                real(dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                           :: j            ! counter
        integer(i4)                           :: nstep        ! nuber of steps
        real(dp)                              :: x
        real(dp)                              :: hh, h6, xh
        real(dp),   dimension( size(ystart) ) :: dy, dyt, dym, y, yt

        y(:) = ystart(:)                        ! load starting values
        nstep = nint( (x2-x1)/h, i4 )           ! find number of steps

        if ( allocated(xout) ) deallocate(xout) ! clear out old stored variables if necessary
        if ( allocated(yout) ) deallocate(yout)
        allocate( xout(nstep+1_i4) )            ! allocate storage for saved values
        allocate( yout(nstep+1_i4, size(ystart)) )

        yout(1,:) = y(:)
        xout(1) = x1
        x = x1

        hh = h*0.5_dp
        h6 = h/6.0_dp
        xh = x+hh

        do j=1, nstep                          ! take nstep steps
            call derivs( x, y, para, dy )
            yt = y + hh*dy                     ! first step
            call derivs( xh, yt, para, dyt )         ! second step
            yt = y + hh*dyt
            call derivs( xh, yt, para, dym )         ! third step
            yt = y + h*dym
            dym = dyt + dym
            call derivs( x+h, yt, para, dyt )        ! fourth step
            y = y + h6*( dy+dyt+2.0_dp*dym )   ! accumulate increments with proper weights
            if ( h .lt. tiny(1._dp) )    stop 'RK4_para_dp --> stepsize not significant!'
            x = x+h
            xout(j+1_i4) = x            ! store intermediate steps
            yout(j+1_i4,:) = y(:)
        end do

    end subroutine RK4_para_dp

    ! ------------------------------------------------------------------

    ! SINGLE PRECISION 4th order RUNGE-KUTTA Adaptive Step-size
    subroutine RK4as_para_sp( ystart, x1, x2, h, derivs, para, xout, yout, &   ! obligatory
        hmin, eps )                                                 ! optional

        use mo_kind,   only: i4, sp
        use mo_nrutil, only: reallocate

        implicit none

        ! Intent IN
        real(sp),                                intent(in) :: x1, x2    ! initial and final time
        real(sp),                                intent(in) :: h         ! guessed step size
        real(sp),   dimension(:),                intent(in) :: ystart    ! initial conditions
        real(sp),                   optional,    intent(in) :: hmin
        real(sp),                   optional,    intent(in) :: eps
        real(sp),   dimension(:),                intent(in) :: para      ! parameters for derivs

        ! Intent OUT
        real(sp),   dimension(:),   allocatable, intent(out) :: xout
        real(sp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)  :: x        ! time
                real(sp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(sp),   dimension(:), intent(in)  :: para     ! parameters for derivs
                real(sp),   dimension(:), intent(out) :: dydx     ! derivatives of y
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
        integer(i4),    parameter   :: MAXstp = 1000000_i4  ! max number of steps
        real(sp),       parameter   :: safety = 0.9_sp, pgrow = -0.2_sp, pshrnk = -0.25_sp
        real(sp) :: errcon

        errcon = exp( (1._sp/pgrow)*log(5._sp/safety) )

        if( present(hmin) ) then
            hminIN = hmin
        else
            hminIN = 0.0_sp
        end if

        if( present(eps) ) then
            epsIN = eps
        else
            epsIN = 1e-6_sp
        end if

        x = x1
        hIN = sign( h, x2-x1 )
        nok = 0_i4
        nbad = 0_i4
        kount = 0_i4
        y(:) = ystart(:)
        dxsav = tiny(1._sp)

        xsav = x-2.0_sp*dxsav                                                   ! assures storage of first step
        nullify( xp, yp )
        allocate( xp(256) )
        allocate( yp(size(xp), size(ystart)) )

        call save_a_step                                                        ! save initial step

        do nstep=1, MAXstp                                                      ! take at most MAXstp steps

            call derivs( x, y, para, dydx )
            yscal(:) = abs( y(:) ) + abs( hIN*dydx(:) ) + tiny(1._sp)           ! scaling used to monitor accuracy
                                                                                ! --> CAN BE MODIFIED...

            if ( abs(x-xsav) .gt. abs(dxsav) ) call save_a_step                 ! store intermediate results

            if ( (x+hIN-x2)*(x+hIN-x1) .gt. 0.0_sp ) hIN = x2-x                 ! if stepsize can overshoot, decrease

            do
                call CashKarpRK(y, dydx, x, hIN, ytemp, yerr, derivs, para)     ! take a step
                errmax = maxval( abs(yerr(:)/yscal(:)) )/epsIN                  ! evaluate accuracy
                if ( errmax .lt. 1.0_sp )   exit                                ! step succeeded
                htemp = safety*hIN*(errmax**pshrnk)                             ! truncation error too large, reduce stepsize
                hIN = sign( max( abs(htemp), 0.1_sp*abs(hIN) ), hIN )           ! no more than a factor of 10
                xnew = x+hIN
                if ( abs(xnew-x) .lt. epsilon(1.0_sp) )  stop 'RK4as_para_sp --> hey !!! stepsize underflow!'
            end do

            if ( errmax .gt. errcon )   then                                    ! compute size of next step
                hnext = safety*hIN*(errmax**pgrow)
            else                                                                ! no more than a factor of 5 increase
                hnext = 5.0_sp*hIN
            end if

            hdid = hIN
            x = x+hIN
            y(:) = ytemp(:)

            if ( abs(hdid-hIN) .gt. epsilon(1.0_sp) )    then
                nbad = nbad+1_i4
            else
                nok = nok+1_i4
            end if

            if ( (x-x2)*(x2-x1) .ge. 0.0_sp )   then            ! are we done?!?!
                call save_a_step                                ! save final step
                allocate( xout(kount) )                         ! allocate storage for outputs
                xout(:) = xp(1:kount)
                allocate( yout(kount, size(yp,2)) )             ! allocate storage for outputs
                yout(:,:) = yp(1:kount, :)
                deallocate( xp, yp )                            ! clear out old stored variables
                return                                          ! normal exit
            end if

            if ( abs(hnext-hminIN) .lt. epsilon(1.0_sp) )   stop 'RK4as_para_sp --> WTF! ...stepsize smaller than minimum!!!'
            hIN = hnext

        end do

        stop 'RK4as_para_sp --> dude, too many steps!!!'

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

    end subroutine RK4as_para_sp


    ! DOUBLE PRECISION 4th order RUNGE-KUTTA Adaptive Step-size
    subroutine RK4as_para_dp( ystart, x1, x2, h, derivs, para, xout, yout, &   ! obligatory
        hmin, eps )                                                 ! optional

        use mo_kind,   only: i4, dp
        use mo_nrutil, only: reallocate

        implicit none

        ! Intent IN
        real(dp),                                intent(in) :: x1, x2    ! initial and final time
        real(dp),                                intent(in) :: h         ! guessed step size
        real(dp),   dimension(:),                intent(in) :: ystart    ! initial conditions
        real(dp),                   optional,    intent(in) :: hmin
        real(dp),                   optional,    intent(in) :: eps
        real(dp),   dimension(:),                intent(in) :: para      ! parameters for derivs

        ! Intent OUT
        real(dp),   dimension(:),   allocatable, intent(out) :: xout
        real(dp),   dimension(:,:), allocatable, intent(out) :: yout

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)  :: x        ! time
                real(dp),   dimension(:), intent(in)  :: y        ! unknowns of the equations
                real(dp),   dimension(:), intent(in)  :: para     ! parameters for derivs
                real(dp),   dimension(:), intent(out) :: dydx     ! derivatives of y
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
        real(dp),       parameter   :: safety = 0.9_dp, pgrow = -0.2_dp, pshrnk = -0.25_dp
        real(dp) :: errcon

        errcon = exp( (1._dp/pgrow)*log(5._dp/safety) )

        if( present(hmin) ) then
            hminIN = hmin
        else
            hminIN = 0.0_dp
        end if

        if( present(eps) ) then
            epsIN = eps
        else
            epsIN = 1e-6_dp
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

            call derivs( x, y, para, dydx )
            yscal(:) = abs( y(:) ) + abs( hIN*dydx(:) ) + tiny(1._dp)   ! scaling used to monitor accuracy --> CAN BE MODIFIED...

            if ( abs(x-xsav) .gt. abs(dxsav) ) call save_a_step         ! store intermediate results

            if ( (x+hIN-x2)*(x+hIN-x1) .gt. 0.0_dp ) hIN = x2-x         ! if stepsize can overshoot, decrease

            do
                call CashKarpRK( y, dydx, x, hIN, ytemp, yerr, derivs, para ) ! take a step
                errmax = maxval( abs(yerr(:)/yscal(:)) )/epsIN          ! evaluate accuracy
                if ( errmax .lt. 1.0_dp )   exit                        ! step succeeded
                htemp = safety*hIN*(errmax**pshrnk)                     ! truncation error too large, reduce stepsize
                hIN = sign( max( abs(htemp), 0.1_dp*abs(hIN) ), hIN )   ! no more than a factor of 10
                xnew = x+hIN
!!test
!write(*,*) hIN
                if ( abs(xnew-x) .lt. epsilon(1._dp) )  stop 'RK4as_para_dp --> hey!!! stepsize underflow!'
            end do

            if ( errmax .gt. errcon )   then                            ! compute size of next step
                hnext = safety*hIN*(errmax**pgrow)
            else                                                        ! no more than a factor of 5 increase
                hnext = 5.0_dp*hIN
            end if

            hdid = hIN
            x = x+hIN
            y(:) = ytemp(:)

            if ( abs(hdid-hIN) .gt. epsilon(1._dp) ) then
                nbad = nbad+1_i4
            else
                nok = nok+1_i4
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

            if ( abs(hnext-hminIN) .lt. epsilon(1.0_dp) )   stop 'RK4as_para_dp --> WTF! ...stepsize smaller than minimum!!!'
            hIN = hnext

        end do

        stop 'RK4as_para_dp --> dude, too many steps!!!'

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

    end subroutine RK4as_para_dp


    ! ------------------------------------------------------------------

    ! SINGLE PRECISION Rosenbrock Method for stiff ode's
    subroutine RBstiff_para_sp( ystart, x1, x2, h, derivs, jacobn, para, xout, yout, &      ! obligatory
        hmin, eps )                                                                         ! optional

        use mo_kind,            only : i4, sp
        use mo_utils,           only : ge
        use mo_nrutil,          only : reallocate

        implicit none

        ! Intent IN
        real(sp),                                intent(in)     :: x1, x2    ! initial and final time
        real(sp),                                intent(in)     :: h         ! guessed step size
        real(sp),   dimension(:),                intent(in)     :: ystart    ! initial conditions
        real(sp),                   optional,    intent(in)     :: hmin
        real(sp),                   optional,    intent(in)     :: eps
        real(sp),   dimension(:),                intent(in)     :: para      ! parameters for derivs

        ! Intent OUT
        real(sp),   dimension(:),   allocatable, intent(out)    :: xout
        real(sp),   dimension(:,:), allocatable, intent(out)    :: yout

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)            :: x        ! time
                real(sp),   dimension(:), intent(in)            :: y        ! unknowns of the equations
                real(sp),   dimension(:), intent(in)            :: para     ! parameters for derivs
                real(sp),   dimension(:), intent(out)           :: dydx     ! derivatives of y
            end subroutine derivs

            subroutine jacobn(x, y, para, dfdx, dfdy)
                use mo_kind, only: sp
                implicit none
                real(sp),                   intent(in)          :: x
                real(sp),   dimension(:),   intent(in)          :: y
                real(sp),   dimension(:),   intent(in)          :: para     ! parameters for derivs
                real(sp),   dimension(:),   intent(out)         :: dfdx
                real(sp),   dimension(:,:), intent(out)         :: dfdy
            end subroutine jacobn
        end interface

        ! Internal variables
        integer(i4)                         :: n    
        integer(i4)                         :: nstep                        ! nuber of steps
        integer(i4)                         :: kount                        ! total number of saved steps
        real(sp)                            :: hIN, hdid, hnext, x, hminIN, epsIN, xsav, dxsav
        real(sp),   dimension(size(ystart)) :: dydx, y, yscal, C

        ! Pointers
        real(sp),   dimension(:),   pointer :: xp
        real(sp),   dimension(:,:), pointer :: yp

        ! parameters
        integer(i4),    parameter           :: MAXstp = 1000000_i4          ! max number of steps


        if( present(hmin) ) then
            hminIN = abs(hmin)
        else
            hminIN = 0.0_sp
        end if

        if( present(eps) ) then
            epsIN = abs(eps)
        else
            epsIN = 1e-6_sp
        end if
        
        C       = 1.0_sp                                                !lower boundary for scaling

        hIN     = sign( h, x2-x1 )
        hnext   = 0.0_sp
        hdid    = 0.0_sp

        kount   = 0_i4

        x       = x1
        y(:)    = ystart(:)
        dxsav   = tiny(1._sp)
        xsav    = x-2.0_sp*dxsav                                        ! assures storage of first step

        nullify( xp, yp )
        allocate( xp(256) )
        allocate( yp(size(xp), size(ystart)) )

        call save_a_step                                                ! save initial step

        do nstep=1, MAXstp                                              ! take at most MAXstp steps

            call derivs( x, y, para, dydx )
            forall( n = 1:size(yscal) ) yscal(n) = max (abs(y(n)),C(n)) ! scaling used to monitor accuracy --> CAN BE MODIFIED...

            if ( abs(x-xsav) .gt. abs(dxsav) ) call save_a_step         ! store intermediate results

            if ( (x+hIN-x2)*(x+hIN-x1) .gt. 0.0_sp ) hIN = x2-x         ! if stepsize can overshoot, decrease

            call RBstep( y, dydx, x, hIN, epsIN, yscal, hdid, hnext, derivs, jacobn, para) !do one Rosenbrock-Step

            if ( ge((x-x2)*(x2-x1), 0.0_sp) )   then                    ! are we done?!?!
                call save_a_step                                        ! save final step
                allocate( xout(kount) )                                 ! allocate storage for outputs
                xout(:) = xp(1:kount)
                allocate( yout(kount, size(yp,2)) )                     ! allocate storage for outputs
                yout(:,:) = yp(1:kount, :)
                deallocate( xp, yp )                                    ! clear out old stored variables
                return                                                  ! normal exit
            end if

            if ( abs(hnext-hminIN) .lt. epsilon(1.0_sp) )   stop 'RBstiff_para_sp --> WTF! ...stepsize smaller than minimum!!!'
            
            hIN     = hnext

        end do

        stop 'RBstiff_para_sp --> dude, too many steps!!!'

    contains

        subroutine save_a_step      ! --> like a macro in C
            kount = kount+1_i4
            if ( kount .gt. size(xp) ) then
                xp=>reallocate( xp, 2*size(xp) )
                yp=>reallocate( yp, size(xp), size(yp,2) )
            end if
            xp(kount) = x
            yp(kount, :) = y(:)
            xsav = x
        end subroutine save_a_step

    end subroutine RBstiff_para_sp

    ! DOUBLE PRECISION Rosenbrock Method for stiff ode's
    subroutine RBstiff_para_dp( ystart, x1, x2, h, derivs, jacobn, para, xout, yout, &      ! obligatory
        hmin, eps )                                                                         ! optional

        use mo_kind,            only : i4, dp
        use mo_utils,           only : ge
        use mo_nrutil,          only : reallocate

        implicit none

        ! Intent IN
        real(dp),                                intent(in)     :: x1, x2    ! initial and final time
        real(dp),                                intent(in)     :: h         ! guessed step size
        real(dp),   dimension(:),                intent(in)     :: ystart    ! initial conditions
        real(dp),                   optional,    intent(in)     :: hmin
        real(dp),                   optional,    intent(in)     :: eps
        real(dp),   dimension(:),                intent(in)     :: para      ! parameters for derivs

        ! Intent OUT
        real(dp),   dimension(:),   allocatable, intent(out)    :: xout
        real(dp),   dimension(:,:), allocatable, intent(out)    :: yout

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)            :: x        ! time
                real(dp),   dimension(:), intent(in)            :: y        ! unknowns of the equations
                real(dp),   dimension(:), intent(in)            :: para     ! parameters for derivs
                real(dp),   dimension(:), intent(out)           :: dydx     ! derivatives of y
            end subroutine derivs

            subroutine jacobn(x, y, para, dfdx, dfdy)
                use mo_kind, only: dp
                implicit none
                real(dp),                   intent(in)          :: x
                real(dp),   dimension(:),   intent(in)          :: y
                real(dp),   dimension(:),   intent(in)          :: para     ! parameters for derivs
                real(dp),   dimension(:),   intent(out)         :: dfdx
                real(dp),   dimension(:,:), intent(out)         :: dfdy
            end subroutine jacobn
        end interface

        ! Internal variables
        integer(i4)                         :: n    
        integer(i4)                         :: nstep                        ! nuber of steps
        integer(i4)                         :: kount                        ! total number of saved steps
        real(dp)                            :: hIN, hdid, hnext, x, hminIN, epsIN, xsav, dxsav
        real(dp),   dimension(size(ystart)) :: dydx, y, yscal, C

        ! Pointers
        real(dp),   dimension(:),   pointer :: xp
        real(dp),   dimension(:,:), pointer :: yp

        ! parameters
        integer(i4),    parameter           :: MAXstp = 1000000_i4          ! max number of steps


        if( present(hmin) ) then
            hminIN = abs(hmin)
        else
            hminIN = 0.0_dp
        end if

        if( present(eps) ) then
            epsIN = abs(eps)
        else
            epsIN = 1e-6_dp
        end if
        
        C       = 1.0_dp                                                !lower boundary for scaling

        hIN     = sign( h, x2-x1 )
        hnext   = 0.0_dp
        hdid    = 0.0_dp

        kount   = 0_i4

        x       = x1
        y(:)    = ystart(:)
        dxsav   = tiny(1._dp)
        xsav    = x-2.0_dp*dxsav                                        ! assures storage of first step

        nullify( xp, yp )
        allocate( xp(256) )
        allocate( yp(size(xp), size(ystart)) )

        call save_a_step                                                ! save initial step

        do nstep=1, MAXstp                                              ! take at most MAXstp steps

            call derivs( x, y, para, dydx )
            forall( n = 1:size(yscal) ) yscal(n) = max (abs(y(n)),C(n)) ! scaling used to monitor accuracy --> CAN BE MODIFIED...

            if ( abs(x-xsav) .gt. abs(dxsav) ) call save_a_step         ! store intermediate results

            if ( (x+hIN-x2)*(x+hIN-x1) .gt. 0.0_dp ) hIN = x2-x         ! if stepsize can overshoot, decrease

            call RBstep( y, dydx, x, hIN, epsIN, yscal, hdid, hnext, derivs, jacobn, para) !do one Rosenbrock-Step

            if ( ge((x-x2)*(x2-x1), 0.0_dp) )   then                    ! are we done?!?!
                call save_a_step                                        ! save final step
                allocate( xout(kount) )                                 ! allocate storage for outputs
                xout(:) = xp(1:kount)
                allocate( yout(kount, size(yp,2)) )                     ! allocate storage for outputs
                yout(:,:) = yp(1:kount, :)
                deallocate( xp, yp )                                    ! clear out old stored variables
                return                                                  ! normal exit
            end if

            if ( abs(hnext-hminIN) .lt. epsilon(1.0_dp) )   stop 'RBstiff_para_dp --> WTF! ...stepsize smaller than minimum!!!'
            
            hIN     = hnext

        end do

        stop 'RBstiff_para_dp --> dude, too many steps!!!'

    contains

        subroutine save_a_step      ! --> like a macro in C
            kount = kount+1_i4
            if ( kount .gt. size(xp) ) then
                xp=>reallocate( xp, 2*size(xp) )
                yp=>reallocate( yp, size(xp), size(yp,2) )
            end if
            xp(kount) = x
            yp(kount, :) = y(:)
            xsav = x
        end subroutine save_a_step

    end subroutine RBstiff_para_dp

    ! ------------------------------------------------------------------
    !   PRIVATE METHODS
    ! ------------------------------------------------------------------

    ! SINGLE PRECISION CASH-KARP RUNGE-KUTTA step
    subroutine CashKarpRK_para_sp( y, dydx, x, h, yout, yerr, derivs, para )
        ! Given values for N variables y and their derivatives dydx known at x, use the fifth-order
        ! Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
        ! the incremented variables as yout. Also return an estimate of the local truncation error
        ! in yout using the embedded fourth order method.

        use mo_kind,   only : i4, sp
        use mo_nrutil, only : assert_eq

        implicit none

        ! Intent IN
        real(sp),   dimension(:), intent(in)  :: y, dydx
        real(sp),                 intent(in)  :: x, h
        real(sp),   dimension(:), intent(in)  :: para     ! parameters for derivs

        ! Intent OUT
        real(sp),   dimension(:), intent(out) :: yout, yerr

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: sp
                implicit none
                real(sp),                 intent(in)  :: x
                real(sp),   dimension(:), intent(in)  :: y
                real(sp),   dimension(:), intent(in)  :: para     ! parameters for derivs
                real(sp),   dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: ndum
        real(sp),   dimension(:),   allocatable :: ak2, ak3, ak4, ak5, ak6, ytemp

        ! parameters
        real(sp),   parameter   :: A2=0.2_sp, A3=0.3_sp, A4=0.6_sp, A5=1._sp, A6=0.875_sp, B21=0.2_sp, B31=3._sp/40._sp, &
            B32=9._sp/40._sp, B41=0.3_sp, B42=-0.9_sp, B43=1.2_sp, B51=-11._sp/54._sp, B52=2.5_sp, &
            B53=-70._sp/27._sp, B54=35._sp/27._sp, B61=1631._sp/55296._sp, B62=175._sp/512._sp, &
            B63=575._sp/13824._sp, B64=44275._sp/110592._sp, B65=253._sp/4096._sp, &
            C1=37._sp/378._sp, C3=250._sp/621._sp, C4=125._sp/594._sp, C6=512._sp/1771._sp, &
            DC1=C1-2825._sp/27648._sp, DC3=C3-18575._sp/48384._sp, &
            DC4=C4-13525._sp/55296._sp, DC5=-277._sp/14336._sp, DC6=C6-.25_sp

        ndum = assert_eq( size(y), size(dydx) ,size(yout) ,size(yerr) , 'CashKarpRK_para_sp' )

        allocate( ak2(ndum), ak3(ndum), ak4(ndum), ak5(ndum), ak6(ndum), ytemp(ndum) )

        ytemp = y + B21*h*dydx                                  ! first step
        call derivs( x+A2*h, ytemp, para, ak2 )                       ! second step
        ytemp = y + h*(B31*dydx+B32*ak2)
        call derivs( x+A3*h, ytemp, para, ak3 )                       ! third step
        ytemp = y + h*(B41*dydx+B42*ak2+B43*ak3)
        call derivs( x+A4*h, ytemp, para, ak4 )                       ! fourth step
        ytemp = y + h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
        call derivs( x+A5*h, ytemp, para, ak5 )                       ! fifth step
        ytemp = y + h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
        call derivs( x+A6*h, ytemp, para, ak6 )                       ! sixth step
        yout = y + h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)             ! accumulate increments with proper weights
        yerr = h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)

    end subroutine CashKarpRK_para_sp


    ! DOUBLE PRECISION CASH-KARP RUNGE-KUTTA step
    subroutine CashKarpRK_para_dp( y, dydx, x, h, yout, yerr, derivs, para )
        ! Given values for N variables y and their derivatives dydx known at x, use the fifth-order
        ! Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
        ! the incremented variables as yout. Also return an estimate of the local truncation error
        ! in yout using the embedded fourth order method.

        use mo_kind,   only: i4, dp
        use mo_nrutil, only: assert_eq

        implicit none

        ! Intent IN
        real(dp),   dimension(:), intent(in)  :: y, dydx
        real(dp),                 intent(in)  :: x, h
        real(dp),   dimension(:), intent(in)  :: para     ! parameters for derivs

        ! Intent OUT
        real(dp),   dimension(:), intent(out) :: yout, yerr

        interface
            subroutine derivs( x, y, para, dydx )
                use mo_kind, only: dp
                implicit none
                real(dp),                 intent(in)  :: x
                real(dp),   dimension(:), intent(in)  :: y
                real(dp),   dimension(:), intent(in)  :: para     ! parameters for derivs
                real(dp),   dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        ! Internal variables
        integer(i4)                             :: ndum
        real(dp),   dimension(:),   allocatable :: ak2, ak3, ak4, ak5, ak6, ytemp

        ! parameters
        real(dp),   parameter   :: A2=0.2_dp, A3=0.3_dp, A4=0.6_dp, A5=1._dp, A6=0.875_dp, B21=0.2_dp, B31=3._dp/40._dp, &
            B32=9._dp/40._dp, B41=0.3_dp, B42=-0.9_dp, B43=1.2_dp, B51=-11._dp/54._dp, B52=2.5_dp, &
            B53=-70._dp/27._dp, B54=35._dp/27._dp, B61=1631._dp/55296._dp, B62=175._dp/512._dp, &
            B63=575._dp/13824._dp, B64=44275._dp/110592._dp, B65=253._dp/4096._dp, &
            C1=37._dp/378._dp, C3=250._dp/621._dp, C4=125._dp/594._dp, C6=512._dp/1771._dp, &
            DC1=C1-2825._dp/27648._dp, DC3=C3-18575._dp/48384._dp, &
            DC4=C4-13525._dp/55296._dp, DC5=-277._dp/14336._dp, DC6=C6-.25_dp

        ndum = assert_eq( size(y), size(dydx) ,size(yout) ,size(yerr) , 'CashKarpRK_para_dp' )

        allocate( ak2(ndum), ak3(ndum), ak4(ndum), ak5(ndum), ak6(ndum), ytemp(ndum) )

        ytemp = y + B21*h*dydx                                  ! first step
        call derivs( x+A2*h, ytemp, para, ak2 )                       ! second step
        ytemp = y + h*(B31*dydx+B32*ak2)
        call derivs( x+A3*h, ytemp, para, ak3 )                       ! third step
        ytemp = y + h*(B41*dydx+B42*ak2+B43*ak3)
        call derivs( x+A4*h, ytemp, para, ak4 )                       ! fourth step
        ytemp = y + h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
        call derivs( x+A5*h, ytemp, para, ak5 )                       ! fifth step
        ytemp = y + h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
        call derivs( x+A6*h, ytemp, para, ak6 )                       ! sixth step
        yout = y + h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)             ! accumulate increments with proper weights
        yerr = h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)

    end subroutine CashKarpRK_para_dp


    ! SINGLE PRECISION ROSENBROCK step
    subroutine RBstep_para_sp(y, dydx, x, htry, eps, yscal, hdid, hnext, derivs, jacobn, para)

        use mo_kind,            only : i4, sp
        use mo_utils,           only : eq, le
        use mo_nrutil,          only : assert_eq, diagadd
        use mo_linear_algebra,  only : solve_linear_equations

        implicit none

        real(sp),       dimension(:),   intent(inout)   :: y
        real(sp),       dimension(:),   intent(in)      :: dydx, yscal
        real(sp),                       intent(inout)   :: x
        real(sp),       dimension(:),   intent(in)      :: para     ! parameters for derivs

        real(sp),                       intent(in)      :: htry, eps
        real(sp),                       intent(out)     :: hdid, hnext

        interface
            subroutine derivs(x, y, para, dydx)
                use mo_kind, only: sp
                implicit none
                real(sp),                   intent(in)  :: x
                real(sp),   dimension(:),   intent(in)  :: y
                real(sp),   dimension(:),   intent(in)  :: para     ! parameters for derivs
                real(sp),   dimension(:),   intent(out) :: dydx
            end subroutine derivs

            subroutine jacobn(x, y, para, dfdx, dfdy)
                use mo_kind, only: sp
                implicit none
                real(sp),                   intent(in)  :: x
                real(sp),   dimension(:),   intent(in)  :: y
                real(sp),   dimension(:),   intent(in)  :: para     ! parameters for derivs
                real(sp),   dimension(:),   intent(out) :: dfdx
                real(sp),   dimension(:,:), intent(out) :: dfdy
            end subroutine jacobn
        end interface

        integer(i4)                                     :: jtry,ndum
!        integer(i4),    dimension(size(y))              :: indx
        real(sp),       dimension(size(y))              :: dfdx,dytmp,err,g1,g2,g3,g4,ysav
        real(sp),       dimension(size(y),size(y))      :: a,dfdy
        real(sp)                                        :: errmax,h,xsav !,d

        integer(i4), parameter  ::  maxtry  = 40_i4
        real(sp), parameter     ::  safety  = 0.9_sp,               grow    = 1.5_sp,           pgrow   = -0.25_sp,&
                                    shrnk   = 0.5_sp,               pshrnk  = -1.0_sp/3.0_sp,   errcon  = 0.1296_sp,&
                                    gam     = 1.0_sp/2.0_sp,&
                                    a21     = 2.0_sp,               a31     = 48.0_sp/25.0_sp,  a32     = 6.0_sp/25.0_sp,&
                                    c21     = -8.0_sp,              c31     = 372.0_sp/25.0_sp, c32     = 12.0_sp/5.0_sp,&
                                    c41     = -112.0_sp/125.0_sp,   c42     = -54._sp/125.0_sp, c43     = -2.0_sp/5.0_sp,&
                                    b1      = 19.0_sp/9.0_sp,       b2      = 1.0_sp/2.0_sp,    b3      = 25.0_sp/108.0_sp,&
                                    b4      = 125.0_sp/108.0_sp,&
                                    e1      = 17.0_sp/54.0_sp,      e2      = 7.0_sp/36.0_sp,   e3      = 0.0_sp,&
                                    e4      = 125.0_sp/108.0_sp,&
                                    c1x     = 1.0_sp/2.0_sp,        c2x     = -3.0_sp/2.0_sp,   c3x     = 121.0_sp/50.0_sp,&
                                    c4x     = 29.0_sp/250.0_sp,&
                                    a2x     = 1.0_sp,               a3x     = 3.0_sp/5.0_sp


        ndum        =   assert_eq(size(y),size(dydx),size(yscal),'stiff')
        xsav        =   x
        ysav(:)     =   y(:)

        call jacobn(xsav,ysav,para,dfdx,dfdy)

        h           =   htry

        do jtry=1,maxtry

            a(:,:)  =   -dfdy(:,:)

            call diagadd(a,1.0_sp/(gam*h))
!            call ludcmp(a,indx,d)

            g1      =   dydx+h*c1x*dfdx
            g1      =   solve_linear_equations(a,g1)

 !           call lubksb(a,indx,g1)

            y       =   ysav+a21*g1
            x       =   xsav+a2x*h

            call derivs(x,y,para,dytmp)

            g2      =   dytmp+h*c2x*dfdx+c21*g1/h
            g2      =   solve_linear_equations(a,g2)

  !          call lubksb(a,indx,g2)

            y       =   ysav+a31*g1+a32*g2
            x       =   xsav+a3x*h

            call derivs(x,y,para,dytmp)

            g3      =   dytmp+h*c3x*dfdx+(c31*g1+c32*g2)/h
            g3      =   solve_linear_equations(a,g3)

   !         call lubksb(a,indx,g3)

            g4      =   dytmp+h*c4x*dfdx+(c41*g1+c42*g2+c43*g3)/h
            g4      =   solve_linear_equations(a,g4)

    !        call lubksb(a,indx,g4)

            y       =   ysav+b1*g1+b2*g2+b3*g3+b4*g4
            err     =   e1*g1+e2*g2+e3*g3+e4*g4
            x       =   xsav+h

            if (eq(x,xsav)) stop 'stepsize not significant in RBstep_para_sp'

            errmax  =   maxval(abs(err/yscal))/eps

            if (le(errmax,1.0_sp)) then
                hdid    =   h
                hnext   =  merge(safety*h*errmax**pgrow, grow*h, errmax > errcon)
                return
            else
                hnext   =   safety*h*errmax**pshrnk
                h       =   sign(max(abs(hnext),shrnk*abs(h)),h)
            end if

        end do

        stop 'exceeded maxtry in RBstep_para_sp'

    end subroutine RBstep_para_sp


    ! DOUBLE PRECISION ROSENBROCK step
    subroutine RBstep_para_dp(y, dydx, x, htry, eps, yscal, hdid, hnext, derivs, jacobn, para)

        use mo_kind,            only : i4, dp
        use mo_utils,           only : eq, le
        use mo_nrutil,          only : assert_eq, diagadd
        use mo_linear_algebra,  only : solve_linear_equations

        implicit none

        real(dp),       dimension(:),   intent(inout)   :: y
        real(dp),       dimension(:),   intent(in)      :: dydx, yscal
        real(dp),                       intent(inout)   :: x
        real(dp),       dimension(:),   intent(in)      :: para     ! parameters for derivs

        real(dp),                       intent(in)      :: htry, eps
        real(dp),                       intent(out)     :: hdid, hnext

        interface
            subroutine derivs(x, y, para, dydx)
                use mo_kind, only: dp
                implicit none
                real(dp),                   intent(in)  :: x
                real(dp),   dimension(:),   intent(in)  :: y
                real(dp),   dimension(:),   intent(in)  :: para     ! parameters for derivs
                real(dp),   dimension(:),   intent(out) :: dydx
            end subroutine derivs

            subroutine jacobn(x, y, para, dfdx, dfdy)
                use mo_kind, only: dp
                implicit none
                real(dp),                   intent(in)  :: x
                real(dp),   dimension(:),   intent(in)  :: y
                real(dp),   dimension(:),   intent(in)  :: para     ! parameters for derivs
                real(dp),   dimension(:),   intent(out) :: dfdx
                real(dp),   dimension(:,:), intent(out) :: dfdy
            end subroutine jacobn
        end interface

        integer(i4)                                     :: jtry,ndum
!        integer(i4),    dimension(size(y))              :: indx
        real(dp),       dimension(size(y))              :: dfdx,dytmp,err,g1,g2,g3,g4,ysav
        real(dp),       dimension(size(y),size(y))      :: a,dfdy
        real(dp)                                        :: errmax,h,xsav !,d

        integer(i4), parameter  ::  maxtry  = 40_i4
        real(dp), parameter     ::  safety  = 0.9_dp,               grow    = 1.5_dp,           pgrow   = -0.25_dp,&
                                    shrnk   = 0.5_dp,               pshrnk  = -1.0_dp/3.0_dp,   errcon  = 0.1296_dp,&
                                    gam     = 1.0_dp/2.0_dp,&
                                    a21     = 2.0_dp,               a31     = 48.0_dp/25.0_dp,  a32     = 6.0_dp/25.0_dp,&
                                    c21     = -8.0_dp,              c31     = 372.0_dp/25.0_dp, c32     = 12.0_dp/5.0_dp,&
                                    c41     = -112.0_dp/125.0_dp,   c42     = -54._dp/125.0_dp, c43     = -2.0_dp/5.0_dp,&
                                    b1      = 19.0_dp/9.0_dp,       b2      = 1.0_dp/2.0_dp,    b3      = 25.0_dp/108.0_dp,&
                                    b4      = 125.0_dp/108.0_dp,&
                                    e1      = 17.0_dp/54.0_dp,      e2      = 7.0_dp/36.0_dp,   e3      = 0.0_dp,&
                                    e4      = 125.0_dp/108.0_dp,&
                                    c1x     = 1.0_dp/2.0_dp,        c2x     = -3.0_dp/2.0_dp,   c3x     = 121.0_dp/50.0_dp,&
                                    c4x     = 29.0_dp/250.0_dp,&
                                    a2x     = 1.0_dp,               a3x     = 3.0_dp/5.0_dp


        ndum        =   assert_eq(size(y),size(dydx),size(yscal),'stiff')
        xsav        =   x
        ysav(:)     =   y(:)

        call jacobn(xsav,ysav,para,dfdx,dfdy)

        h           =   htry

        do jtry=1,maxtry

            a(:,:)  =   -dfdy(:,:)

            call diagadd(a,1.0_dp/(gam*h))
!            call ludcmp(a,indx,d)

            g1      =   dydx+h*c1x*dfdx
            g1      =   solve_linear_equations(a,g1)

 !           call lubksb(a,indx,g1)

            y       =   ysav+a21*g1
            x       =   xsav+a2x*h

            call derivs(x,y,para,dytmp)

            g2      =   dytmp+h*c2x*dfdx+c21*g1/h
            g2      =   solve_linear_equations(a,g2)

  !          call lubksb(a,indx,g2)

            y       =   ysav+a31*g1+a32*g2
            x       =   xsav+a3x*h

            call derivs(x,y,para,dytmp)

            g3      =   dytmp+h*c3x*dfdx+(c31*g1+c32*g2)/h
            g3      =   solve_linear_equations(a,g3)

   !         call lubksb(a,indx,g3)

            g4      =   dytmp+h*c4x*dfdx+(c41*g1+c42*g2+c43*g3)/h
            g4      =   solve_linear_equations(a,g4)

    !        call lubksb(a,indx,g4)

            y       =   ysav+b1*g1+b2*g2+b3*g3+b4*g4
            err     =   e1*g1+e2*g2+e3*g3+e4*g4
            x       =   xsav+h

            if (eq(x,xsav)) stop 'stepsize not significant in RBstep_para_dp'

            errmax  =   maxval(abs(err/yscal))/eps

            if (le(errmax,1.0_dp)) then
                hdid    =   h
                hnext   =  merge(safety*h*errmax**pgrow, grow*h, errmax > errcon)
                return
            else
                hnext   =   safety*h*errmax**pshrnk
                h       =   sign(max(abs(hnext),shrnk*abs(h)),h)
            end if

        end do

        stop 'exceeded maxtry in RBstep_para_dp'

    end subroutine RBstep_para_dp

    ! ------------------------------------------------------------------


end module mo_ode_solver
