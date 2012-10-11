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

    use mo_kind, only : i4, sp, dp, lgt

    implicit none

    private

    public :: EULERdumb     ! Euler method with equal time-steps increments
    public :: RKdumb        ! Fourth-order Runge-Kutta method with equal time-steps increments
    public :: RKas          ! Runge-Kutta method with adaptive stepsize control

    ! Interfaces for single and double precision routines; sort alphabetically
    interface EULER
        module procedure EULER_sp, EULER_dp
    end interface

    interface EULERdumb
        module procedure EULERdumb_sp, EULERdumb_dp
    end interface

    interface RK4
        module procedure RK4_sp, RK4_dp
    end interface

    interface RKas
        module procedure RKas_sp, RKas_dp
    end interface

    interface RKck
        module procedure RKck_sp, RKck_dp
    end interface

    interface RKdumb
        module procedure RKdumb_sp, RKdumb_dp
    end interface

    interface RKqs
        module procedure RKqs_sp, RKqs_dp
    end interface

    interface header
        module procedure header_sp, header_dp
    end interface

    ! storage variables (MODULE VARIABLES)
    integer(i4) :: nok, nbad                            ! number of good and bad (but retried and fixed) steps taken
    integer(i4) :: kount                                ! total number of saved steps
    logical(lgt), save :: dump = .false.                ! if set to true, then results are stored in a file
    logical(lgt), save :: save_steps = .false.          ! if set to true, then intermediate values are stored
                                                        ! in xp_sp/dp and yp_sp/dp at intervals greater than dxsav
    ! SINGLE PRECISION
    real(sp), dimension(:), allocatable :: xx_sp        ! varible for storing intermediate steps
    real(sp), dimension(:,:), allocatable :: yy_sp      ! varible for storing intermediate steps
    real(sp), dimension(:), pointer :: xp_sp
    real(sp), dimension(:,:), pointer :: yp_sp
    ! DOUBLE PRECISION
    real(dp), dimension(:), allocatable :: xx_dp
    real(dp), dimension(:,:), allocatable :: yy_dp
    real(dp), dimension(:), pointer :: xp_dp
    real(dp), dimension(:,:), pointer :: yp_dp

    ! ------------------------------------------------------------------

contains

    ! ------------------------------------------------------------------

    !     NAME
    !         EULERdumb, RKdumb

    !     PURPOSE
    !         Integration of Ordinary Differential Equations using Euler or fourth-order Runge-Kutta methods
    !         with fixed time-steps increments.
    !         Starting from N initial values vstart known at x1, use Euler or Runge-Kutta to advance nstep equal increments to x2.
    !         Results are stored in the variables xout and yout.

    !     CALLING SEQUENCE
    !         call EULERdumb( vstart, x1, x2, h, derivs, Tout, xout, yout, fileName )
    !         call RKdumb( vstart, x1, x2, h, derivs, Tout, xout, yout, fileName )
    !
    !         --> see also example in test directory --> test/test_mo_ode_solver

    !     INDENT(IN)
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

    !     INDENT(INOUT)
    !         none

    !     INDENT(OUT)
    !         none

    !     INDENT(IN), OPTIONAL
    !         character(*)  ::  fileName    ... name of the file in which all the results will be stored
    !                                           DEFAULT: no output files will be created
    !
    !     INDENT(INOUT), OPTIONAL
    !         none

    !     INDENT(OUT), OPTIONAL
    !         real(sp/dp)                                  ::  Tout     ... print out the elapsed CPU time
    !                                                                       DEFAULT: no outputs will be created
    !
    !         real(sp/dp),  dimension(:),   allocatable    ::  xout     ... storage for outputs (time)
    !                                                                       DEFAULT: no outputs will be created
    !
    !         real(sp/dp),  dimension(:,:), allocatable    ::  yout     ... storage for outputs (incremented variables)
    !                                                                       DEFAULT: no outputs will be created

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
    !         Modified, ....

    ! SINGLE PRECISION
    subroutine EULERdumb_sp( vstart, x1, x2, h, derivs, &       ! obligatory
        Tout, xout, yout,          &       ! optional out
        fileName                   &       ! optional in
        )

        use mo_kind, only : i4, sp
        use mo_nrutil, only : nrerror

        implicit none

        real(sp), dimension(:), intent(in) :: vstart
        real(sp), intent(in) :: x1, x2, h
        character(*), intent(in), optional :: fileName
        real(sp), intent(out), optional :: Tout
        real(sp), dimension(:), allocatable, intent(out), optional :: xout
        real(sp), dimension(:,:), allocatable, intent(out), optional :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : sp
                implicit none
                real(sp), intent(in) :: x                       ! time
                real(sp), dimension(:), intent(in) :: y         ! unknowns of the equations
                real(sp), dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        integer(i4) :: k, nstep
        real(sp) :: x, T1, T2
        real(sp), dimension( size(vstart) ) :: dv, v, vout

        call header( x1, x2, h, vstart, 'EULER method' )    ! header with settings informations

        if( present(fileName) ) dump = .true.

        v(:) = vstart(:)                        ! load starting values
        nstep = int( (x2-x1)/h, i4 )            ! find number of steps

        if ( allocated(xx_sp) ) deallocate(xx_sp)   ! clear out old stored variables if necessary
        if ( allocated(yy_sp) ) deallocate(yy_sp)
        allocate( xx_sp(nstep+1) )                  ! allocate storage for saved values
        allocate( yy_sp(size(vstart), nstep+1) )

        yy_sp(:,1) = v(:)
        xx_sp(1) = x1
        x = x1
        T1 = 0._sp
        T2 = T1

        if( present(Tout) ) call cpu_time(T1)
        do k=1,nstep                            ! take nstep steps
            call derivs( x, v, dv )
            call EULER( v, dv, x, h, vout, derivs )
            v(:) = vout(:)
            if (x+h == x) call nrerror( 'stepsize not significant in EULERdumb' )
            x = x+h
            xx_sp(k+1) = x                      ! store intermediate steps
            yy_sp(:,k+1) = v(:)
        end do
        if( present(Tout) ) then
            call cpu_time(T2)
            Tout = T2 - T1
            write(*,*) 'Elapsed CPU time EULER = ', Tout
        end if

        if(dump) then
            open( unit = 10, file = fileName, action = "write", status = "replace" )      ! dumping results
            do k=1, size(xx_sp)
                write(10,*) xx_sp(k), yy_sp(:,k)
            end do
            close( unit = 10 )
            write(*,*) ''
            write(*,*) fileName, ' was successfully saved!'
        end if

        if( present(xout) ) then
            allocate( xout(nstep+1) )                   ! allocate storage for outputs
            xout(:) = xx_sp(:)
        end if
        if( present(yout) ) then
            allocate( yout(size(vstart), nstep+1) )     ! allocate storage for outputs
            yout(:,:) = yy_sp(:,:)
        end if

    end subroutine EULERdumb_sp


    ! DOUBLE PRECISION
    subroutine EULERdumb_dp( vstart, x1, x2, h, derivs, &       ! obligatory
        Tout, xout, yout,          &       ! optional out
        fileName                   &       ! optional in
        )

        use mo_kind, only : i4, dp
        use mo_nrutil, only : nrerror

        implicit none

        real(dp), dimension(:), intent(in) :: vstart
        real(dp), intent(in) :: x1, x2, h
        character(*), intent(in), optional :: fileName
        real(dp), intent(out), optional :: Tout
        real(dp), dimension(:), allocatable, intent(out), optional :: xout
        real(dp), dimension(:,:), allocatable, intent(out), optional :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : dp
                implicit none
                real(dp), intent(in) :: x                       ! time
                real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
                real(dp), dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        integer(i4) :: k, nstep
        real(dp) :: x, T1, T2
        real(dp), dimension( size(vstart) ) :: dv, v, vout

        call header( x1, x2, h, vstart, 'EULER method' )    ! header with settings informations

        if( present(fileName) ) dump = .true.

        v(:) = vstart(:)                        ! load starting values
        nstep = int( (x2-x1)/h, i4 )            ! find number of steps

        if ( allocated(xx_dp) ) deallocate(xx_dp)   ! clear out old stored variables if necessary
        if ( allocated(yy_dp) ) deallocate(yy_dp)
        allocate( xx_dp(nstep+1) )                  ! allocate storage for saved values
        allocate( yy_dp(size(vstart), nstep+1) )

        yy_dp(:,1) = v(:)
        xx_dp(1) = x1
        x = x1
        T1 = 0._dp
        T2 = T1

        if( present(Tout) ) call cpu_time(T1)
        do k=1,nstep                            ! take nstep steps
            call derivs( x, v, dv )
            call EULER( v, dv, x, h, vout, derivs )
            v(:) = vout(:)
            if (x+h == x) call nrerror( 'stepsize not significant in EULERdumb' )
            x = x+h
            xx_dp(k+1) = x                      ! store intermediate steps
            yy_dp(:,k+1) = v(:)
        end do
        if( present(Tout) ) then
            call cpu_time(T2)
            Tout = T2 - T1
            write(*,*) 'Elapsed CPU time EULER = ', Tout
        end if

        if(dump) then
            open( unit = 10, file = fileName, action = "write", status = "replace" )      ! dumping results
            do k=1, size(xx_dp)
                write(10,*) xx_dp(k), yy_dp(:,k)
            end do
            close( unit = 10 )
            write(*,*) ''
            write(*,*) fileName, ' was successfully saved!'
        end if

        if( present(xout) ) then
            allocate( xout(nstep+1) )                   ! allocate storage for outputs
            xout(:) = xx_dp(:)
        end if
        if( present(yout) ) then
            allocate( yout(size(vstart), nstep+1) )     ! allocate storage for outputs
            yout(:,:) = yy_dp(:,:)
        end if

    end subroutine EULERdumb_dp


    ! SINGLE PRECISION
    subroutine RKdumb_sp( vstart, x1, x2, h, derivs, &       ! obligatory
        Tout, xout, yout,          &       ! optional out
        fileName                   &       ! optional in
        )

        use mo_kind, only : i4, sp
        use mo_nrutil, only : nrerror

        implicit none

        real(sp), dimension(:), intent(in) :: vstart
        real(sp), intent(in) :: x1, x2, h
        character(*), intent(in), optional :: fileName
        real(sp), intent(out), optional :: Tout
        real(sp), dimension(:), allocatable, intent(out), optional :: xout
        real(sp), dimension(:,:), allocatable, intent(out), optional :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : sp
                implicit none
                real(sp), intent(in) :: x                       ! time
                real(sp), dimension(:), intent(in) :: y         ! unknowns of the equations
                real(sp), dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        integer(i4) :: k, nstep
        real(sp) :: x, T1, T2
        real(sp), dimension( size(vstart) ) :: dv, v, vout

        call header( x1, x2, h, vstart, 'RUNGE-KUTTA method' )    ! header with settings informations

        if( present(fileName) ) dump = .true.

        v(:) = vstart(:)                        ! load starting values
        nstep = int( (x2-x1)/h, i4 )            ! find number of steps

        if ( allocated(xx_sp) ) deallocate(xx_sp)   ! clear out old stored variables if necessary
        if ( allocated(yy_sp) ) deallocate(yy_sp)
        allocate( xx_sp(nstep+1) )                  ! allocate storage for saved values
        allocate( yy_sp(size(vstart), nstep+1) )

        yy_sp(:,1) = v(:)
        xx_sp(1) = x1
        x = x1
        T1 = 0._sp
        T2 = T1

        if( present(Tout) ) call cpu_time(T1)
        do k=1,nstep                            ! take nstep steps
            call derivs( x, v, dv )
            call RK4( v, dv, x, h, vout, derivs )
            v(:) = vout(:)
            if (x+h == x) call nrerror( 'stepsize not significant in RKdumb' )
            x = x+h
            xx_sp(k+1) = x                      ! store intermediate steps
            yy_sp(:,k+1) = v(:)
        end do
        if( present(Tout) ) then
            call cpu_time(T2)
            Tout = T2 - T1
            write(*,*) 'Elapsed CPU time RUNGE KUTTA = ', Tout
        end if

        if(dump) then
            open( unit = 10, file = fileName, action = "write", status = "replace" )      ! dumping results
            do k=1, size(xx_sp)
                write(10,*) xx_sp(k), yy_sp(:,k)
            end do
            close( unit = 10 )
            write(*,*) ''
            write(*,*) fileName, ' was successfully saved!'
        end if

        if( present(xout) ) then
            allocate( xout(nstep+1) )                   ! allocate storage for outputs
            xout(:) = xx_sp(:)
        end if
        if( present(yout) ) then
            allocate( yout(size(vstart), nstep+1) )     ! allocate storage for outputs
            yout(:,:) = yy_sp(:,:)
        end if

    end subroutine RKdumb_sp


    ! DOUBLE PRECISION
    subroutine RKdumb_dp( vstart, x1, x2, h, derivs, &       ! obligatory
        Tout, xout, yout,          &       ! optional out
        fileName                   &       ! optional in
        )

        use mo_kind, only : i4, dp
        use mo_nrutil, only : nrerror

        implicit none

        real(dp), dimension(:), intent(in) :: vstart
        real(dp), intent(in) :: x1, x2, h
        character(*), intent(in), optional :: fileName
        real(dp), intent(out), optional :: Tout
        real(dp), dimension(:), allocatable, intent(out), optional :: xout
        real(dp), dimension(:,:), allocatable, intent(out), optional :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : dp
                implicit none
                real(dp), intent(in) :: x                       ! time
                real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
                real(dp), dimension(:), intent(out) :: dydx     ! derivatives of y
            end subroutine derivs
        end interface

        integer(i4) :: k, nstep
        real(dp) :: x, T1, T2
        real(dp), dimension( size(vstart) ) :: dv, v, vout

        call header( x1, x2, h, vstart, 'RUNGE-KUTTA method' )    ! header with settings informations

        if( present(fileName) ) dump = .true.

        v(:) = vstart(:)                        ! load starting values
        nstep = int( (x2-x1)/h, i4 )            ! find number of steps

        if ( allocated(xx_dp) ) deallocate(xx_dp)   ! clear out old stored variables if necessary
        if ( allocated(yy_dp) ) deallocate(yy_dp)
        allocate( xx_dp(nstep+1) )                  ! allocate storage for saved values
        allocate( yy_dp(size(vstart), nstep+1) )

        yy_dp(:,1) = v(:)
        xx_dp(1) = x1
        x = x1
        T1 = 0._dp
        T2 = T1

        if( present(Tout) ) call cpu_time(T1)
        do k=1,nstep                            ! take nstep steps
            call derivs( x, v, dv )
            call RK4( v, dv, x, h, vout, derivs )
            v(:) = vout(:)
            if (x+h == x) call nrerror( 'stepsize not significant in RKdumb' )
            x = x+h
            xx_dp(k+1) = x                      ! store intermediate steps
            yy_dp(:,k+1) = v(:)
        end do
        if( present(Tout) ) then
            call cpu_time(T2)
            Tout = T2 - T1
            write(*,*) 'Elapsed CPU time RUNGE KUTTA = ', Tout
        end if

        if(dump) then
            open( unit = 10, file = fileName, action = "write", status = "replace" )      ! dumping results
            do k=1, size(xx_dp)
                write(10,*) xx_dp(k), yy_dp(:,k)
            end do
            close( unit = 10 )
            write(*,*) ''
            write(*,*) fileName, ' was successfully saved!'
        end if

        if( present(xout) ) then
            allocate( xout(nstep+1) )                   ! allocate storage for outputs
            xout(:) = xx_dp(:)
        end if
        if( present(yout) ) then
            allocate( yout(size(vstart), nstep+1) )     ! allocate storage for outputs
            yout(:,:) = yy_dp(:,:)
        end if

    end subroutine RKdumb_dp


    ! ------------------------------------------------------------------

    !     NAME
    !         RKas

    !     PURPOSE
    !         Integration of Ordinary Differential Equations using Runge-Kutta with adaptive stepsize control.
    !         Integrate the array of starting values ystart from x1 to x2 with accuracy eps, storing intermediate results
    !         in the module variables. h1 should be set as a guessed first stepsize,
    !         hmin as the minimum allowed stepsize (can be zero).
    !         On output ystart is replaced by values at the end of the integration interval.

    !     CALLING SEQUENCE
    !         call RKas_sp( ystart, x1, x2, h1, derivs, Tout, xout, yout, fileName, hminIN, epsIN )
    !
    !         --> see also example in test directory --> test/test_mo_ode_solver

    !     INDENT(IN)
    !         real(sp/dp),  dimension(:),   ::  ystart          ... initial conditions
    !                                                               N inital values known at the time x1 (given N ODEs)
    !
    !         real(sp/dp)                   ::  x1              ... initial time
    !
    !         real(sp/dp)                   ::  x2              ... final time
    !
    !         real(sp/dp)                   ::  h1              ... guessed first stepsize
    !
    !         interface                     ::  derivs_sp/dp    ... returns derivatives dydx of y at x

    !     INDENT(INOUT)
    !         none

    !     INDENT(OUT)
    !         none

    !     INDENT(IN), OPTIONAL
    !         character(*)  ::  fileName    ... name of the file in which all the results will be stored
    !                                           DEFAULT: no output files will be created
    !
    !         real(sp/dp)   ::  hminIN      ... minimum allowed stepsize (can be zero)
    !                                           DEFAULT: 0.0
    !
    !         real(sp/dp)   ::  epsIN       ... accuracy (overall tolerance level)
    !                                           DEFAULT: 10.0**(-6.0)
    !
    !     INDENT(INOUT), OPTIONAL
    !         none

    !     INDENT(OUT), OPTIONAL
    !         real(sp/dp)                                  ::  Tout     ... print out the elapsed CPU time
    !                                                                       DEFAULT: no outputs will be created
    !
    !         real(sp/dp),  dimension(:),   allocatable    ::  xout     ... storage for outputs (time)
    !                                                                       DEFAULT: no outputs will be created
    !
    !         real(sp/dp),  dimension(:,:), allocatable    ::  yout     ... storage for outputs (incremented variables)
    !                                                                       DEFAULT: no outputs will be created

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
    !         Modified, ....


    ! SINGLE PRECISION
    subroutine RKas_sp( ystart, x1, x2, h1, derivs, Tout, xout, yout, fileName, hminIN, epsIN )
        ! subroutine RKas_sp( ystart, x1, x2, h1, derivs, rkqs, Tout, xout, yout, fileName, hminIN, epsIN )

        use mo_kind, only : i4, sp
        use mo_nrutil, only : nrerror, reallocate

        implicit none

        real(sp), intent(in) :: x1, x2, h1
        real(sp), dimension(:), intent(in) :: ystart
        real(sp), intent(in), optional :: hminIN, epsIN
        character(*), intent(in), optional :: fileName
        real(sp), intent(out), optional :: Tout
        real(sp), dimension(:), allocatable, intent(out), optional :: xout
        real(sp), dimension(:,:), allocatable, intent(out), optional :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : sp
                implicit none
                real(sp), intent(in) :: x
                real(sp), dimension(:), intent(in) :: y
                real(sp), dimension(:), intent(out) :: dydx
            end subroutine derivs

        !            subroutine rkqs( y, dydx, x, htry, eps, yscal, hdid, hnext, derivs )
        !                use mo_kind, only : sp
        !                real(sp), dimension(:), intent(inout) :: y
        !                real(sp), dimension(:), intent(in) :: dydx, yscal
        !                real(sp), intent(inout) :: x
        !                real(sp), intent(in) :: htry, eps
        !                real(sp), intent(out) :: hdid, hnext
        !                interface
        !                    subroutine derivs( x, y, dydx )
        !                        use mo_kind, only : sp
        !                        implicit none
        !                        real(sp), intent(in) :: x
        !                        real(sp), dimension(:), intent(in) :: y
        !                        real(sp), dimension(:), intent(out) :: dydx
        !                    end subroutine derivs
        !                end interface
        !            end subroutine rkqs

        end interface

        integer(i4) :: nstp, k
        real(sp) :: h, hdid, hnext, x, hmin, eps, T1, T2, xsav, dxsav
        real(sp), dimension(size(ystart)) :: dydx, y, yscal

        integer(i4), parameter :: MAXstp = 10000_i4

        if( present(hminIN) ) then
            hmin = hminIN
        else
            hmin = 0.0_sp
        end if

        if( present(epsIN) ) then
            eps = epsIN
        else
            eps = 10._sp**(-6._sp)
        end if

        call header( x1, x2, h1, ystart, 'RUNGE-KUTTA with adaptive stepsize control', hmin, eps )   ! header with settings
                                                                                                     ! informations
        if( present(xout) .or. present(yout) .or. present(fileName) ) save_steps = .true.

        x = x1
        h = sign( h1, x2-x1 )
        nok = 0_i4
        nbad = 0_i4
        kount = 0_i4
        y(:) = ystart(:)
        T1 = 0._sp
        T2 = T1
        xsav = tiny(1._sp)
        dxsav = tiny(1._sp)

        if ( save_steps ) then
            xsav = x-2.0_sp*dxsav                    ! assures storage of first step
            nullify( xp_sp, yp_sp )
            allocate( xp_sp(256) )
            allocate( yp_sp(size(ystart), size(xp_sp)) )
        end if

        if ( save_steps ) call save_a_step              ! save initial step

        if( present(Tout) ) call cpu_time(T1)
        do nstp=1, MAXstp                               ! take at most MAXstp steps
            call derivs( x, y, dydx )
            yscal(:) = abs( y(:) ) + abs( h*dydx(:) ) + tiny(1._sp)     ! scaling used to monitor accuracy --> CAN BE MODIFIED...

            if ( save_steps .and. (abs(x-xsav) > abs(dxsav)) ) &     ! store intermediate results
                call save_a_step

            if ( (x+h-x2)*(x+h-x1) > 0.0_sp ) h = x2-x                  ! if stepsize can overshoot, decrease
            call RKqs( y, dydx, x, h, eps, yscal, hdid, hnext, derivs ) ! <-- if you pass it as argument
                                                                        !     you have to call sp or dp!!!

            if ( hdid == h ) then
                nok = nok+1
            else
                nbad = nbad+1
            end if

            if ( (x-x2)*(x2-x1) >= 0.0_sp ) then        ! are we done?!?!
                if( present(Tout) ) then
                    call cpu_time(T2)
                    Tout = T2 - T1
                    write(*,*) 'Elapsed CPU time RUNGE KUTTA as = ', Tout
                end if
                if ( save_steps ) then
                    call save_a_step                    ! save final step

                    if( present(fileName) ) then
                        open( unit = 10, file = fileName, action = "write", status = "replace" )  ! dumping results
                        do k=1, kount
                            write(10,*) xp_sp(k), yp_sp(:,k)
                        end do
                        close( unit = 10 )
                        write(*,*) ''
                        write(*,*) fileName, ' was successfully saved!'
                    end if

                    if( present(xout) ) then
                        allocate( xout(kount) )                 ! allocate storage for outputs
                        xout(:) = xp_sp(1:kount)
                    end if
                    if( present(yout) ) then
                        allocate( yout(size(yp_sp,1),kount) )   ! allocate storage for outputs
                        yout(:,:) = yp_sp(:,1:kount)
                    end if

                    deallocate( xp_sp, yp_sp )                  ! clear out old stored variables
                end if

                return                                          ! normal exit
            end if

            if ( abs(hnext) < hmin ) &
                call nrerror( 'WTF! ...stepsize smaller than minimum in RKas!!!' )
            h = hnext
        end do

        call nrerror( 'dude, too many steps in RKas!!!' )

    contains

        subroutine save_a_step      ! --> like a macro in C
            kount = kount+1
            if ( kount > size(xp_sp) ) then
                xp_sp=>reallocate( xp_sp, 2*size(xp_sp) )
                yp_sp=>reallocate( yp_sp, size(yp_sp,1), size(xp_sp) )
            end if
            xp_sp(kount) = x
            yp_sp(:, kount) = y(:)
            xsav = x
        end subroutine save_a_step

    end subroutine RKas_sp


    subroutine RKas_dp( ystart, x1, x2, h1, derivs, Tout, xout, yout, fileName, hminIN, epsIN )     ! DOUBLE PRECISION
        ! subroutine RKas_dp( ystart, x1, x2, h1, derivs, rkqs, Tout, xout, yout, fileName, hminIN, epsIN )

        use mo_kind, only : i4, dp
        use mo_nrutil, only : nrerror, reallocate

        implicit none

        real(dp), intent(in) :: x1, x2, h1
        real(dp), intent(in), optional :: epsIN, hminIN
        character(*), intent(in), optional :: fileName
        real(dp), intent(out), optional :: Tout
        real(dp), dimension(:), intent(in) :: ystart
        real(dp), dimension(:), allocatable, intent(out), optional :: xout
        real(dp), dimension(:,:), allocatable, intent(out), optional :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : dp
                implicit none
                real(dp), intent(in) :: x
                real(dp), dimension(:), intent(in) :: y
                real(dp), dimension(:), intent(out) :: dydx
            end subroutine derivs

        !		subroutine rkqs( y, dydx, x, htry, eps, yscal, hdid, hnext, derivs )
        !			use mo_kind, only : dp
        !			real(dp), dimension(:), intent(inout) :: y
        !			real(dp), dimension(:), intent(in) :: dydx, yscal
        !			real(dp), intent(inout) :: x
        !			real(dp), intent(in) :: htry, eps
        !			real(dp), intent(out) :: hdid, hnext
        !		interface
        !			subroutine derivs( x, y, dydx )
        !			use mo_kind, only : dp
        !			implicit none
        !			real(dp), intent(in) :: x
        !			real(dp), dimension(:), intent(in) :: y
        !			real(dp), dimension(:), intent(out) :: dydx
        !			end subroutine derivs
        !		end interface
        !		end subroutine rkqs

        end interface

        integer(i4) :: nstp, k
        real(dp) :: h, hdid, hnext, x, hmin, eps, T1, T2, xsav, dxsav
        real(dp), dimension(size(ystart)) :: dydx, y, yscal

        integer(i4), parameter :: MAXstp = 10000_i4

        if( present(hminIN) ) then
            hmin = hminIN
        else
            hmin = 0.0_dp
        end if

        if( present(epsIN) ) then
            eps = epsIN
        else
            eps = 10._dp**(-6._dp)
        end if

        call header( x1, x2, h1, ystart, 'RUNGE-KUTTA with adaptive stepsize control', hmin, eps )   ! header with settings
                                                                                                     ! informations
        if( present(xout) .or. present(yout) .or. present(fileName) ) save_steps = .true.

        x = x1
        h = sign( h1, x2-x1 )
        nok = 0_i4
        nbad = 0_i4
        kount = 0_i4
        y(:) = ystart(:)
        T1 = 0._dp
        T2 = T1
        xsav = tiny(1._dp)
        dxsav = tiny(1._dp)

        if ( save_steps ) then
            xsav = x-2.0_dp*dxsav                    ! assures storage of first step
            nullify( xp_dp, yp_dp )
            allocate( xp_dp(256) )
            allocate( yp_dp(size(ystart), size(xp_dp)) )
        end if

        if ( save_steps ) call save_a_step              ! save initial step

        if( present(Tout) ) call cpu_time(T1)
        do nstp=1, MAXstp                               ! take at most MAXstp steps
            call derivs( x, y, dydx )
            yscal(:) = abs( y(:) ) + abs( h*dydx(:) ) + tiny(1._dp)     ! scaling used to monitor accuracy --> CAN BE MODIFIED...

            if ( save_steps .and. (abs(x-xsav) > abs(dxsav)) ) &     ! store intermediate results
                call save_a_step

            if ( (x+h-x2)*(x+h-x1) > 0.0_dp ) h = x2-x                  ! if stepsize can overshoot, decrease
            call RKqs( y, dydx, x, h, eps, yscal, hdid, hnext, derivs ) ! <-- if you pass it as argument
                                                                        !     you have to call sp or dp!!!

            if ( hdid == h ) then
                nok = nok+1
            else
                nbad = nbad+1
            end if

            if ( (x-x2)*(x2-x1) >= 0.0_dp ) then        ! are we done?!?!
                if( present(Tout) ) then
                    call cpu_time(T2)
                    Tout = T2 - T1
                    write(*,*) 'Elapsed CPU time RUNGE KUTTA as = ', Tout
                end if
                if ( save_steps ) then
                    call save_a_step                    ! save final step

                    if( present(fileName) ) then
                        open( unit = 10, file = fileName, action = "write", status = "replace" )  ! dumping results
                        do k=1, kount
                            write(10,*) xp_dp(k), yp_dp(:,k)
                        end do
                        close( unit = 10 )
                        write(*,*) ''
                        write(*,*) fileName, ' was successfully saved!'
                    end if

                    if( present(xout) ) then
                        allocate( xout(kount) )                 ! allocate storage for outputs
                        xout(:) = xp_dp(1:kount)
                    end if
                    if( present(yout) ) then
                        allocate( yout(size(yp_dp,1),kount) )   ! allocate storage for outputs
                        yout(:,:) = yp_dp(:,1:kount)
                    end if

                    deallocate( xp_dp, yp_dp )                  ! clear out old stored variables
                end if

                return                                          ! normal exit
            end if

            if ( abs(hnext) < hmin ) &
                call nrerror( 'WTF! ...stepsize smaller than minimum in RKas!!!' )
            h = hnext
        end do

        call nrerror( 'dude, too many steps in RKas!!!' )

    contains

        subroutine save_a_step      ! --> like a macro in C
            kount = kount+1
            if ( kount > size(xp_dp) ) then
                xp_dp=>reallocate( xp_dp, 2*size(xp_dp) )
                yp_dp=>reallocate( yp_dp, size(yp_dp,1), size(xp_dp) )
            end if
            xp_dp(kount) = x
            yp_dp(:, kount) = y(:)
            xsav = x
        end subroutine save_a_step

    end subroutine RKas_dp


    ! ============================================================================
    !	PRIVATE METHODS
    ! ============================================================================

                ! EULER		SINGLE PRECISION
    subroutine EULER_sp( y, dydx, x, h, yout , derivs )

        use mo_kind, only : i4, sp
        use mo_nrutil, only : assert_eq

        implicit none
        real(sp), dimension(:), intent(in) :: y, dydx
        real(sp), intent(in) :: x, h
        real(sp), dimension(:), intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : sp
                implicit none
                real(sp), intent(in) :: x
                real(sp), dimension(:), intent(in) :: y
                real(sp), dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        integer(i4) :: ndum
        real(sp), dimension(:), allocatable :: dyt, yt

        ndum = assert_eq( size(y), size(dydx), size(yout), 'EULER_sp' )

        allocate( dyt(ndum), yt(ndum) )

        yt = y
        call derivs( x, yt, dyt )
        yout = yt + h*dyt

    end subroutine EULER_sp


                ! EULER		DOUBLE PRECISION
    subroutine EULER_dp( y, dydx, x, h, yout , derivs )

        use mo_kind, only : i4, dp
        use mo_nrutil, only : assert_eq

        implicit none
        real(dp), dimension(:), intent(in) :: y, dydx
        real(dp), intent(in) :: x, h
        real(dp), dimension(:), intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : dp
                implicit none
                real(dp), intent(in) :: x
                real(dp), dimension(:), intent(in) :: y
                real(dp), dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        integer(i4) :: ndum
        real(dp), dimension(:), allocatable :: dyt, yt

        ndum = assert_eq( size(y), size(dydx), size(yout), 'EULER_dp' )

        allocate( dyt(ndum), yt(ndum) )

        yt = y
        call derivs( x, yt, dyt )
        yout = yt + h*dyt

    end subroutine EULER_dp


                ! RUNGE-KUTTA 4	SINGLE PRECISION
    subroutine RK4_sp( y, dydx, x, h, yout , derivs )
        ! Given values for the N variables y and their derivatives dydx known at x, use the fourth-order
        ! Runge-Kutta method to advance the solution over an interval h and return the incremented
        ! variables as yout, which need not be a distinct array from y.
        ! y, dydx and yout are all of length N.

        use mo_kind, only : i4, sp
        use mo_nrutil, only : assert_eq

        implicit none
        real(sp), dimension(:), intent(in) :: y, dydx
        real(sp), intent(in) :: x, h
        real(sp), dimension(:), intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : sp
                implicit none
                real(sp), intent(in) :: x
                real(sp), dimension(:), intent(in) :: y
                real(sp), dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        integer(i4) :: ndum
        real(sp) :: h6, hh, xh
        real(sp), dimension(:), allocatable :: dym, dyt, yt

        ndum = assert_eq( size(y), size(dydx), size(yout), 'RK4_sp' )

        allocate( dym(ndum), dyt(ndum), yt(ndum) )

        hh = h*0.5_sp
        h6 = h/6.0_sp
        xh = x+hh

        yt = y + hh*dydx                            ! first step
        call derivs( xh, yt, dyt )                  ! second step
        yt = y + hh*dyt
        call derivs( xh, yt, dym )                  ! third step
        yt = y + h*dym
        dym = dyt + dym
        call derivs( x+h, yt, dyt )                 ! fourth step
        yout = y + h6*( dydx+dyt+2.0_sp*dym )       ! accumulate increments with proper weights

    end subroutine RK4_sp


                ! RUNGE-KUTTA 4	DOUBLE PRECISION
    subroutine RK4_dp( y, dydx, x, h, yout , derivs )

        use mo_kind, only : i4, dp
        use mo_nrutil, only : assert_eq

        implicit none
        real(dp), dimension(:), intent(in) :: y, dydx
        real(dp), intent(in) :: x, h
        real(dp), dimension(:), intent(out) :: yout

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : dp
                implicit none
                real(dp), intent(in) :: x
                real(dp), dimension(:), intent(in) :: y
                real(dp), dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        integer(i4) :: ndum
        real(dp) :: h6, hh, xh
        real(dp), dimension(:), allocatable :: dym, dyt, yt

        ndum = assert_eq( size(y), size(dydx), size(yout), 'RK4_dp' )

        allocate( dym(ndum), dyt(ndum), yt(ndum) )

        hh = h*0.5_dp
        h6 = h/6.0_dp
        xh = x+hh

        yt = y + hh*dydx                            ! first step
        call derivs( xh, yt, dyt )                  ! second step
        yt = y + hh*dyt
        call derivs( xh, yt, dym )                  ! third step
        yt = y + h*dym
        dym = dyt + dym
        call derivs( x+h, yt, dyt )                 ! fourth step
        yout = y + h6*( dydx+dyt+2.0_dp*dym )       ! accumulate increments with proper weights

    end subroutine RK4_dp


                ! RUNGE-KUTTA Fifth order  (adjust stepsize)        SINGLE PRECISION
    subroutine RKqs_sp( y, dydx, x, htry, eps, yscal, hdid, hnext, derivs )
        ! Fifth order Runge-Kutta step with monitoring of local truncation error to ensure accuracy
        ! and adjust stepsize. Input are the dependent variable vector y and its derivative dydx at
        ! the starting value of the independent variable x. Also input are the stepsize to be attempted
        ! htry, the required accuracy eps, and the vector yscal against which the error is scaled.
        ! y, dydx, and yscal are all of the same length. On output, y and x are replaced by their new values,
        ! hdid is the stepsize that was actually accomplished, and hnext is the estimated next stepsize.

        use mo_kind, only : i4, sp
        use mo_nrutil, only : assert_eq, nrerror

        implicit none

        real(sp), dimension(:), intent(in) :: dydx, yscal
        real(sp), intent(in) :: htry, eps
        real(sp), dimension(:), intent(inout) :: y
        real(sp), intent(inout) :: x
        real(sp), intent(out) :: hdid, hnext

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : sp
                implicit none
                real(sp), intent(in) :: x
                real(sp), dimension(:), intent(in) :: y
                real(sp), dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        integer(i4) :: ndum
        real(sp) :: errmax, h, htemp, xnew
        real(sp), dimension(:), allocatable :: yerr, ytemp

        real(sp), parameter :: safety=0.9_sp, pgrow=-0.2_sp, pshrnk=-0.25_sp, &
            errcon = (5._sp/safety)**(1._sp/pgrow)

        ndum = assert_eq( size(y), size(dydx), size(yscal), 'rkqs' )

        allocate( yerr(ndum), ytemp(ndum) )

        h = htry            ! set stepsize to the initial trial value

        do
            call RKck( y, dydx, x, h, ytemp, yerr, derivs )     ! take a step
            errmax = maxval( abs(yerr(:)/yscal(:)) )/eps        ! evaluate accuracy
            if ( errmax <= 1.0_sp ) exit                        ! step succeeded
            htemp = safety*h*(errmax**pshrnk)                   ! truncation error too large, reduce stepsize
            h = sign( max( abs(htemp), 0.1_sp*abs(h) ), h )     ! no more than a factor of 10
            xnew = x+h
            if ( xnew == x ) call nrerror( 'hey!!! stepsize underflow in RKqs' )
        end do                              ! go back for another try

        if ( errmax > errcon ) then         ! compute size of next step
            hnext = safety*h*(errmax**pgrow)
        else                                ! no more than a factor of 5 increase
            hnext = 5.0_sp*h
        end if

        hdid = h
        x = x+h
        y(:) = ytemp(:)

    end subroutine RKqs_sp


                ! RUNGE-KUTTA Fifth order (adjust stepsize)		DOUBLE PRECISION
    subroutine RKqs_dp( y, dydx, x, htry, eps, yscal, hdid, hnext, derivs )

        use mo_kind, only : i4, dp
        use mo_nrutil, only : assert_eq, nrerror

        implicit none

        real(dp), dimension(:), intent(in) :: dydx, yscal
        real(dp), intent(in) :: htry, eps
        real(dp), dimension(:), intent(inout) :: y
        real(dp), intent(inout) :: x
        real(dp), intent(out) :: hdid, hnext

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : dp
                implicit none
                real(dp), intent(in) :: x
                real(dp), dimension(:), intent(in) :: y
                real(dp), dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        integer(i4) :: ndum
        real(dp) :: errmax, h, htemp, xnew
        real(dp), dimension(:), allocatable :: yerr, ytemp

        real(dp), parameter :: safety=0.9_dp, pgrow=-0.2_dp, pshrnk=-0.25_dp, &
            errcon = (5._dp/safety)**(1._dp/pgrow)

        ndum = assert_eq( size(y), size(dydx), size(yscal), 'rkqs' )

        allocate( yerr(ndum), ytemp(ndum) )

        h = htry            ! set stepsize to the initial trial value

        do
            call RKck( y, dydx, x, h, ytemp, yerr, derivs )     ! take a step
            errmax = maxval( abs(yerr(:)/yscal(:)) )/eps        ! evaluate accuracy
            if ( errmax <= 1.0_dp ) exit                        ! step succeeded
            htemp = safety*h*(errmax**pshrnk)                   ! truncation error too large, reduce stepsize
            h = sign( max( abs(htemp), 0.1_dp*abs(h) ), h )     ! no more than a factor of 10
            xnew = x+h
            if ( xnew == x ) call nrerror( 'hey!!! stepsize underflow in RKqs' )
        end do                              ! go back for another try

        if ( errmax > errcon ) then         ! compute size of next step
            hnext = safety*h*(errmax**pgrow)
        else                                ! no more than a factor of 5 increase
            hnext = 5.0_dp*h
        end if

        hdid = h
        x = x+h
        y(:) = ytemp(:)

    end subroutine RKqs_dp


                ! CASH-KARP RUNGE-KUTTA step	SINGLE PRECISION
    subroutine RKck_sp( y, dydx, x, h, yout, yerr, derivs )
        ! Given values for N variables y and their derivatives dydx known at x, use the fifth-order
        ! Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
        ! the incremented variables as yout. Also return an estimate of the local truncation error
        ! in yout using the embedded fourth order method.

        use mo_kind, only : i4, sp
        use mo_nrutil, only : assert_eq

        implicit none

        real(sp), dimension(:), intent(in) :: y, dydx
        real(sp), intent(in) :: x, h
        real(sp), dimension(:), intent(out) :: yout, yerr

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : sp
                implicit none
                real(sp), intent(in) :: x
                real(sp), dimension(:), intent(in) :: y
                real(sp), dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        integer(i4) :: ndum
        real(sp), dimension(:), allocatable :: ak2, ak3, ak4, ak5, ak6, ytemp

        real(sp), parameter :: A2=.2_sp, A3=.3_sp, A4=.6_sp, A5=1._sp, A6=.875_sp, B21=.2_sp, B31=3._sp/40._sp, &
            B32=9._sp/40._sp, B41=.3_sp, B42=-.9_sp, B43=1.2_sp, B51=-11._sp/54._sp, B52=2.5_sp, &
            B53=-70._sp/27._sp, B54=35._sp/27._sp, B61=1631._sp/55296._sp, B62=175._sp/512._sp, &
            B63=575._sp/13824._sp, B64=44275._sp/110592._sp, B65=253._sp/4096._sp, &
            C1=37._sp/378._sp, C3=250._sp/621._sp, C4=125._sp/594._sp, C6=512._sp/1771._sp, &
            DC1=C1-2825._sp/27648._sp, DC3=C3-18575._sp/48384._sp, &
            DC4=C4-13525._sp/55296._sp, DC5=-277._sp/14336._sp, DC6=C6-.25_sp

        ndum = assert_eq( size(y), size(dydx) ,size(yout) ,size(yerr) , 'rkck' )

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

    end subroutine RKck_sp


                ! CASH-KARP RUNGE-KUTTA step	DOUBLE PRECISION
    subroutine RKck_dp( y, dydx, x, h, yout, yerr, derivs )

        use mo_kind, only : i4, dp
        use mo_nrutil, only : assert_eq

        implicit none

        real(dp), dimension(:), intent(in) :: y, dydx
        real(dp), intent(in) :: x, h
        real(dp), dimension(:), intent(out) :: yout, yerr

        interface
            subroutine derivs( x, y, dydx )
                use mo_kind, only : dp
                implicit none
                real(dp), intent(in) :: x
                real(dp), dimension(:), intent(in) :: y
                real(dp), dimension(:), intent(out) :: dydx
            end subroutine derivs
        end interface

        integer(i4) :: ndum
        real(dp), dimension(:), allocatable :: ak2, ak3, ak4, ak5, ak6, ytemp

        real(dp), parameter :: A2=.2_dp, A3=.3_dp, A4=.6_dp, A5=1._dp, A6=.875_dp, B21=.2_dp, B31=3._dp/40._dp, &
            B32=9._dp/40._dp, B41=.3_dp, B42=-.9_dp, B43=1.2_dp, B51=-11._dp/54._dp, B52=2.5_dp, &
            B53=-70._dp/27._dp, B54=35._dp/27._dp, B61=1631._dp/55296._dp, B62=175._dp/512._dp, &
            B63=575._dp/13824._dp, B64=44275._dp/110592._dp, B65=253._dp/4096._dp, &
            C1=37._dp/378._dp, C3=250._dp/621._dp, C4=125._dp/594._dp, C6=512._dp/1771._dp, &
            DC1=C1-2825._dp/27648._dp, DC3=C3-18575._dp/48384._dp, &
            DC4=C4-13525._dp/55296._dp, DC5=-277._dp/14336._dp, DC6=C6-.25_dp

        ndum = assert_eq( size(y), size(dydx) ,size(yout) ,size(yerr) , 'rkck' )

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

    end subroutine RKck_dp

    subroutine header_sp( x1, x2, h, xi, MethodName, hmin, eps )

        real(sp), intent(in) :: x1, x2, h
        real(sp), dimension(:), intent(in) :: xi
        character(*), intent(in) :: MethodName
        real(sp), intent(in), optional :: hmin, eps

        write(*,*) ''
        write(*,*) '//////////////////////////////////////////////////////////////////////'
        write(*,*) '//////////////////          ODE SOLVER           /////////////////////'
        write(*,*) '//////////////////////////////////////////////////////////////////////'
        write(*,*) '                 --> ', MethodName
        write(*,*) ''
        write(*,*) '--------------SETTINGS: START--------------'
        write(*,*) 'number of equations nEq = ', size(xi)
        write(*,*) 'intial time x1 = ', x1
        write(*,*) 'final time x2 = ', x2
        write(*,*) 'incremental time step h = ', h
        write(*,*) 'initial conditions xi = ', xi(:)
        if( present(hmin) ) write(*,*) 'minimum allowed stepsize hmin = ', hmin
        if( present(eps) ) write(*,*) 'accuracy eps = ', eps
        write(*,*) '--------------SETTINGS: END----------------'
        write(*,*) ''
    end subroutine header_sp

    subroutine header_dp( x1, x2, h, xi, MethodName, hmin, eps )

        real(dp), intent(in) :: x1, x2, h
        real(dp), dimension(:), intent(in) :: xi
        character(*), intent(in) :: MethodName
        real(dp), intent(in), optional :: hmin, eps

        write(*,*) ''
        write(*,*) '//////////////////////////////////////////////////////////////////////'
        write(*,*) '//////////////////          ODE SOLVER           /////////////////////'
        write(*,*) '//////////////////////////////////////////////////////////////////////'
        write(*,*) '                  -->', MethodName
        write(*,*) ''
        write(*,*) '--------------SETTINGS: START--------------'
        write(*,*) 'number of equations nEq = ', size(xi)
        write(*,*) 'intial time x1 = ', x1
        write(*,*) 'final time x2 = ', x2
        write(*,*) 'incremental time step h = ', h
        write(*,*) 'initial conditions xi = ', xi(:)
        if( present(hmin) ) write(*,*) 'minimum allowed stepsize hmin = ', hmin
        if( present(eps) ) write(*,*) 'accuracy eps = ', eps
        write(*,*) '--------------SETTINGS: END----------------'
        write(*,*) ''
    end subroutine header_dp

end module mo_ode_solver
