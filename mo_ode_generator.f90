!> \file mo_ode_generator.f90

!> \brief Generating and solving systems of ODE.

!> \details Given N reactants generates & solves all the corresponding ODE system.

!> \authors Giovanni Dalmasso
!> \date Jul 2013

module mo_ode_generator

    ! This module generates & solves all the corresponding ODE system of N given reactants.

    ! Written Giovanni Dalmasso, Jul 2013

    ! License
    ! -------
    ! This file is part of the UFZ Fortran library.

    ! The UFZ Fortran library is free software: you can redistribute it and/or modify
    ! it under the terms of the GNU Lesser General Public License as published by
    ! the Free Software Foundation, either version 3 of the License, or
    ! (at your option) any later version.

    ! The UFZ Fortran library is distributed in the hope that it will be useful,
    ! but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    ! GNU Lesser General Public License for more details.

    ! You should have received a copy of the GNU Lesser General Public License
    ! along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
    ! If not, see <http://www.gnu.org/licenses/>.

    ! Copyright 2013 Giovanni Dalmasso

    use mo_kind, only : dp

    implicit none

    public ::  ode_generator    ! Generating and solving systems of linear ODEs.

    ! ------------------------------------------------------------------

    !     NAME
    !         ode_generator

    !     PURPOSE
    !>        \brief Generating and solving systems of ODEs.
    !
    !>        \details  Given N reactants generates & solves all the corresponding ODE system.
    !
    !     INTENT(IN)
    !>        \param[in] "integer(i4)  :: n_reac"               number of reactants involved in the network
    !>        \param[in] "real(dp),    :: dt"                   initial size time step
    !>        \param[in] "real(dp),    :: t_start"              initial time
    !>        \param[in] "real(dp),    :: t_end"                final time
    !>        \param[in] "real(dp),    :: x_initial(:)"         initial condition for each reactant
    !>        \param[in] "real(dp),    :: para(:,:)"            array of parameters: \n
    !>                                                          par(j,k) is the reaction rate to convert
    !>                                                          reactant j into reactant k (not all may be used)\n
    !>                                                          par(j,j) = feedback loop\n
    !>                                                          if reaction is not existing par(j,k) = 0.0\n
    !>                                                          reaction rate has to be positive
    !
    !     INTENT(INOUT)
    !         None
    !
    !     INTENT(OUT)
    !>        \param[out] "real(dp),   :: ode_results(:,:,:)"   ode_results(i,j,k) concentration of reactant (k-1)
    !>                                                          at the jth time point in network permutation i\n
    !>                                                          jth time point is ode_results(i,j,1)\n
    !>                                                          dim_1 = number of permutations\n
    !>                                                          dim_2 = number of time points\n
    !>                                                          dim_3 = number of reactants + 1 \n
    !>                                                                  (dim_3=1 : time points)\n
    !>                                                                  (dim_3=i : reactant i-1)\n
    !
    !     INTENT(IN), OPTIONAL
    !>        \param[in] "real(dp),    optional :: nodata_value"choosing the value to no-data\n
    !>                                                          DEFAULT: -9999.9_dp
    !>        \param[in] "integer(i8), optional :: jPerm"       choosing one single permutation\n
    !>                                                          ode_results(1,:,:) is returned\n
    !>                                                          DEFAULT: all
    !>        \param[in] "logical,     optional :: printflag"   flag for printing results on screen\n
    !>                                                          DEFAULT: false
    !
    !     INTENT(INOUT), OPTIONAL
    !         None
    !
    !     INTENT(OUT), OPTIONAL
    !         None
    !
    !     RESTRICTIONS
    !>       \note No sink/ source terms. This means a reactant can not be produced or degraded.\n
    !>             No prior information can be passed. This means no reaction can be fixed or discarded
    !>             in each of the generated networks from the beginning.
    !>             Meaning that each generated network has to contain a specific reation or not.\n
    !>             Number of reactants less or equal then 7. This is due to the fact that the number of
    !>             permutations will be too large to be stored in integer i8.
    !
    !     EXAMPLE
    !         -> see example in test directory
    !
    !     LITERATURE
    !
    !     HISTORY
    !>        \author Giovanni Dalmasso
    !>        \date Jul 2013
    !         Modified Giovanni Dalmasso, Apr 2014 - added printflag
    !                                              - added variable to store results
    !                                              - renamed variables
    ! ------------------------------------------------------------------
    interface ode_generator
        module procedure ode_generator_dp
    end interface

    ! Private method
    interface get_odes
        module procedure get_odes_dp
    end interface

    interface set_alpha
        module procedure set_alpha_dp
    end interface

    interface set_beta
        module procedure set_beta_dp
    end interface

    interface get_binary
        module procedure get_binary_i4, get_binary_i8
    end interface

    interface get_binaryPer
        module procedure get_binaryPer_dp
    end interface

    private

    !! module variables
    real(dp),   dimension(:,:), allocatable ::  alpha   ! on/off variables matrix
    real(dp),   dimension(:,:), allocatable ::  beta    ! parameters matrix

contains

    subroutine ode_generator_dp( n_reac, dt, t_start, t_end, x_initial, para,   &   ! IN
            ode_results,                                                        &   ! OUT
            nodata_value, jPerm, printflag                                      &   ! optional IN
            )

        use mo_kind,            only    :   i4, i8, dp
        use mo_ode_solver,      only    :   RK4as

        implicit none

        !! Intent IN
        integer(i4),                                    intent(in)  ::  n_reac      ! number of reactants involved in the network
        real(dp),                                       intent(in)  ::  dt          ! initial size time step
        real(dp),                                       intent(in)  ::  t_start     ! initial time
        real(dp),                                       intent(in)  ::  t_end       ! final time
        real(dp),   dimension(:),                       intent(in)  ::  x_initial   ! initial conditions
        real(dp),   dimension(:,:),                     intent(in)  ::  para        ! array of paramenters

        !! Intent IN optional
        real(dp),       optional,                       intent(in)  ::  nodata_value! choosing the value to no-data
        integer(i8),    optional,                       intent(in)  ::  jPerm       ! choosing one single permutation
        logical,        optional,                       intent(in)  ::  printflag   ! flag for printing results on screen

        !! Intent OUT
        real(dp),   dimension(:,:,:),   allocatable,    intent(out) ::  ode_results
        !                                                                   ! ode_results(i,j,k) concentration of reactant (k-1)
        !                                                                   ! at the jth time point in network permutation i\n
        !                                                                   ! jth time point is ode_results(i,j,1)\n
        !                                                                   ! dim_1 = number of permutations\n
        !                                                                   ! dim_2 = number of time points\n
        !                                                                   ! dim_3 = number of reactants + 1 \n
        !                                                                   !         (dim_3=1 : time points)\n
        !                                                                   !         (dim_3=i : reactant i-1)

        !! Internal variables
        integer(i4)                                 ::  ii, jj            ! counters
        integer(i4)                                 ::  dim_new, dim_old  ! dimensions needed for re-allocation
        integer(i8)                                 ::  iPer              ! current permutation
        integer(i8)                                 ::  nPer              ! number of permutations
        real(dp)                                    ::  nodata_valueIn    ! no-data value
        real(dp)                                    ::  dt_min            ! step size
        real(dp),   dimension(:),       allocatable ::  xout              ! time
        real(dp),   dimension(:,:),     allocatable ::  yout              ! concentration per reactant
        real(dp),   dimension(:,:,:),   allocatable ::  results_temp      ! temp variable
        logical                                     ::  printflagIn

        !! Parameters
        real(dp),   parameter       :: eps = 1e-9_dp            ! epsilon for the ODE solver

        dt_min = dt*1e-11_dp                                    ! minimum allowed stepsize
        nPer = 2_i8**(int(n_reac*n_reac,i4))                    ! number of permutations

        ! check reaction rates
        if ( any(para .lt. 0.0_dp) ) then
            write(*,*) 'mo_ode_generator: Reaction rates need to be non-negative!'
            stop
        end if

        ! seeting no-data value
        nodata_valueIn = -9999.9_dp
        if ( present(nodata_value) )   nodata_valueIn = nodata_value

        printflagIn = .false.
        if ( present(printflag) )   printflagIn = printflag

        !! set up Betas in the ode generator
        call set_beta( para )

        if ( present(jPerm) )    then
            !! ONLY ONE PERMUTATION (jPerm)

            if( jPerm .gt. nPer ) stop 'mo_ode_generator: ode_generator_dp: maximum number of permutation exceeded!!!'
            if( jPerm .lt. 1_i8 ) stop 'mo_ode_generator: ode_generator_dp: jPerm should be greater than 0!!!'

            call set_alpha( n_reac, jPerm - 1_i8 )
            call RK4as( x_initial, t_start, t_end, dt, get_odes_dp, xout, yout, hmin=dt_min, eps=eps )

            !! filling the results in ode_results
            if( allocated(ode_results) )   deallocate(ode_results)
            allocate( ode_results( 1_4, size(xout), n_reac+1_i4 ) )

            ode_results(1,:,1) = xout
            forall ( ii = 2:n_reac+1_i4 )    ode_results(1,:,ii) = yout(:,ii-1_i4)

            !! writing the solutions of the current system of ODEs
            if ( printflagIn )  then
                write(*,*) ''
                write(*,*) 'permutation num -->', jPerm
                write(*,*) '        t                   y(t,j)    ( j=1 ,...,', n_reac,')'
                do jj=1, size(xout)
                    write(*,*) xout(jj), yout(jj,:)
                end do
                write(*,*) ''
            end if

        else
            !! LOOPING OVER ALL THE POSSIBLE PERMUTATIONS

            dim_old = 0_i4

            !! looping over all the permutations
            do iPer=0_i8, nPer-1_i8

                call set_alpha( n_reac, iPer )
                call RK4as( x_initial, t_start, t_end, dt, get_odes_dp, xout, yout, hmin=dt_min, eps=eps )

                dim_new = size(xout)

                if ( dim_old .lt. dim_new )   then
                    dim_old = dim_new

                    if ( .not.( allocated(ode_results) ) )   then

                        !! 1st time allocation
                        allocate( ode_results( nPer, dim_new, n_reac+1_i4 ) )
                        allocate( results_temp( nPer, dim_new, n_reac+1_i4 ) )

                        !! save results
                        ode_results(iPer+1_i8,:,1) = xout
                        forall ( jj = 2:n_reac+1_i4 )    ode_results(iPer+1_i8,:,jj) = yout(:,jj-1_i4)

                    else

                        if( allocated(results_temp) )   deallocate(results_temp)
                        allocate( results_temp( size(ode_results,1), size(ode_results,2), size(ode_results,3) ) )

                        !! save a copy of the old array
                        results_temp(:,:,:) = ode_results(:,:,:)

                        !! dealloc + realloc
                        if( allocated(ode_results) )   deallocate(ode_results) ! to be sure!
                        allocate( ode_results( nPer, dim_new, n_reac+1_i4 ) )

                        ode_results(:,:,:) = -9999.9_dp

                        !! copy old results back
                        ode_results(1:size(results_temp,1), 1:size(results_temp,2), 1:size(results_temp,3)) = results_temp(:,:,:)

                        !! save new results
                        ode_results(iPer+1_i8,:,1) = xout
                        forall ( jj = 2:n_reac+1_i4 )    ode_results(iPer+1_i8,:,jj) = yout(:,jj-1_i4)

                    end if

                else

                    !! save  results
                    ode_results(iPer+1_i8,1:dim_new,1) = xout
                    forall ( jj = 2:n_reac+1_i4 )    ode_results(iPer+1_i8,1:dim_new,jj) = yout(:,jj-1_i4)

                end if


                !! writing the solutions of the current system of ODEs
                if ( printflagIn )  then
                    write(*,*) ''
                    write(*,*) 'num permutation -->', iPer+1_i4
                    write(*,*) '        t                   y(t,j)    ( j=1 ,...,', n_reac,')'
                    do jj=1, size(xout)
                        write(*,*) xout(jj), yout(jj,:)
                    end do
                    write(*,*) ''
                end if

            end do

        end if

    end subroutine ode_generator_dp


    ! ============================================================================
    !   PRIVATE METHODS
    ! ============================================================================

    !! GENERATE ODES
    subroutine get_odes_dp( x, y, dydx )

        use mo_kind,    only    :   i4, dp
        use mo_nrutil,  only    :   assert_eq

        implicit none

        !! Intent IN
        real(dp),                   intent(in)  ::  x       ! space variable
        real(dp), dimension(:),     intent(in)  ::  y       ! variables of the equations

        !! Intent OUT
        real(dp),   dimension(:),   intent(out) ::  dydx    ! derivative of y

        !! Internal variables
        integer(i4)     ::  j, k    ! counters
        integer(i4)     ::  nVar    ! number of equations

        !! dummy line
        dydx(1) = x
        dydx(:) = huge(1.0_dp)

        !! number of equations & checkin dimension
        nVar = assert_eq( size(y), size(alpha,1), size(beta,1), 'get_odes_dp')

        !! generating odes...
        do j=1, nVar
            do k=1, nVar
                if ( k .eq. j ) cycle   ! skipping when k = j
                dydx(j) = alpha(k,j)*beta(k,j)*y(k) - alpha(j,k)*beta(j,k)*y(j) + alpha(j,j)*beta(j,j)*y(j)
            end do
        end do

    end subroutine get_odes_dp


    !! SET ALPHA on/off variables matrix
    subroutine set_alpha_dp( nVar, line )

        use mo_kind,            only    :   i4, i8

        implicit none

        !! Intent IN
        integer(i4),    intent(in)                  :: nVar     ! number of variables
        integer(i8),    intent(in)                  :: line     ! line of set_alpha_dp

        if (.not. allocated(alpha)) then
           allocate(alpha(nVar,nVar))
        else
           deallocate(alpha)
           allocate(alpha(nVar,nVar))
        endif
        alpha =  get_binaryPer(nVar, line)

    end subroutine set_alpha_dp


    !! SET BETA paramenters matrix
    subroutine set_beta_dp( para )

        use mo_kind,            only    :   dp

        implicit none

        !! Intent IN
        real(dp),   dimension(:,:),   intent(in)  :: para     ! matrix of parameters

        if( allocated(beta) )   deallocate(beta)
        allocate( beta(size(para,1),size(para,1)) )

        beta = para

    end subroutine set_beta_dp


    !! decimal into a binary
    function get_binary_i4( decimal, nDigits )

        use mo_kind,    only    : i4, dp

        implicit none

        !! Intent IN
        integer(i4),    intent(in)                  :: decimal          ! decimal number to convert
        integer(i4),    intent(in)                  :: nDigits          ! number of digits required

        !! OUTPUT
        integer(i4),    dimension(:),   allocatable :: get_binary_i4    ! binary number

        !! Internal variables
        integer(i4)                     :: j            ! counters
        integer(i4)                     :: decimalIn

        !! checking the nuber of digits required for the binary
        if ( decimal .ne. 0_i4 .and.    &
             nDigits .lt. 4_i4*( int(log(real(decimal,dp))/log(16.0_dp), i4) ) )    stop 'get_binary_i4 --> to less digits!!!'

        allocate( get_binary_i4(nDigits) )

        get_binary_i4 = 0_i4

        !! converting decimal into binary
        if ( decimal .eq. 0_i4 )    then
            get_binary_i4 = 0_i4
        else
            decimalIn = decimal
            do j=1, nDigits
                get_binary_i4(j) = mod(decimalIn,2)
                decimalIn = decimalIn/2_i4
            end do
            !! inverting array
            get_binary_i4(nDigits:1:-1) = get_binary_i4(1:nDigits)
        end if

    end function get_binary_i4


    !! decimal into a binary
    function get_binary_i8( decimal, nDigits )

        use mo_kind,    only    : i4, i8, dp

        implicit none

        !! Intent IN
        integer(i8),    intent(in)                  :: decimal          ! decimal number to convert
        integer(i4),    intent(in)                  :: nDigits          ! number of digits required

        !! OUTPUT
        integer(i8),    dimension(:),   allocatable :: get_binary_i8    ! binary number

        !! Internal variables
        integer(i8)                     :: j            ! counters
        integer(i8)                     :: decimalIn

        !! checking the nuber of digits required for the binary
        if ( decimal .ne. 0_i8 .and.    &
             nDigits .lt. 4_i4*( int(log(real(decimal,dp))/log(16.0_dp), i4) ) )    stop 'get_binary_i8 --> to less digits!!!'

        allocate( get_binary_i8(nDigits) )

        get_binary_i8 = 0_i8

        !! converting decimal into binary
        if ( decimal .eq. 0_i8 )    then
            get_binary_i8 = 0_i8
        else
            decimalIn = decimal
            do j=1_i8, nDigits
                get_binary_i8(j) = mod(decimalIn,2_i8)
                decimalIn = decimalIn/2_i8
            end do
            !! inverting array
            get_binary_i8(nDigits:1:-1) = get_binary_i8(1:nDigits)
        end if

    end function get_binary_i8


    !! binary permutations
    function get_binaryPer_dp( nVar, line )

        use mo_kind,    only    :   i4, i8, dp

        implicit none

        !! Intent IN
        integer(i4),    intent(in)               :: nVar          ! number of digits I need
        integer(i8),    intent(in)               :: line          ! which permutation I want --> decimal number to convert

        !! OUTPUT
        real(dp),   dimension(nVar, nVar)  :: get_binaryPer_dp    ! binary number

        !! Internal variables
        !real(dp), dimension(:), allocatable  ::  get_binTemp
        real(dp), dimension(nVar*nVar) ::  get_binTemp

        get_binTemp = real( get_binary(line, nVar*nVar), dp )
        get_binaryPer_dp(:,:) = reshape( get_binTemp(:), (/ nVar, nVar /) )

    end function get_binaryPer_dp

end module mo_ode_generator
