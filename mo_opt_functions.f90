MODULE mo_opt_functions

  ! This modules provides test functions for minimisation routines

  ! Written, Jul 2012 Matthias Cuntz

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
  ! along with the UFZ Fortran library. If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2012 Matthias Cuntz

  use mo_kind, only: i4, dp

  IMPLICIT NONE

  PRIVATE

  ! ------------------------------------------------------------------
  ! test_min package of John Burkardt
  PUBLIC :: quadratic                         ! Simple quadratic, (x-2)^2+1.
  PUBLIC :: quadratic_exponential             ! Quadratic plus exponential, x^2 + e^(-x).
  PUBLIC :: quartic                           ! Quartic, x^4 + 2x^2 + x + 3.
  PUBLIC :: steep_valley                      ! Steep valley, e^x + 1/(100x).
  PUBLIC :: steep_valley2                     ! Steep valley, e^x - 2x + 1/(100x) - 1/(1000000x^2)
  PUBLIC :: dying_snake                       ! The dying snake, ( x + sin(x) ) * e^(-x^2).
  PUBLIC :: thin_pole                         ! The "Thin Pole", x^2+1+log((pi-x)^2)/pi^4
  PUBLIC :: oscillatory_parabola              ! The oscillatory parabola
  PUBLIC :: cosine_combo                      ! The cosine combo
  PUBLIC :: abs1                              ! 1 + |3x-1|
  ! ------------------------------------------------------------------
  ! test_opt package of John Burkardt
  PUBLIC :: fletcher_powell_helical_valley    ! The Fletcher-Powell helical valley function, N = 3.
  PUBLIC :: biggs_exp6                        ! The Biggs EXP6 function, N = 6.
  PUBLIC :: gaussian                          ! The Gaussian function, N = 3.
  PUBLIC :: powell_badly_scaled               ! The Powell badly scaled function, N = 2.
  PUBLIC :: box_3dimensional                  ! The Box 3-dimensional function, N = 3.
  PUBLIC :: variably_dimensioned              ! The variably dimensioned function, 1 <= N.
  PUBLIC :: watson                            ! The Watson function, 2 <= N.
  PUBLIC :: penalty1                          ! The penalty function #1, 1 <= N.
  PUBLIC :: penalty2                          ! The penalty function #2, 1 <= N.
  PUBLIC :: brown_badly_scaled                ! The Brown badly scaled function, N = 2.
  PUBLIC :: brown_dennis                      ! The Brown and Dennis function, N = 4.
  PUBLIC :: gulf_rd                           ! The Gulf R&D function, N = 3.
  PUBLIC :: trigonometric                     ! The trigonometric function, 1 <= N.
  PUBLIC :: ext_rosenbrock_parabolic_valley   ! The extended Rosenbrock parabolic valley function, 1 <= N.
  PUBLIC :: ext_powell_singular_quartic       ! The extended Powell singular quartic function, 4 <= N.
  PUBLIC :: beale                             ! The Beale function, N = 2.
  PUBLIC :: wood                              ! The Wood function, N = 4.
  PUBLIC :: chebyquad                         ! The Chebyquad function, 1 <= N.
  PUBLIC :: leon_cubic_valley                 ! Leon''s cubic valley function, N = 2.
  PUBLIC :: gregory_karney_tridia_matrix      ! Gregory and Karney''s Tridiagonal Matrix Function, 1 <= N.
  PUBLIC :: hilbert                           ! The Hilbert function, 1 <= N.
  PUBLIC :: de_jong_f1                        ! The De Jong Function F1, N = 3.
  PUBLIC :: de_jong_f2                        ! The De Jong Function F2, N = 2.
  PUBLIC :: de_jong_f3                        ! The De Jong Function F3 (discontinuous), N = 5.
  PUBLIC :: de_jong_f4                        ! The De Jong Function F4 (Gaussian noise), N = 30.
  PUBLIC :: de_jong_f5                        ! The De Jong Function F5, N = 2.
  PUBLIC :: schaffer_f6                       ! The Schaffer Function F6, N = 2.
  PUBLIC :: schaffer_f7                       ! The Schaffer Function F7, N = 2.
  PUBLIC :: goldstein_price_polynomial        ! The Goldstein Price Polynomial, N = 2.
  PUBLIC :: branin_rcos                       ! The Branin RCOS Function, N = 2.
  PUBLIC :: shekel_sqrn5                      ! The Shekel SQRN5 Function, N = 4.
  PUBLIC :: shekel_sqrn7                      ! The Shekel SQRN7 Function, N = 4.
  PUBLIC :: shekel_sqrn10                     ! The Shekel SQRN10 Function, N = 4.
  PUBLIC :: six_hump_camel_back_polynomial    ! The Six-Hump Camel-Back Polynomial, N = 2.
  PUBLIC :: shubert                           ! The Shubert Function, N = 2.
  PUBLIC :: stuckman                          ! The Stuckman Function, N = 2.
  PUBLIC :: easom                             ! The Easom Function, N = 2.
  PUBLIC :: bohachevsky1                      ! The Bohachevsky Function #1, N = 2.
  PUBLIC :: bohachevsky2                      ! The Bohachevsky Function #2, N = 2.
  PUBLIC :: bohachevsky3                      ! The Bohachevsky Function #3, N = 2.
  PUBLIC :: colville_polynomial               ! The Colville Polynomial, N = 4.
  PUBLIC :: powell3d                          ! The Powell 3D function, N = 3.
  PUBLIC :: himmelblau                        ! The Himmelblau function, N = 2.
  ! ------------------------------------------------------------------
  ! Miscellaneous functions
  PUBLIC :: griewank                          ! Griewank function, N = 2 or N = 10.
  PUBLIC :: rosenbrock                        ! Rosenbrock parabolic valley function, N = 2.
  PUBLIC :: sphere_model                      ! Sphere model, N >= 1.
  PUBLIC :: rastrigin                         ! Rastrigin function, N >= 2.
  PUBLIC :: schwefel                          ! Schwefel function, N >= 2.
  PUBLIC :: ackley                            ! Ackley function, N >= 2.
  PUBLIC :: michalewicz                       ! Michalewicz function, N >= 2.
  PUBLIC :: booth                             ! Booth function, N = 2.
  PUBLIC :: hump                              ! Hump function, N = 2.
  PUBLIC :: levy                              ! Levy function, N >= 2.
  PUBLIC :: matyas                            ! Matyas function, N = 2.
  PUBLIC :: perm                              ! Perm function, N = 4.
  PUBLIC :: power_sum                         ! Power sum function, N = 4.
  ! ------------------------------------------------------------------
  ! test_optimization package of John Burkardt - inputs are x(m,n) and output f(n), e.g. compare
  ! rosenbrock = 100.0_dp * (x(2)-x(1)**2)**2 + (1.0_dp-x(1))**2
  ! rosenbrock_2d(j) = sum((1.0_dp-x(1:m,j))**2) + sum((x(2:m,j)-x(1:m-1,j))**2)
  PUBLIC :: sphere_model_2d                   !  The sphere model, (M,N).
  PUBLIC :: axis_parallel_hyper_ellips_2d     !  The axis-parallel hyper-ellipsoid function, (M,N).
  PUBLIC :: rotated_hyper_ellipsoid_2d        !  The rotated hyper-ellipsoid function, (M,N).
  PUBLIC :: rosenbrock_2d                     !  Rosenbrock''s valley, (M,N).
  PUBLIC :: rastrigin_2d                      !  Rastrigin''s function, (M,N).
  PUBLIC :: schwefel_2d                       !  Schwefel''s function, (M,N).
  PUBLIC :: griewank_2d                       !  Griewank''s function, (M,N).
  PUBLIC :: power_sum_2d                      !  The power sum function, (M,N).
  PUBLIC :: ackley_2d                         !  Ackley''s function, (M,N).
  PUBLIC :: michalewicz_2d                    !  Michalewicz''s function, (M,N).
  PUBLIC :: drop_wave_2d                      !  The drop wave function, (M,N).
  PUBLIC :: deceptive_2d                      !  The deceptive function, (M,N).

CONTAINS

  ! ------------------------------------------------------------------
  !
  !  Simple quadratic, (x-2)^2+1
  !  Solution: x = 2.0
  !  With Brent method:
  !   A,  X*,  B:   1.9999996       2.0000000       2.0000004
  !  FA, FX*, FB:   1.0000000       1.0000000       1.0000000
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    25 February 2002
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function quadratic( x )

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: quadratic

    if (size(x,1) .gt. 1_i4) stop 'quadratic: Input has to be array of size 1'
    quadratic = ( x(1) - 2.0_dp ) * ( x(1) - 2.0_dp ) + 1.0_dp

  end function quadratic

  ! ------------------------------------------------------------------
  !
  !  Quadratic plus exponential, x^2 + e^(-x)
  !  Solution: x = 0.35173370
  !  With Brent method:
  !   A,  X*,  B:  0.35173337      0.35173370      0.35173404
  !  FA, FX*, FB:  0.82718403      0.82718403      0.82718403
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 February 2002
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    LE Scales,
  !    Introduction to Non-Linear Optimization,
  !    Springer, 1985.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function quadratic_exponential( x )

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: quadratic_exponential

    if (size(x,1) .gt. 1_i4) stop 'quadratic_exponential: Input has to be array of size 1'
    quadratic_exponential = x(1) * x(1) + exp ( - x(1) )

  end function quadratic_exponential

  ! ------------------------------------------------------------------
  !
  !  Quartic, x^4 + 2x^2 + x + 3
  !  Solution: x = -0.23673291
  !  With Brent method:
  !   A,  X*,  B: -0.23673324     -0.23673291     -0.23673257
  !  FA, FX*, FB:   2.8784928       2.8784928       2.8784928
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 February 2002
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    LE Scales,
  !    Introduction to Non-Linear Optimization,
  !    Springer, 1985.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function quartic( x )

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: quartic

    if (size(x,1) .gt. 1_i4) stop 'quartic: Input has to be array of size 1'
    quartic = ( ( x(1) * x(1) + 2.0_dp ) * x(1) + 1.0_dp ) * x(1) + 3.0_dp

  end function quartic

  ! ------------------------------------------------------------------
  !
  !  Steep valley, e^x + 1/(100x)
  !  Solution:
  !       if x >  0.0 :   x = 0.95344636E-01
  !       if x < -0.1 :   x = -8.99951
  !  Search domain: x <= -0.1
  !
  !  With Brent method:
  !   A,  X*,  B:  0.95344301E-01  0.95344636E-01  0.95344971E-01
  !  FA, FX*, FB:   1.2049206       1.2049206       1.2049206
  !
  !  Discussion:
  !
  !    This function has a pole at x = 0,
  !    near which
  !       f -> negative infinity for x -> 0-0 (left) and
  !       f -> positive infinity for x -> 0+0 (right)
  !    f has a local maximum at x ~ -0.105412 .
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 February 2002
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    LE Scales,
  !    Introduction to Non-Linear Optimization,
  !    Springer, 1985.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function steep_valley( x )

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: steep_valley

    if (size(x,1) .gt. 1_i4) stop 'steep_valley: Input has to be array of size 1'
    steep_valley = exp ( x(1) ) + 0.01_dp / x(1)

    steep_valley = steep_valley + 0.0009877013_dp

  end function steep_valley

  ! ------------------------------------------------------------------
  !
  !  Steep valley2, e^x - 2x + 1/(100x) - 1/(1000000x^2)
  !
  !  Solution: x = 0.70320487
  !  Search domain: 0.001 <= x
  !
  !  With Brent method:
  !   A,  X*,  B:  0.70320453      0.70320487      0.70320521
  !  FA, FX*, FB:  0.62802572      0.62802572      0.62802572
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 February 2002
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    LE Scales,
  !    Introduction to Non-Linear Optimization,
  !    Springer, 1985.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function steep_valley2( x )

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: steep_valley2

    if (size(x,1) .gt. 1_i4) stop 'steep_valley2: Input has to be array of size 1'
    steep_valley2 = exp( x(1) ) - 2.0_dp * x(1) + 0.01_dp / x(1) - 0.000001_dp / x(1) / x(1)

  end function steep_valley2

  ! ------------------------------------------------------------------
  !
  !  The dying snake, ( x + sin(x) ) * e^(-x^2)
  !  Solution: x = -0.67957876
  !  With Brent method:
  !   A,  X*,  B: -0.67957911     -0.67957876     -0.67957842
  !  FA, FX*, FB: -0.82423940     -0.82423940     -0.82423940
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 February 2002
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization Without Derivatives,
  !    Prentice Hall 1973,
  !    Reprinted Dover, 2002
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function dying_snake(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: dying_snake

    if (size(x,1) .gt. 1_i4) stop 'dying_snake: Input has to be array of size 1'
    dying_snake = ( x(1) + sin ( x(1) ) ) * exp ( - x(1) * x(1) )

    dying_snake = dying_snake + 0.8242393985_dp

  end function dying_snake

  ! ------------------------------------------------------------------
  !
  !  The "Thin Pole", 3x^2+1+log((pi-x)^2)/pi^4
  !  Solution:
  !     x    = 0.00108963
  !     f(x) = 1.0235
  !
  !  With Brent method:
  !   A,  X*,  B:   2.0000000       2.0000007       2.0000011
  !  FA, FX*, FB:   13.002719       13.002727       13.002732
  !
  !  Discussion:
  !
  !    This function looks positive, but has a pole at x = pi,
  !    near which f -> negative infinity, and has two zeroes nearby.
  !    None of this will show up computationally.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 February 2003
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Arnold Krommer, Christoph Ueberhuber,
  !    Numerical Integration on Advanced Systems,
  !    Springer, 1994, pages 185-186.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function thin_pole(x)

    use mo_constants, only: pi_dp
    use mo_utils,     only: eq

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: thin_pole

    if (size(x,1) .gt. 1_i4) stop 'thin_pole: Input has to be array of size 1'
    if ( eq(x(1),pi_dp) ) then
       thin_pole = - 10000.0_dp
    else
       thin_pole = 3.0_dp * x(1) * x(1) + 1.0_dp + ( log ( ( x(1) - pi_dp ) * ( x(1) - pi_dp ) ) ) / pi_dp**4
    end if

  end function thin_pole

  ! ------------------------------------------------------------------
  !
  !  The oscillatory parabola x^2 - 10*sin(x^2-3x+2)
  !  Solution:
  !     x    = 0.146623
  !     f(x) = -9.97791
  !  With Brent method:
  !   A,  X*,  B:  -1.3384524      -1.3384521      -1.3384517
  !  FA, FX*, FB:  -8.1974224      -8.1974224      -8.1974224
  !
  !  Discussion:
  !
  !    This function is oscillatory, with many local minima.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    25 January 2008
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function oscillatory_parabola(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: oscillatory_parabola

    if (size(x,1) .gt. 1_i4) stop 'oscillatory_parabola: Input has to be array of size 1'
    oscillatory_parabola = x(1) * x(1) - 10.0_dp * sin ( x(1) * x(1) - 3.0_dp * x(1) + 2.0_dp )

    oscillatory_parabola = oscillatory_parabola + 9.9779149346_dp

  end function oscillatory_parabola

  ! ------------------------------------------------------------------
  !
  !  The cosine combo cos(x)+5cos(1.6x)-2cos(2x)+5cos(4.5x)+7cos(9x)
  !  Solution:
  !     x    = -21.9443 + 62.831853 * k   , k = Integer
  !     x    =  21.9443 - 62.831853 * k   , k = Integer
  !     f(x) = -14.6772
  !
  !  With Brent method:
  !   A,  X*,  B:   1.0167817       1.0167821       1.0167824
  !  FA, FX*, FB:  -6.2827509      -6.2827509      -6.2827509
  !
  !  Discussion:
  !
  !    This function is symmetric, oscillatory, and has many local minima.
  !
  !    The function has a local minimum at 1.0167817.
  !
  !    The global optimum which function value -14.6772
  !    appears infinite times.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 February 2009
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Isabel Beichl, Dianne O'Leary, Francis Sullivan,
  !    Monte Carlo Minimization and Counting: One, Two, Too Many,
  !    Computing in Science and Engineering,
  !    Volume 9, Number 1, January/February 2007.
  !
  !    Dianne O'Leary,
  !    Scientific Computing with Case Studies,
  !    SIAM, 2008,
  !    ISBN13: 978-0-898716-66-5,
  !    LC: QA401.O44.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function cosine_combo(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: cosine_combo

    if (size(x,1) .gt. 1_i4) stop 'cosine_combo: Input has to be array of size 1'
    cosine_combo =   cos (         x(1) ) &
         + 5.0_dp * cos ( 1.6_dp * x(1) ) &
         - 2.0_dp * cos ( 2.0_dp * x(1) ) &
         + 5.0_dp * cos ( 4.5_dp * x(1) ) &
         + 7.0_dp * cos ( 9.0_dp * x(1) )

    cosine_combo = cosine_combo + 14.6771885214_dp

  end function cosine_combo

  ! ------------------------------------------------------------------
  !
  !  abs1, 1 + |3x-1|
  !  Solution: x = 1./3.
  !  With Brent method:
  !   A,  X*,  B:  0.33333299      0.33333351      0.33333385
  !  FA, FX*, FB:   1.0000010       1.0000005       1.0000015
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 February 2012
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument of the objective function.
  !

  function abs1(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: abs1

    if (size(x,1) .gt. 1_i4) stop 'abs1: Input has to be array of size 1'
    abs1 = 1.0_dp + abs ( 3.0_dp * x(1) - 1.0_dp )

  end function abs1

  ! ------------------------------------------------------------------
  !
  !  R8_AINT truncates an R8 argument to an integer.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 October 2011
  !
  !  Author:
  !
  !    John Burkardt.
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, the argument.
  !
  !    Output, real(dp) VALUE, the truncated version of X.
  !

  function r8_aint( x )

    implicit none

    real(dp) :: r8_aint
    real(dp) :: val
    real(dp) :: x

    if ( x < 0.0_dp ) then
       val = -int( abs ( x ) )
    else
       val =  int( abs ( x ) )
    end if

    r8_aint = val

  end function r8_aint

  !*****************************************************************************80
  !
  !! NORMAL_01_SAMPLE samples the standard Normal PDF.
  !
  !  Discussion:
  !
  !    The standard normal distribution has mean 0 and standard
  !    deviation 1.
  !
  !  Method:
  !
  !    The Box-Muller method is used.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Output, real(dp) X, a sample of the PDF.
  !

  subroutine normal_01_sample ( x )

    use mo_constants, only: pi_dp
    use mo_utils, only: le

    implicit none

    integer(i4), save :: iset = -1
    real(dp) v1
    real(dp) v2
    real(dp) x
    real(dp), save :: xsave = 0.0_dp

    if ( iset == -1 ) then
       call random_seed ( )
       iset = 0
    end if

    if ( iset == 0 ) then

       call random_number ( harvest = v1 )

       if ( le(v1,0.0_dp) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'NORMAL_01_SAMPLE - Fatal error!'
          write ( *, '(a)' ) '  V1 <= 0.'
          write ( *, '(a,g14.6)' ) '  V1 = ', v1
          stop
       end if

       call random_number ( harvest = v2 )

       if ( le(v2,0.0_dp) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'NORMAL_01_SAMPLE - Fatal error!'
          write ( *, '(a)' ) '  V2 <= 0.'
          write ( *, '(a,g14.6)' ) '  V2 = ', v2
          stop
       end if

       x = sqrt ( - 2.0_dp * log ( v1 ) ) * cos ( 2.0_dp * pi_dp * v2 )

       xsave = sqrt ( - 2.0_dp * log ( v1 ) ) * sin ( 2.0_dp * PI_dp * v2 )

       iset = 1

    else

       x = xsave
       iset = 0

    end if

    return
  end subroutine normal_01_sample

  ! ------------------------------------------------------------------
  !
  ! The Fletcher-Powell helical valley function, N = 3.
  ! Solution: x(1:3) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function fletcher_powell_helical_valley(x)

    use mo_constants, only: pi_dp

    implicit none

    !    integer(i4) :: n

    real(dp) :: fletcher_powell_helical_valley
    real(dp) :: th
    real(dp), dimension(:), intent(in) :: x

    if ( 0.0_dp < x(1) ) then
       th = 0.5_dp * atan ( x(2) / x(1) ) / pi_dp
    else if ( x(1) < 0.0_dp ) then
       th = 0.5_dp * atan ( x(2) / x(1) ) / pi_dp + 0.5_dp
    else if ( 0.0_dp < x(2) ) then
       th = 0.25_dp
    else if ( x(2) < 0.0_dp ) then
       th = - 0.25_dp
    else
       th = 0.0_dp
    end if
    !call p01_th ( x, th )

    fletcher_powell_helical_valley = 100.0_dp * ( x(3) - 10.0_dp * th )**2 &
         + 100.0_dp * ( sqrt ( x(1) * x(1) + x(2) * x(2) ) - 1.0_dp )**2 &
         + x(3) * x(3)

  end function fletcher_powell_helical_valley

  ! ------------------------------------------------------------------
  !
  ! The Biggs EXP6 function, N = 6.
  ! Solution: x(1:6) = (/ 1.0_dp, 10.0_dp, 1.0_dp, 5.0_dp, 4.0_dp, 3.0_dp /)
  !           at f(x*) = 0.0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function biggs_exp6(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: c
    real(dp) :: biggs_exp6
    real(dp) :: fi
    integer(i4) ::i
    real(dp), dimension(:), intent(in) :: x

    biggs_exp6 = 0.0_dp

    do i = 1, 13

       c = - real ( i, dp ) / 10.0_dp

       fi = x(3)     * exp ( c * x(1) )         - x(4) * exp ( c * x(2) ) &
            + x(6)     * exp ( c * x(5) )         -        exp ( c ) &
            + 5.0_dp  * exp ( 10.0_dp * c ) - 3.0_dp  * exp ( 4.0_dp * c )

       biggs_exp6 = biggs_exp6 + fi * fi

    end do

  end function biggs_exp6

  ! ------------------------------------------------------------------
  !
  ! The Gaussian function, N = 3.
  ! Search domain: -Pi <= xi <= Pi, i = 1, 2, 3.
  ! Solution:
  !     x(1:n) = (/ 0.39895613783875655_dp, 1.0000190844878036_dp, 0.0_dp /)
  !     at f(x*) = 0.0
  !     found with Mathematica
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function gaussian(x)

    implicit none

    !    integer(i4) :: n

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: gaussian

    integer(i4) :: i
    real(dp)    :: t
    real(dp)    :: y(15)

    y(1:15) = (/ 0.0009_dp, 0.0044_dp, 0.0175_dp, 0.0540_dp, 0.1295_dp, &
         0.2420_dp, 0.3521_dp, 0.3989_dp, 0.3521_dp, 0.2420_dp, &
         0.1295_dp, 0.0540_dp, 0.0175_dp, 0.0044_dp, 0.0009_dp /)

    gaussian = 0.0_dp

    do i = 1, 15

       ! avoiding underflow
       t =  - 0.5_dp * x(2) * &
            ( 3.5_dp - 0.5_dp * real ( i - 1, dp ) - x(3) )**2
       if ( t .lt. -709._dp ) then
          t = -y(i)
       else
          t = x(1) * exp ( t ) - y(i)
       end if

       gaussian = gaussian + t * t

    end do

  end function gaussian

  ! ------------------------------------------------------------------
  !
  ! The Powell badly scaled function, N = 2.
  ! Solution: x(1:2) = (/ 1.098159D-05, 9.106146_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function powell_badly_scaled(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: powell_badly_scaled
    real(dp) :: f1
    real(dp) :: f2
    real(dp), dimension(:), intent(in) :: x

    f1 = 10000.0_dp * x(1) * x(2) - 1.0_dp
    f2 = exp ( - x(1) ) + exp ( - x(2) ) - 1.0001_dp

    powell_badly_scaled = f1 * f1 + f2 * f2

  end function powell_badly_scaled

  ! ------------------------------------------------------------------
  !
  ! The Box 3-dimensional function, N = 3.
  ! Solution: x(1:3) = (/ 1.0_dp, 10.0_dp, 1.0_dp /)
  ! seems to be not the only solution
  !
  !  Discussion:
  !
  !    The function is formed by the sum of squares of 10 separate terms.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function box_3dimensional(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: c
    real(dp) :: box_3dimensional
    real(dp) :: fi
    integer(i4) ::i
    real(dp), dimension(:), intent(in) :: x

    box_3dimensional = 0.0_dp

    do i = 1, 10

       c = - real ( i, dp ) / 10.0_dp

       fi = exp ( c * x(1) ) - exp ( c * x(2) ) - x(3) * &
            ( exp ( c ) - exp ( 10.0_dp * c ) )

       box_3dimensional = box_3dimensional + fi * fi

    end do

  end function box_3dimensional

  ! ------------------------------------------------------------------
  !
  ! The variably dimensioned function, 1 <= N.
  ! Solution: x(1:n) = 1.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function variably_dimensioned(x)

    implicit none

    integer(i4) :: n

    real(dp) :: variably_dimensioned
    real(dp) :: f1
    real(dp) :: f2
    integer(i4) ::i
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    f1 = 0.0_dp
    do i = 1, n
       f1 = f1 + real ( i, dp ) * ( x(i) - 1.0_dp )
    end do

    f2 = 0.0_dp
    do i = 1, n
       f2 = f2 + ( x(i) - 1.0_dp )**2
    end do

    variably_dimensioned = f1 * f1 * ( 1.0_dp + f1 * f1 ) + f2

  end function variably_dimensioned

  ! ------------------------------------------------------------------
  !
  ! The Watson function, 2 <= N.
  ! Solution: n==6: x(1:n) = (/ -0.015725_dp, 1.012435_dp, -0.232992_dp,
  !                             1.260430_dp, -1.513729_dp, 0.992996_dp /)
  !           n==9  x(1:n) = (/ -0.000015_dp, 0.999790_dp, 0.014764_dp,
  !                              0.146342_dp, 1.000821_dp, -2.617731_dp,
  !                              4.104403_dp, -3.143612_dp, 1.052627_dp /)
  !           else unknown
  !
  !  Discussion:
  !
  !    For N = 9, the problem of minimizing the Watson function is
  !    very ill conditioned.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function watson(x)

    implicit none

    integer(i4) :: n

    real(dp) :: d
    real(dp) :: watson
    integer(i4) ::i
    integer(i4) ::j
    real(dp) :: s1
    real(dp) :: s2
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    watson = 0.0_dp
    do i = 1, 29

       s1 = 0.0_dp
       d = 1.0_dp
       do j = 2, n
          s1 = s1 + real ( j - 1, dp ) * d * x(j)
          d = d * real ( i, dp ) / 29.0_dp
       end do

       s2 = 0.0_dp
       d = 1.0_dp
       do j = 1, n
          s2 = s2 + d * x(j)
          d = d * real ( i, dp ) / 29.0_dp
       end do

       watson = watson + ( s1 - s2 * s2 - 1.0_dp )**2

    end do

    watson = watson + x(1) * x(1) + ( x(2) - x(1) * x(1) - 1.0_dp )**2

  end function watson

  ! ------------------------------------------------------------------
  !
  ! The penalty function #1, 1 <= N.
  ! Solution:
  !    with Mathematica (numerical results)
  !    x(1:1) = 0.5000049998750047_dp
  !    f(x)   = 0.000002499975000499985_dp
  !
  !    x(1:2) = 0.35355985481744073_dp
  !    f(x)   = 0.000008357780799989139_dp
  !
  !    x(1:3) = 0.2886234818387535_dp
  !    f(x)   = 0.000015179340383244187_dp
  !
  !    x(1:4) = 0.2500074995875379_dp
  !    f(x)   = 0.000022499775008999372_dp
  !
  !    x(1:5) = 0.2236145612000511_dp
  !    f(x)   = 0.000030139018845277502_dp
  !
  !    x(1:6) = 0.20413210344548943_dp
  !    f(x)   = 0.00003800472253975947_dp
  !
  !    x(1:7) = 0.18899034607915027_dp
  !    f(x)   = 0.00004604202648884626_dp
  !
  !    x(1:8) = 0.1767849268724064_dp
  !    f(x)   = 0.00005421518662591646_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function penalty1(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: penalty1

    integer(i4)         :: n
    real(dp), parameter :: ap = 0.00001_dp
    real(dp)            :: t1
    real(dp)            :: t2

    n = size(x)
    t1 = - 0.25_dp + sum ( x(1:n)**2 )

    t2 = sum ( ( x(1:n) - 1.0_dp )**2 )

    penalty1 = ap * t2 + t1 * t1

  end function penalty1

  ! ------------------------------------------------------------------
  !
  ! The penalty function #2, 1 <= N.
  ! Solution:
  !    with Mathematica (numerical results)
  !
  !    x(1:1) = (/ 0.7914254791661984_dp /)
  !    f(x)   = 0.4893952147007766_dp
  !
  !    x(1:2) = (/ 0.2000002053277478_dp, 0.95916622198851_dp /)
  !    f(x)   = 0.000000806639004111886_dp
  !
  !    x(1:3) = (/ 0.19999990827773015_dp, 0.521451491987707_dp, 0.5798077888356243_dp /)
  !    f(x)   = 0.0000031981885430064677_dp
  !
  !    x(1:4) = (/ 0.19999930547325262_dp, 0.19165582942317613_dp, 0.47960388408453586_dp, &
  !                0.519390092116238_dp /)
  !    f(x)   = 0.000009376294404291533_dp
  !
  !    x(1:5) = (/ 0.19999832542354226_dp, 0.09277573718478778_dp, 0.20939661885734498_dp, &
  !                0.4467231497703405_dp, 0.4846761802898935_dp /)
  !    f(x)   = 0.000021387590981336432_dp
  !
  !    x(1:6) = (/ 0.19999687534275826_dp, 0.05202606784766213_dp, 0.11181980463664613_dp, &
  !                0.21122219156860741_dp, 0.4237150785565095_dp, 0.45116217266910164_dp /)
  !    f(x)   = 0.00004193126194907272_dp
  !
  !    x(1:7) = (/ 0.1999947753716143_dp,  0.03223698574640245_dp, 0.06616491022124847_dp, &
  !                0.11800087933141483_dp, 0.2086993295142654_dp, 0.4018822793193398_dp, &
  !                0.4272124489678876_dp /)
  !    f(x)   = 0.00007445577533975632_dp
  !
  !    x(1:8) = (/ 0.19999196971962244_dp, 0.02124196902699145_dp, 0.04164127214921979_dp, &
  !                0.0721284757968418_dp, 0.11998769382443852_dp, 0.20359198648845672_dp, &
  !                0.3808503545095562_dp,  0.41039238853670246_dp /)
  !    f(x)   = 0.0001233351431976925_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function penalty2(x)

    implicit none

    integer(i4) :: n

    real(dp), parameter :: ap = 0.00001_dp
    real(dp) :: d2
    real(dp) :: penalty2
    integer(i4) ::j
    real(dp) :: s1
    real(dp) :: s2
    real(dp) :: s3
    real(dp) :: t1
    real(dp) :: t2
    real(dp) :: t3
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    t1 = -1.0_dp
    t2 = 0.0_dp
    t3 = 0.0_dp
    d2 = 1.0_dp
    s2 = 0.0_dp
    do j = 1, n
       t1 = t1 + real ( n - j + 1, dp ) * x(j)**2
       s1 = exp ( x(j) / 10.0_dp )
       if ( 1 < j ) then
          s3 = s1 + s2 - d2 * ( exp ( 0.1_dp ) + 1.0_dp )
          t2 = t2 + s3 * s3
          t3 = t3 + ( s1 - 1.0_dp / exp ( 0.1_dp ) )**2
       end if
       s2 = s1
       d2 = d2 * exp ( 0.1_dp )
    end do

    penalty2 = ap * ( t2 + t3 ) + t1 * t1 + ( x(1) - 0.2_dp )**2

  end function penalty2

  ! ------------------------------------------------------------------
  !
  ! The Brown badly scaled function, N = 2.
  ! Solution: x(1:2) = (/ 1.0D+06, 2.0D-06 /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function brown_badly_scaled(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: brown_badly_scaled
    real(dp), dimension(:), intent(in) :: x

    brown_badly_scaled = ( x(1) - 1000000.0_dp )**2 &
         + ( x(2) - 0.000002_dp )**2 &
         + ( x(1) * x(2) - 2.0_dp )**2

  end function brown_badly_scaled

  ! ------------------------------------------------------------------
  !
  ! The Brown and Dennis function, N = 4.
  ! Solution:
  !    x(1:4) = (/ -11.5944399230_dp, 13.2036300657_dp, &
  !                -0.4034395329_dp,   0.2367787297_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function brown_dennis(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: c
    real(dp) :: brown_dennis
    real(dp) :: f1
    real(dp) :: f2
    integer(i4) ::i
    real(dp), dimension(:), intent(in) :: x

    brown_dennis = 0.0_dp

    do i = 1, 20

       c = real ( i, dp ) / 5.0_dp
       f1 = x(1) + c * x(2) - exp ( c )
       f2 = x(3) + sin ( c ) * x(4) - cos ( c )

       brown_dennis = brown_dennis + f1**4 + 2.0_dp * f1 * f1 * f2 * f2 + f2**4

    end do

    brown_dennis = brown_dennis - 85822.2016263563_dp

  end function brown_dennis

  ! ------------------------------------------------------------------
  !
  ! The Gulf R&D function, N = 3.
  ! Solution: x(1:3) = (/ 50.0_dp, 25.0_dp, 1.5_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function gulf_rd(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: arg
    real(dp) :: gulf_rd
    integer(i4) ::i
    real(dp) :: r
    real(dp) :: t
    real(dp), dimension(:), intent(in) :: x
    real(dp) :: sqrtHuge

    sqrtHuge = Huge(1.0_dp)**0.5_dp -1.0_dp

    gulf_rd = 0.0_dp
    do i = 1, 99

       arg = real ( i, dp ) / 100.0_dp
       r = abs ( ( - 50.0_dp * log ( arg ) )**( 2.0_dp / 3.0_dp ) &
            + 25.0_dp - x(2) )

       ! print*, 'arg         = ',arg
       ! print*, 'r           = ',r
       ! print*, 'x(1)        = ',x(1)
       ! print*, 'x(2)        = ',x(2)
       ! print*, 'x(3)        = ',x(3)

       ! avoiding underflow
       if ( x(3)*Log(r) .gt. -708.-dp ) then
          print*, '-exp(x(3)*Log(r)) = ',-exp(x(3)*Log(r))
          print*, 'x(1)              = ',x(1)
          if ( -exp(x(3)*Log(r)) / x(1) .lt. -708._dp) then
             t = -arg
          else
             if ( -r**x(3) / x(1) .gt. 708._dp) then
                t = 1000000._dp - arg
             else
                t = exp ( - r**x(3) / x(1) ) - arg
             end if
          end if
       else
          t = -arg
       end if

       if ( abs(t) .gt. sqrtHuge ) then
          ! under/overflow case
          t = sqrtHuge
       end if

       if ( Huge(1.0_dp) -gulf_rd .gt. t*t ) then
          ! usual case
          gulf_rd = gulf_rd + t * t
       else
          ! overflow case
          gulf_rd = Huge(1.0_dp)
       end if

    end do

  end function gulf_rd

  ! ------------------------------------------------------------------
  !
  ! The trigonometric function, 1 <= N.
  ! Solution: x(1:n) = 0.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function trigonometric(x)

    implicit none

    integer(i4) :: n

    real(dp) :: trigonometric
    integer(i4) ::j
    real(dp) :: s1
    real(dp) :: t
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    s1 = sum ( cos ( x(1:n) ) )

    trigonometric = 0.0_dp
    do j = 1, n
       t = real ( n + j, dp ) - sin ( x(j) ) &
            - s1 - real ( j, dp ) * cos ( x(j) )
       trigonometric = trigonometric + t * t
    end do

  end function trigonometric

  ! ------------------------------------------------------------------
  !
  ! The extended Rosenbrock parabolic valley function, 1 <= N.
  ! Solution: x(1:n) = 1.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function ext_rosenbrock_parabolic_valley(x)

    implicit none

    integer(i4) :: n

    real(dp) :: ext_rosenbrock_parabolic_valley
    integer(i4) ::j
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    ext_rosenbrock_parabolic_valley = 0.0_dp
    do j = 1, n
       if ( mod ( j, 2 ) == 1 ) then
          ext_rosenbrock_parabolic_valley = ext_rosenbrock_parabolic_valley + ( 1.0_dp - x(j) )**2
       else
          ext_rosenbrock_parabolic_valley = ext_rosenbrock_parabolic_valley + 100.0_dp * ( x(j) - x(j-1) * x(j-1) )**2
       end if
    end do

  end function ext_rosenbrock_parabolic_valley

  ! ------------------------------------------------------------------
  !
  ! The extended Powell singular quartic function, 4 <= N.
  ! Solution: x(1:n) = 0.0_dp
  !
  !  Discussion:
  !
  !    The Hessian matrix is doubly singular at the minimizer,
  !    suggesting that most optimization routines will experience
  !    a severe slowdown in convergence.
  !
  !    The problem is usually only defined for N being a multiple of 4.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function ext_powell_singular_quartic(x)

    implicit none

    integer(i4) :: n

    real(dp) :: ext_powell_singular_quartic
    real(dp) :: f1
    real(dp) :: f2
    real(dp) :: f3
    real(dp) :: f4
    integer(i4) ::j
    real(dp), dimension(:), intent(in) :: x
    real(dp) :: xjp1
    real(dp) :: xjp2
    real(dp) :: xjp3

    n = size(x)
    ext_powell_singular_quartic = 0.0_dp

    do j = 1, n, 4

       if ( j + 1 <= n ) then
          xjp1 = x(j+1)
       else
          xjp1 = 0.0_dp
       end if

       if ( j + 2 <= n ) then
          xjp2 = x(j+2)
       else
          xjp2 = 0.0_dp
       end if

       if ( j + 3 <= n ) then
          xjp3 = x(j+3)
       else
          xjp3 = 0.0_dp
       end if

       f1 = x(j) + 10.0_dp * xjp1

       if ( j + 1 <= n ) then
          f2 = xjp2 - xjp3
       else
          f2 = 0.0_dp
       end if

       if ( j + 2 <= n ) then
          f3 = xjp1 - 2.0_dp * xjp2
       else
          f3 = 0.0_dp
       end if

       if ( j + 3 <= n ) then
          f4 = x(j) - xjp3
       else
          f4 = 0.0_dp
       end if

       ext_powell_singular_quartic = ext_powell_singular_quartic +            f1 * f1 &
            +  5.0_dp * f2 * f2 &
            +            f3 * f3 * f3 * f3 &
            + 10.0_dp * f4 * f4 * f4 * f4

    end do

  end function ext_powell_singular_quartic

  ! ------------------------------------------------------------------
  !
  ! The Beale function, N = 2.
  ! Solution: x(1:2) = (/ 3.0_dp, 0.5_dp /)
  ! Search domain: -4.5 <= xi <= 4.5, i = 1, 2.
  !
  !  Discussion:
  !
  !    Range according to:
  !       http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/
  !       Hedar_files/TestGO_files/Page288.htm
  !
  !    This function has a valley approaching the line X(2) = 1.
  !
  !    The function has a global minimizer:
  !
  !      X(*) = ( 3.0, 0.5 ), F(X*) = 0.0.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 January 2008
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Evelyn Beale,
  !    On an Iterative Method for Finding a Local Minimum of a Function
  !    of More than One Variable,
  !    Technical Report 25,
  !    Statistical Techniques Research Group,
  !    Princeton University, 1958.
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function beale(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: beale
    real(dp) :: f1
    real(dp) :: f2
    real(dp) :: f3
    real(dp), dimension(:), intent(in) :: x

    f1 = 1.5_dp   - x(1) * ( 1.0_dp - x(2)    )
    f2 = 2.25_dp  - x(1) * ( 1.0_dp - x(2) * x(2) )
    f3 = 2.625_dp - x(1) * ( 1.0_dp - x(2) * x(2) * x(2) )

    beale = f1 * f1 + f2 * f2 + f3 * f3

  end function beale

  ! ------------------------------------------------------------------
  !
  ! The Wood function, N = 4.
  ! Solution: x(1:4) = (/ 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    06 January 2008
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function wood(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: wood
    real(dp) :: f1
    real(dp) :: f2
    real(dp) :: f3
    real(dp) :: f4
    real(dp) :: f5
    real(dp) :: f6
    real(dp), dimension(:), intent(in) :: x

    f1 = x(2) - x(1) * x(1)
    f2 = 1.0_dp - x(1)
    f3 = x(4) - x(3) * x(3)
    f4 = 1.0_dp - x(3)
    f5 = x(2) + x(4) - 2.0_dp
    f6 = x(2) - x(4)

    wood = 100.0_dp * f1 * f1 &
         +             f2 * f2 &
         +  90.0_dp * f3 * f3 &
         +             f4 * f4 &
         +  10.0_dp * f5 * f5 &
         +   0.1_dp * f6 * f6

  end function wood

  ! ------------------------------------------------------------------
  !
  ! The Chebyquad function, 1 <= N.
  ! Solution: n==2: x(1:2) = (/ 0.2113249_dp, 0.7886751_dp /)
  !           n==4: x(1:4) = (/ 0.1026728_dp, 0.4062037_dp, 0.5937963_dp, 0.8973272_dp /)
  !           n==6: x(1:6) = (/ 0.066877_dp, 0.288741_dp, 0.366682_dp,
  !                             0.633318_dp, 0.711259_dp, 0.933123_dp /)
  !           n==8: x(1:8) = (/ 0.043153_dp, 0.193091_dp, 0.266329_dp, 0.500000_dp, &
  !                             0.500000_dp, 0.733671_dp, 0.806910_dp, 0.956847_dp /)
  !           else unknown
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function chebyquad(x)

    implicit none

    integer(i4) :: n

    real(dp) :: chebyquad
    real(dp), dimension(:), intent(in) :: x

    real(dp), dimension(size(x)) :: fvec
    integer(i4) :: i
    integer(i4) :: j
    real(dp) :: t
    real(dp) :: t1
    real(dp) :: t2
    real(dp) :: th

    !
    !  Compute FVEC.
    !
    n = size(x)
    fvec(1:n) = 0.0_dp
    do j = 1, n
       t1 = 1.0_dp
       t2 = 2.0_dp * x(j) - 1.0_dp
       t = 2.0_dp * t2
       do i = 1, n
          fvec(i) = fvec(i) + t2
          th = t * t2 - t1
          t1 = t2
          t2 = th
       end do
    end do

    do i = 1, n
       fvec(i) = fvec(i) / real ( n, dp )
       if ( mod ( i, 2 ) == 0 ) then
          fvec(i) = fvec(i) + 1.0_dp / ( real ( i, dp )**2 - 1.0_dp )
       end if
    end do
    !call p18_fvec ( n, x, fvec )
    !
    !  Compute F.
    !
    chebyquad = sum ( fvec(1:n)**2 )

  end function chebyquad

  ! ------------------------------------------------------------------
  !
  ! The Leon''s cubic valley function, N = 2.
  ! Solution: x(1:2) = (/ 1.0_dp, 1.0_dp /)
  !
  !  Discussion:
  !
  !    The function is similar to Rosenbrock's.  The "valley" follows
  !    the curve X(2) = X(1)**3.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !    A Leon,
  !    A Comparison of Eight Known Optimizing Procedures,
  !    in Recent Advances in Optimization Techniques,
  !    edited by Abraham Lavi, Thomas Vogl,
  !    Wiley, 1966.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function leon_cubic_valley(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: leon_cubic_valley
    real(dp) :: f1
    real(dp) :: f2
    real(dp), dimension(:), intent(in) :: x

    f1 = x(2) - x(1) * x(1) * x(1)
    f2 = 1.0_dp - x(1)

    leon_cubic_valley = 100.0_dp * f1 * f1 &
         +             f2 * f2

  end function leon_cubic_valley

  ! ------------------------------------------------------------------
  !
  ! The Gregory and Karney''s Tridiagonal Matrix Function, 1 <= N.
  ! Solution: forall(i=1:n) x(i) = real(n+1-i, dp)
  !
  !  Discussion:
  !
  !    The function has the form
  !      f = x'*A*x - 2*x(1)
  !    where A is the (-1,2,-1) tridiagonal matrix, except that A(1,1) is 1.
  !    The minimum value of F(X) is -N, which occurs for
  !      x = ( n, n-1, ..., 2, 1 ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    20 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Prentice Hall, 1973,
  !    Reprinted by Dover, 2002.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function gregory_karney_tridia_matrix(x)

    implicit none

    integer(i4) :: n

    real(dp) :: gregory_karney_tridia_matrix
    integer(i4) ::i
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    gregory_karney_tridia_matrix = x(1) * x(1) + 2.0_dp * sum ( x(2:n)**2 )

    do i = 1, n-1
       gregory_karney_tridia_matrix = gregory_karney_tridia_matrix - 2.0_dp * x(i) * x(i+1)
    end do

    gregory_karney_tridia_matrix = gregory_karney_tridia_matrix - 2.0_dp * x(1)

    gregory_karney_tridia_matrix = gregory_karney_tridia_matrix + real(n,dp)

  end function gregory_karney_tridia_matrix

  ! ------------------------------------------------------------------
  !
  ! The Hilbert function, 1 <= N.
  ! Solution: x(1:n) = 0.0_dp
  !
  !  Discussion:
  !
  !    The function has the form
  !      f = x'*A*x
  !    where A is the Hilbert matrix.  The minimum value
  !    of F(X) is 0, which occurs for
  !      x = ( 0, 0, ..., 0 ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    20 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Richard Brent,
  !    Algorithms for Minimization with Derivatives,
  !    Dover, 2002,
  !    ISBN: 0-486-41998-3,
  !    LC: QA402.5.B74.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function hilbert(x)

    implicit none

    integer(i4) :: n

    real(dp) :: hilbert
    integer(i4) ::i
    integer(i4) ::j
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    hilbert = 0.0_dp

    do i = 1, n
       do j = 1, n
          hilbert = hilbert + x(i) * x(j) / real ( i + j - 1, dp )
       end do
    end do

  end function hilbert

  ! ------------------------------------------------------------------
  !
  ! The De Jong Function F1, N = 3.
  ! Solution: x(1:n) = 0.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function de_jong_f1(x)

    implicit none

    integer(i4) :: n

    real(dp) :: de_jong_f1
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    de_jong_f1 = sum ( x(1:n)**2 )

  end function de_jong_f1

  ! ------------------------------------------------------------------
  !
  ! The De Jong Function F2, N = 2.
  ! Solution: x(1:n) = 1.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function de_jong_f2(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: de_jong_f2
    real(dp), dimension(:), intent(in) :: x

    de_jong_f2 = 100.0_dp * ( x(1) * x(1) - x(2) )**2 + ( 1.0_dp - x(1) )**2

  end function de_jong_f2

  ! ------------------------------------------------------------------
  !
  ! The De Jong Function F3 (discontinuous), N = 5.
  ! Solution:
  !    x(1:n) depends on search space
  !           = sum of integer part of left boundary values
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function de_jong_f3(x)

    implicit none

    integer(i4) :: n

    real(dp) :: de_jong_f3
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    ! Original: conversion to int only possible up to i8, else "Floating invalid operation - aborting"
    ! de_jong_f3 = real ( sum ( int ( x(1:n) ) ), dp )

    de_jong_f3 = sum ( aint ( x(1:n) ) )

  end function de_jong_f3

  ! ------------------------------------------------------------------
  !
  ! The De Jong Function F4 (Gaussian noise), N = 30.
  ! Solution: x(1:n) = 0.0_dp
  !
  !  Discussion:
  !
  !    The function includes Gaussian noise, multiplied by a parameter P.
  !
  !    If P is zero, then the function is a proper function, and evaluating
  !    it twice with the same argument will yield the same results.
  !    Moreover, P25_G and P25_H are the correct gradient and hessian routines.
  !
  !    If P is nonzero, then evaluating the function at the same argument
  !    twice will generally yield two distinct values; this means the function
  !    is not even a well defined mathematical function, let alone continuous;
  !    the gradient and hessian are not correct.  And yet, at least for small
  !    values of P, it may be possible to approximate the minimizer of the
  !    (underlying well-defined ) function.
  !
  !    The value of the parameter P is by default 1.  The user can manipulate
  !    this value by calling P25_P_GET or P25_P_SET.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function de_jong_f4(x)

    implicit none

    integer(i4) :: n

    real(dp) :: de_jong_f4
    real(dp) :: gauss
    integer(i4) ::i
    real(dp) :: p
    real(dp), dimension(:), intent(in) :: x

    real(dp), save :: p_save = 1.0_dp

    n = size(x)
    p = p_save
    !call p25_p_get ( p )

    call normal_01_sample ( gauss )

    de_jong_f4 = p * gauss
    do i = 1, n
       de_jong_f4 = de_jong_f4 + real ( i, dp ) * x(i)**4
    end do

  end function de_jong_f4

  ! ------------------------------------------------------------------
  !
  ! The De Jong Function F5, N = 2.
  ! Solution: x(1:n) = 0.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function de_jong_f5(x)

    implicit none

    !    integer(i4) :: n

    integer(i4) ::a1
    integer(i4) ::a2
    real(dp) :: de_jong_f5
    real(dp) :: fi
    real(dp) :: fj
    integer(i4) ::j
    integer(i4) ::j1
    integer(i4) ::j2
    integer(i4), parameter :: jroot = 5
    integer(i4), parameter :: k = 500
    real(dp), dimension(:), intent(in) :: x

    fi = real(k,dp)

    do j=1, jroot * jroot

       j1 = mod(j-1_i4, jroot) + 1_i4
       a1 = -32_i4 + j1 * 16_i4

       j2 = (j-1_i4) / jroot
       a2 = -32_i4 + j2 * 16_i4

       fj = real(j,dp) + (x(1) - real(a1,dp))**6 + (x(2) - real(a2,dp))**6

       fi = fi + 1.0_dp / fj

    end do

    de_jong_f5 = 1.0_dp / fi

    de_jong_f5 = de_jong_f5 - 0.0019996667_dp

  end function de_jong_f5

  ! ------------------------------------------------------------------
  !
  ! The Schaffer Function F6, N = 2.
  ! Solution: x(1:n) = (/ 0.0_dp, 0.0_dp /)
  !
  !  Discussion:
  !
  !    F can be regarded as a function of R = SQRT ( X(1)^2 + X(2)^2 ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function schaffer_f6(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: a
    real(dp) :: b
    real(dp) :: schaffer_f6
    real(dp) :: r
    real(dp), dimension(:), intent(in) :: x

    r = sqrt ( x(1)**2 + x(2)**2 )

    a = ( 1.0_dp + 0.001_dp * r**2 )**( -2 )

    b = ( sin ( r ) )**2 - 0.5_dp

    schaffer_f6 = 0.5_dp + a * b

  end function schaffer_f6

  ! ------------------------------------------------------------------
  !
  ! The Schaffer Function F7, N = 2.
  ! Solution: x(1:n) = (/ 0.0_dp, 0.0_dp /)
  !
  !  Discussion:
  !
  !    Note that F is a function of R^2 = X(1)^2 + X(2)^2
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    08 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function schaffer_f7(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: schaffer_f7
    real(dp) :: r
    real(dp), dimension(:), intent(in) :: x

    r = sqrt ( x(1)**2 + x(2)**2 )

    schaffer_f7 = sqrt ( r ) * ( 1.0_dp + ( sin ( 50.0_dp * r**0.2_dp ) )**2 )

  end function schaffer_f7

  ! ------------------------------------------------------------------
  !
  ! The Goldstein Price Polynomial, N = 2.
  ! Solution:
  !    x(1:n) = (/ 0.0_dp, -1.0_dp /)
  !    f(x)   = 3.0_dp
  !    http://www.pg.gda.pl/~mkwies/dyd/geadocu/fcngold.html
  !
  !  Discussion:
  !
  !    Note that F is a polynomial in X.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    08 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function goldstein_price_polynomial(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: a
    real(dp) :: b
    real(dp) :: c
    real(dp) :: d
    real(dp) :: goldstein_price_polynomial
    real(dp), dimension(:), intent(in) :: x

    a = x(1) + x(2) + 1.0_dp

    b = 19.0_dp - 14.0_dp * x(1) + 3.0_dp * x(1) * x(1) - 14.0_dp * x(2) &
         + 6.0_dp * x(1) * x(2) + 3.0_dp * x(2) * x(2)

    c = 2.0_dp * x(1) - 3.0_dp * x(2)

    d = 18.0_dp - 32.0_dp * x(1) + 12.0_dp * x(1) * x(1) + 48.0_dp * x(2) &
         - 36.0_dp * x(1) * x(2) + 27.0_dp * x(2) * x(2)

    goldstein_price_polynomial = ( 1.0_dp + a * a * b ) * ( 30.0_dp + c * c * d )

    goldstein_price_polynomial = goldstein_price_polynomial - 3.0_dp

  end function goldstein_price_polynomial

  ! ------------------------------------------------------------------
  !
  ! The Branin RCOS Function, N = 2.
  ! Solution: 1st solution: x(1:n) = (/ -pi, 12.275_dp /)
  !           2nd solution: x(1:n) = (/  pi,  2.275_dp /)
  !           3rd solution: x(1:n) = (/ 9.42478_dp, 2.475_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function branin_rcos(x)

    use mo_constants, only: pi_dp

    implicit none

    !    integer(i4) :: n

    real(dp), parameter :: a = 1.0_dp
    real(dp) :: b
    real(dp) :: c
    real(dp), parameter :: d = 6.0_dp
    real(dp), parameter :: e = 10.0_dp
    real(dp) :: branin_rcos
    real(dp) :: ff
    real(dp), dimension(:), intent(in) :: x

    b = 5.1_dp / ( 4.0_dp * pi_dp**2 )
    c = 5.0_dp / pi_dp
    ff = 1.0_dp / ( 8.0_dp * pi_dp )

    branin_rcos = a * ( x(2) - b * x(1)**2 + c * x(1) - d )**2 &
         + e * ( 1.0_dp - ff ) * cos ( x(1) ) + e

    branin_rcos = branin_rcos - 0.3978873577_dp

  end function branin_rcos

  ! ------------------------------------------------------------------
  !
  ! The Shekel SQRN5 Function, N = 4.
  ! Solution: x(1:n) = (/ 4.0000371429_dp, 4.0001315700_dp, 4.0000379073_dp, 4.0001323857_dp /)
  !
  !  Discussion:
  !
  !    The minimal function value is -10.1527236935_dp.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function shekel_sqrn5(x)

    implicit none

    integer(i4), parameter :: m = 5
    integer(i4) :: n

    real(dp), parameter, dimension ( 4, m ) :: a = reshape ( &
         (/ 4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, &
         1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
         8.0_dp, 8.0_dp, 8.0_dp, 8.0_dp, &
         6.0_dp, 6.0_dp, 6.0_dp, 6.0_dp, &
         3.0_dp, 7.0_dp, 3.0_dp, 7.0_dp /), (/ 4, m /) )
    real(dp), save, dimension ( m ) :: c = &
         (/ 0.1_dp, 0.2_dp, 0.2_dp, 0.4_dp, 0.6_dp /)
    real(dp) :: shekel_sqrn5
    integer(i4) ::j
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    shekel_sqrn5 = 0.0_dp
    do j = 1, m
       shekel_sqrn5 = shekel_sqrn5 - 1.0_dp / ( c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 ) )
    end do

    shekel_sqrn5 = shekel_sqrn5 + 10.1527236935_dp

  end function shekel_sqrn5

  ! ------------------------------------------------------------------
  !
  ! The Shekel SQRN7 Function, N = 4.
  ! Solution: x(1:n) = (/ 4.0005729560_dp, 4.0006881764_dp, 3.9994902225_dp, 3.9996048794_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function shekel_sqrn7(x)

    implicit none

    integer(i4), parameter :: m = 7
    integer(i4) :: n

    real(dp), parameter, dimension ( 4, m ) :: a = reshape ( &
         (/ 4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, &
         1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
         8.0_dp, 8.0_dp, 8.0_dp, 8.0_dp, &
         6.0_dp, 6.0_dp, 6.0_dp, 6.0_dp, &
         3.0_dp, 7.0_dp, 3.0_dp, 7.0_dp, &
         2.0_dp, 9.0_dp, 2.0_dp, 9.0_dp, &
         5.0_dp, 5.0_dp, 3.0_dp, 3.0_dp /), (/ 4, m /) )
    real(dp), save, dimension ( m ) :: c = &
         (/ 0.1_dp, 0.2_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.6_dp, 0.3_dp /)
    real(dp) :: shekel_sqrn7
    integer(i4) ::j
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    shekel_sqrn7 = 0.0_dp
    do j = 1, m
       shekel_sqrn7 = shekel_sqrn7 - 1.0_dp / ( c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 ) )
    end do

    shekel_sqrn7 = shekel_sqrn7 + 10.4024645722_dp

  end function shekel_sqrn7

  ! ------------------------------------------------------------------
  !
  ! The Shekel SQRN10 Function, N = 4.
  ! Solution: x(1:n) = (/ 4.0007465727_dp, 4.0005916919_dp, 3.9996634360_dp, 3.9995095935_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function shekel_sqrn10(x)

    implicit none

    integer(i4), parameter :: m = 10
    integer(i4) :: n

    real(dp), parameter, dimension ( 4, m ) :: a = reshape ( &
         (/ 4.0, 4.0, 4.0, 4.0, &
         1.0, 1.0, 1.0, 1.0, &
         8.0, 8.0, 8.0, 8.0, &
         6.0, 6.0, 6.0, 6.0, &
         3.0, 7.0, 3.0, 7.0, &
         2.0, 9.0, 2.0, 9.0, &
         5.0, 5.0, 3.0, 3.0, &
         8.0, 1.0, 8.0, 1.0, &
         6.0, 2.0, 6.0, 2.0, &
         7.0, 3.6, 7.0, 3.6 /), (/ 4, m /) )

    real(dp), save, dimension ( m ) :: c = &
         (/ 0.1_dp, 0.2_dp, 0.2_dp, 0.4_dp, 0.6_dp, &
         0.6_dp, 0.3_dp, 0.7_dp, 0.5_dp, 0.5_dp /)
    real(dp) :: shekel_sqrn10
    integer(i4) ::j
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    shekel_sqrn10 = 0.0_dp
    do j = 1, m
       shekel_sqrn10 = shekel_sqrn10 - 1.0_dp / ( c(j) + sum ( ( x(1:n) - a(1:n,j) )**2 ) )
    end do

    shekel_sqrn10 = shekel_sqrn10 + 10.5359339075_dp

  end function shekel_sqrn10

  ! ------------------------------------------------------------------
  !
  ! The Six-Hump Camel-Back Polynomial, N = 2.
  ! Solution: 1st solution: x(1:n) = (/ -0.0898_dp,  0.7126_dp /)
  !           2nd solution: x(1:n) = (/  0.0898_dp, -0.7126_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function six_hump_camel_back_polynomial(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: six_hump_camel_back_polynomial
    real(dp), dimension(:), intent(in) :: x

    six_hump_camel_back_polynomial = ( 4.0_dp - 2.1_dp * x(1)**2 + x(1)**4 / 3.0_dp ) * x(1)**2 &
         + x(1) * x(2) + 4.0_dp * ( x(2)**2 - 1.0_dp ) * x(2)**2

    six_hump_camel_back_polynomial = six_hump_camel_back_polynomial + 1.0316284229_dp

  end function six_hump_camel_back_polynomial

  ! ------------------------------------------------------------------
  !
  ! The Shubert Function, N = 2.
  ! Solution: x(1:n) = (/ -7.0835064076515595_dp, -7.708313735499347_dp /)
  !
  !  Discussion:
  !
  !    For -10 <= X(I) <= 10, the function has 760 local minima, 18 of which
  !    are global minima, with minimum value -186.73090883102378.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !    Bruno Shubert,
  !    A sequential method seeking the global maximum of a function,
  !    SIAM Journal on Numerical Analysis,
  !    Volume 9, pages 379-388, 1972.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function shubert(x)

    implicit none

    integer(i4) :: n

    real(dp) :: shubert
    real(dp) :: factor
    integer(i4) ::i
    integer(i4) ::k
    real(dp) :: k_r8
    real(dp), dimension(:), intent(in) :: x

    n = size(x)
    shubert = 1.0_dp

    do i = 1, n
       factor = 0.0_dp
       do k = 1, 5
          k_r8 = real ( k, dp )
          factor = factor + k_r8 * cos ( ( k_r8 + 1.0_dp ) * x(i) + k_r8 )
       end do
       shubert = shubert * factor
    end do

    shubert = shubert + 186.7309088310_dp

  end function shubert

  ! ------------------------------------------------------------------
  !
  ! The Stuckman Function, N = 2.
  ! Solution: Only iterative solution; Check reference.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 October 2004
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function stuckman(x)

    use mo_utils, only: eq, le

    implicit none

    !    integer(i4) :: n

    real(dp) :: a1
    real(dp) :: a2
    real(dp) :: b
    real(dp) :: stuckman
    real(dp) :: m1
    real(dp) :: m2
    real(dp) :: r11
    real(dp) :: r12
    real(dp) :: r21
    real(dp) :: r22
    real(dp), dimension(:), intent(in) :: x

    real(dp), save :: b_save = 0.0_dp
    real(dp), save :: m1_save = 0.0_dp
    real(dp), save :: m2_save = 0.0_dp
    real(dp), save :: r11_save = 0.0_dp
    real(dp), save :: r12_save = 0.0_dp
    real(dp), save :: r21_save = 0.0_dp
    real(dp), save :: r22_save = 0.0_dp

    !call p36_p_get ( b, m1, m2, r11, r12, r21, r22, seed )
    b = b_save
    m1 = m1_save
    m2 = m2_save
    r11 = r11_save
    r12 = r12_save
    r21 = r21_save
    r22 = r22_save

    a1 = r8_aint ( abs ( x(1) - r11 ) ) + r8_aint ( abs ( x(2) - r21 ) )
    a2 = r8_aint ( abs ( x(1) - r12 ) ) + r8_aint ( abs ( x(2) - r22 ) )

    if ( le(x(1),b) ) then
       if ( eq(a1,0.0_dp) ) then
          stuckman = r8_aint ( m1 )
       else
          stuckman = r8_aint ( m1 * sin ( a1 ) / a1 )
       end if
    else
       if ( eq(a2,0.0_dp) ) then
          stuckman = r8_aint ( m2 )
       else
          stuckman = r8_aint ( m2 * sin ( a2 ) / a2 )
       end if
    end if

  end function stuckman

  ! ------------------------------------------------------------------
  !
  ! The Easom Function, N = 2.
  ! Solution: x(1:n) = (/ pi, pi /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function easom(x)

    use mo_constants, only: pi_dp

    implicit none

    !    integer(i4) :: n

    real(dp) :: arg
    real(dp) :: easom
    real(dp), dimension(:), intent(in) :: x

    arg = - ( x(1) - pi_dp )**2 - ( x(2) - pi_dp )**2
    easom = - cos ( x(1) ) * cos ( x(2) ) * exp ( arg )
    easom = easom + 1.0_dp

  end function easom

  ! ------------------------------------------------------------------
  !
  ! The Bohachevsky Function #1, N = 2.
  ! Solution: x(1:n) = (/ 0.0_dp, 0.0_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function bohachevsky1(x)

    use mo_constants, only: pi_dp

    implicit none

    !    integer(i4) :: n

    real(dp) :: bohachevsky1
    real(dp), dimension(:), intent(in) :: x

    bohachevsky1 =           x(1) * x(1) - 0.3_dp * cos ( 3.0_dp * pi_dp * x(1) ) &
         + 2.0_dp * x(2) * x(2) - 0.4_dp * cos ( 4.0_dp * pi_dp * x(2) ) &
         + 0.7_dp

  end function bohachevsky1

  ! ------------------------------------------------------------------
  !
  ! The Bohachevsky Function #2, N = 2.
  ! Solution: x(1:n) = (/ 0.0_dp, 0.0_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function bohachevsky2(x)

    use mo_constants, only: pi_dp

    implicit none

    !    integer(i4) :: n

    real(dp) :: bohachevsky2
    real(dp), dimension(:), intent(in) :: x

    bohachevsky2 = x(1) * x(1) + 2.0_dp * x(2) * x(2) &
         - 0.3_dp * cos ( 3.0_dp * pi_dp * x(1) ) &
         * cos ( 4.0_dp * pi_dp * x(2) ) + 0.3_dp

  end function bohachevsky2

  ! ------------------------------------------------------------------
  !
  ! The Bohachevsky Function #3, N = 2.
  ! Solution:
  !     x(1:2) = (/ 0.0_dp, 0.0_dp /)
  !     f(x)   = 0.0_dp
  !
  !  Discussion:
  !    J. Burkhardt:
  !       There is a typo in the reference.  I'm just guessing at the correction.
  !    J. Mai:
  !       Typo in function found.
  !       see: http://www-optima.amp.i.kyoto-u.ac.jp/member/student/
  !            hedar/Hedar_files/TestGO_files/Page595.htm
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function bohachevsky3(x)

    use mo_constants, only: pi_dp

    implicit none

    !    integer(i4) :: n

    real(dp) :: bohachevsky3
    real(dp), dimension(:), intent(in) :: x

    bohachevsky3 = x(1)**2 + 2.0_dp * x(2)**2 &
         - 0.3_dp * cos ( 3.0_dp * pi_dp * x(1)  &
         + 4.0_dp * pi_dp * x(2) ) + 0.3_dp

  end function bohachevsky3

  ! ------------------------------------------------------------------
  !
  ! The Colville Polynomial, N = 4.
  ! Solution: x(1:n) = (/ 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 January 2001
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Zbigniew Michalewicz,
  !    Genetic Algorithms + Data Structures = Evolution Programs,
  !    Third Edition,
  !    Springer Verlag, 1996,
  !    ISBN: 3-540-60676-9,
  !    LC: QA76.618.M53.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function colville_polynomial(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: colville_polynomial
    real(dp), dimension(:), intent(in) :: x

    colville_polynomial = 100.0_dp * ( x(2) - x(1)**2 )**2 &
         + ( 1.0_dp - x(1) )**2 &
         + 90.0_dp * ( x(4) - x(3)**2 )**2 &
         + ( 1.0_dp - x(3) )**2 &
         + 10.1_dp * ( ( x(2) - 1.0_dp )**2 + ( x(4) - 1.0_dp )**2 ) &
         + 19.8_dp * ( x(2) - 1.0_dp ) * ( x(4) - 1.0_dp )

  end function colville_polynomial

  ! ------------------------------------------------------------------
  !
  ! The Powell 3D function, N = 3.
  ! Solution: x(1:n) = (/ 1.0_dp, 1.0_dp, 1.0_dp /)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 March 2002
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    MJD Powell,
  !    An Efficient Method for Finding the Minimum of a Function of
  !    Several Variables Without Calculating Derivatives,
  !    Computer Journal,
  !    Volume 7, Number 2, pages 155-162, 1964.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function powell3d(x)

    use mo_constants, only: pi_dp
    use mo_utils, only: eq

    implicit none

    !    integer(i4) :: n

    real(dp) :: arg
    real(dp) :: powell3d
    real(dp) :: term
    real(dp), dimension(:), intent(in) :: x

    if ( eq(x(2),0.0_dp) ) then
       term = 0.0_dp
    else
       !arg = ( x(1) + 2.0_dp * x(2) + x(3) ) / x(2)
       ! changed according to original paper of Powell (1964)
       arg = ( x(1) + x(3) ) / x(2) - 2.0_dp

       term = arg*arg
       if ( term .lt. 708._dp ) then
          term = exp ( - term )
       else
          ! avoid underflow
          term = 0.0_dp
       end if
    end if

    powell3d = 3.0_dp &
         - 1.0_dp / ( 1.0_dp + ( x(1) - x(2) )**2 ) &
         - sin ( 0.5_dp * pi_dp * x(2) * x(3) ) &
         - term

  end function powell3d

  ! ------------------------------------------------------------------
  !
  ! The Himmelblau function, N = 2.
  ! Solution: x(1:2) = (/ 3.0_dp, 2.0_dp /)
  !
  !  Discussion:
  !
  !    This function has 4 global minima:
  !
  !      X* = (  3,        2       ), F(X*) = 0.
  !      X* = (  3.58439, -1.84813 ), F(X*) = 0.
  !      X* = ( -3.77934, -3.28317 ), F(X*) = 0.
  !      X* = ( -2.80512,  3.13134 ), F(X*) = 0.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 January 2008
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    David Himmelblau,
  !    Applied Nonlinear Programming,
  !    McGraw Hill, 1972,
  !    ISBN13: 978-0070289215,
  !   LC: T57.8.H55.
  !
  !  Parameters:
  !
  !    Input, real(dp) :: X(N), the argument of the objective function.
  !

  function himmelblau(x)

    implicit none

    !    integer(i4) :: n

    real(dp) :: himmelblau
    real(dp), dimension(:), intent(in) :: x

    himmelblau = ( x(1)**2 + x(2) - 11.0_dp )**2 &
         + ( x(1) + x(2)**2 - 7.0_dp )**2

  end function himmelblau

  ! ------------------------------------------------------------------
  !
  ! The Griewank Function, N=2 or N=10.
  ! Solution: x(1:n) = 0.0_dp

  !
  ! Coded originally by Q Duan. Edited for incorporation into Fortran DDS algorithm by
  ! Bryan Tolson, Nov 2005.
  !
  ! Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  ! I/O Variable definitions:
  !     nopt     -  the number of decision variables
  !     x_values -      an array of decision variable values (size nopt)
  !     fvalue   -      the value of the objective function with x_values as input

  function griewank(x_values)

    use mo_kind, only: i4, dp

    implicit none

    real(dp), dimension(:), intent(in)  :: x_values
    real(dp) :: griewank

    integer(i4) :: nopt
    integer(i4) :: j
    real(dp)    :: d, u1, u2

    nopt = size(x_values)
    if (nopt .eq. 2) then
       d = 200.0_dp
    else
       d = 4000.0_dp
    end if
    u1 = sum(x_values**2) / d
    u2 = 1.0_dp
    do j=1, nopt
       u2 = u2 * cos(x_values(j)/sqrt(real(j,dp)))
    end do
    griewank = u1 - u2 + 1.0_dp
    !
  end function griewank

  ! ------------------------------------------------------------------
  !
  !  The Rosenbrock parabolic value function, N = 2.
  !  Solution: x(1:n) = 1.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 February 2008
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    R ONeill,
  !    Algorithm AS 47:
  !    Function Minimization Using a Simplex Procedure,
  !    Applied Statistics,
  !    Volume 20, Number 3, 1971, pages 338-345.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(2), the argument.
  !

  function rosenbrock(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: rosenbrock

    real(dp) :: fx
    real(dp) :: fx1
    real(dp) :: fx2

    fx1 = x(2) - x(1) * x(1)
    fx2 = 1.0_dp - x(1)

    fx = 100.0_dp * fx1 * fx1 + fx2 * fx2

    rosenbrock = fx

  end function rosenbrock

  ! ------------------------------------------------------------------
  !
  !  The Sphere model, N >= 1.
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function sphere_model(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: sphere_model

    integer(i4) :: n
    integer(i4) :: j

    n = size(x)
    sphere_model = 0.0_dp
    do j=1, n
       sphere_model = sphere_model + x(j)**2
    enddo

  end function sphere_model

  ! ------------------------------------------------------------------
  !
  !  The Rastrigin function, N >= 2.
  !  Solution: x(1:n) = 0.0_dp
  !  Search domain: -5.12 <= xi <= 5.12
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function rastrigin(x)

    use mo_constants, only: pi_dp

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: rastrigin

    integer(i4) :: n
    integer(i4) :: j

    n = size(x,dim=1)
    rastrigin = 0.0_dp
    do j=1, n
       rastrigin = rastrigin+ (x(j)**2 - 10.0_dp*cos(2.0_dp*pi_dp*x(j)))
    enddo
    rastrigin = rastrigin + 10.0_dp * real(n,dp)

    ! FUNCTN04 of SCEUA F77 source code
    !   Bound: X1=[-1,1], X2=[-1,1]
    !   Global Optimum: -2, (0,0)
    ! n = size(x,dim=1)
    ! rastrigin = 0.0_dp
    ! do j=1, n
    !    rastrigin = rastrigin+ (x(j)**2 - cos(18.0_dp*x(j)))
    ! enddo

  end function rastrigin

  ! ------------------------------------------------------------------
  !
  !  The Schwefel function, N >= 2.
  !  Solution: x(1:n) = 1.0_dp
  !  Solution: x(1:n) = 420.968746_dp   ( see (2) and (3) )
  !
  !  Author:
  !
  !    (1) Matlab code by A. Hedar
  !    (2) http://www.aridolan.com/ga/gaa/Schwefel.html
  !    (3) http://www.it.lut.fi/ip/evo/functions/node10.html
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function schwefel(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: schwefel

    integer(i4) :: n

    n = size(x)
    schwefel = sum(-x*sin(sqrt(abs(x)))) + 418.982887272433799807913601398_dp*real(n,dp)

  end function schwefel

  ! ------------------------------------------------------------------
  !
  !  The Ackley function, N >= 2.
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function ackley(x)

    use mo_constants, only: pi_dp

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: ackley

    integer(i4) :: n
    real(dp), parameter :: a = 20.0_dp
    real(dp), parameter :: b = 0.2_dp
    real(dp), parameter :: c = 2.0_dp*pi_dp
    real(dp) :: s1, s2

    n = size(x)
    s1 = sum(x**2)
    s2 = sum(cos(c*x))
    ackley = -a * exp(-b*sqrt(1.0_dp/real(n,dp)*s1)) - exp(1.0_dp/real(n,dp)*s2) + a + exp(1.0_dp)

  end function ackley

  ! ------------------------------------------------------------------
  !
  !  The Michalewicz function, N >= 2.
  !  Search domain: x restricted to (0, Pi)
  !  Solution:
  !     numerical, so far best found
  !     x(1:2)  = (/ 2.2029262967_dp, 1.5707721052_dp /)
  !     f(x)    = -1.8013033793_dp
  !
  !     numerical, so far best found
  !     x(1:5)  = (/ 2.2029054902_dp, 1.5707963436_dp, 1.2849915892_dp, 1.9230584622_dp, 1.7204697668_dp /)
  !     f(x)    = -4.6876581791_dp
  !     known from literature:
  !     f(x)    = -4.687_dp
  !
  !     x(1:10) = (/ 2.2029055201726093_dp, 2.10505573543129_dp, &
  !                  2.2193332517481035_dp, 1.9230584698663626_dp, &
  !                  0.9966770382174085_dp, 2.0274797779024762_dp, &
  !                  1.7114837946034247_dp, 1.3605717365168617_dp, &
  !                  1.2828240065709524_dp, 1.5707963267948966_dp /)
  !     f(x)    = -7.209069703423156_dp
  !     known from literature:
  !     f(x)    = -9.66_dp
  !
  !  Discussion:
  !    The Michalewicz function is a multimodal test function (n! local optima).
  !    The parameter p defines the "steepness" of the valleys or edges. Larger p leads to
  !    more difficult search. For very large p the function behaves like a needle in the
  !    haystack (the function values for points in the space outside the narrow peaks give
  !    very little information on the location of the global optimum).
  !    http://www.geatbx.com/docu/fcnindex-01.html#P150_6749
  !
  !    http://www.scribd.com/doc/74351406/12/Michalewicz''s-function
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function michalewicz(x)

    use mo_constants, only: pi_dp

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp)                           :: michalewicz

    integer(i4) :: n
    integer(i4) :: j
    real(dp)    :: tmp
    integer(i4), parameter :: p = 20

    n = size(x)
    michalewicz = 0.0_dp
    do j=1, n
       ! michalewicz = michalewicz + sin(x(j)) * sin(x(j)**2 * (real(j,dp)/pi_dp))**p
       tmp = x(j)*x(j) * (real(j,dp)/pi_dp)
       tmp = sin(tmp)
       if (abs(tmp) .lt. 1E-15) then
          tmp = 0.0_dp
       else
          tmp = tmp**p
       end if
       tmp = sin(x(j)) * tmp
       michalewicz = michalewicz + tmp
    end do
    michalewicz = -michalewicz

    select case(n)
    case(2_i4)
       michalewicz = michalewicz + 1.8013033793_dp
    case(5_i4)
       michalewicz = michalewicz + 4.6876581791_dp
    end select

  end function michalewicz

  ! ------------------------------------------------------------------
  !
  !  The Booth function, N = 2.
  !  Solution: x(1:2) = (/ 1.0_dp, 3.0_dp /)
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function booth(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: booth

    booth = (x(1)+2.0_dp*x(2)-7.0_dp)**2 + (2.0_dp*x(1)+x(2)-5.0_dp)**2

  end function booth

  ! ------------------------------------------------------------------
  !
  !  The Hump function, N = 2.
  !
  !  Search Domain:
  !     -5.0_dp <= xi <= 5.0_dp
  !
  !  Solution:
  !     x(1:2) = (/ 0.08984201310031806_dp , -0.7126564030207396_dp /)     OR
  !     x(1:2) = (/ -0.08984201310031806_dp , 0.7126564030207396_dp /)
  !     f(x)   = 0.0_dp
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Literature:
  !
  !    http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1621.htm
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function hump(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: hump

    hump = 1.0316285_dp + 4.0_dp*x(1)**2 - 2.1_dp*x(1)**4 + x(1)**6 / 3.0_dp + x(1)*x(2) - 4.0_dp*x(2)**2 + 4.0_dp*x(2)**4

  end function hump

  ! ------------------------------------------------------------------
  !
  !  The Levy function, N >= 2.
  !  Solution: x(1:n) = 1.0_dp
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function levy(x)

    use mo_constants, only: pi_dp

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: levy

    integer(i4) :: n
    integer(i4) :: i
    real(dp), dimension(size(x)) :: z

    n = size(x)
    z = 1.0_dp+(x-1.0_dp)/4.0_dp
    levy = sin(pi_dp*z(1))**2
    do i=1, n-1
       levy = levy + (z(i)-1.0_dp)**2 * (1.0_dp+10.0_dp*(sin(pi_dp*z(i)+1.0_dp))**2)
    end do
    levy = levy + (z(n)-1.0_dp)**2 * (1.0_dp+(sin(2.0_dp*pi_dp*z(n)))**2)

  end function levy

  ! ------------------------------------------------------------------
  !
  !  The Matyas function, N = 2.
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function matyas(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: matyas

    matyas = 0.26_dp * (x(1)**2 + x(2)**2) -0.48_dp*x(1)*x(2)

  end function matyas

  ! ------------------------------------------------------------------
  !
  !  The Perm function, N >= 4.
  !  Solution: forall(i=1:n) x(i) = real(i,dp)
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function perm(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: perm

    integer(i4) :: n
    integer(i4) :: j, k
    real(dp) :: s_in

    n = size(x)
    perm = 0.0_dp
    do k=1, n
       s_in = 0.0_dp
       do j=1, n
          s_in = s_in + (real(j,dp)**k + 0.5_dp) * ((x(j)/real(j,dp))**k - 1.0_dp)
       enddo
       perm = perm + s_in**2
    enddo


  end function perm

  ! ------------------------------------------------------------------
  !
  !  The Power sum function, N = 4.
  !  Solution:
  !      x    = (/ 1.0000653551_dp, 2.0087089520_dp, 1.9912589253_dp, 2.9999732609_dp /)
  !      f(x) = 0.0000000001
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function power_sum(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: power_sum

    integer(i4) :: n
    integer(i4) :: k
    real(dp), dimension(4) :: b

    n = size(x)
    b = (/ 8.0_dp, 18.0_dp, 44.0_dp, 114.0_dp /)
    power_sum = 0.0_dp
    do k=1, n
       power_sum = power_sum + (sum(x**k) - b(k))**2
    end do

  end function power_sum

  ! ------------------------------------------------------------------
  !
  !  Sphere model, N = 1.
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Author:
  !
  !    Matlab code by A. Hedar
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Parameters:
  !
  !    Input, real(dp) X(N), the argument.
  !

  function sphere(x)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: sphere

    integer(i4) :: n
    integer(i4) :: j

    n = size(x)
    sphere = 0.0_dp
    do j=1, n
       sphere = sphere + x(j)**2
    enddo

  end function sphere

  ! ------------------------------------------------------------------
  !
  !  The sphere model
  !  N = 2
  !  Solution: x(1:n) = 1.0_dp
  !
  !  Discussion:
  !
  !    The function is continuous, convex, and unimodal.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 January 2012
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Hugues Bersini, Marco Dorigo, Stefan Langerman, Gregory Seront,
  !    Luca Gambardella,
  !    Results of the first international contest on evolutionary optimisation,
  !    In Proceedings of 1996 IEEE International Conference on Evolutionary
  !    Computation,
  !    IEEE Press, pages 611-615, 1996.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function sphere_model_2d(x)

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: sphere_model_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: j

    m = size(x,1)
    n = size(x,2)
    do j = 1, n
       sphere_model_2d(j) = sum( ( x(1:m,j) - 1.0_dp ) ** 2 )
    end do

  end function sphere_model_2d

  ! ------------------------------------------------------------------
  !
  !  The axis-parallel hyper-ellipsoid function
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Discussion:
  !
  !    This function is also known as the weighted sphere model.
  !
  !    The function is continuous, convex, and unimodal.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 December 2011
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Marcin Molga, Czeslaw Smutnicki,
  !    Test functions for optimization needs.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function axis_parallel_hyper_ellips_2d(x)

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: axis_parallel_hyper_ellips_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: j
    real(dp) :: y(size(x,1))


    m = size(x,1)
    n = size(x,2)
    forall(j=1:m) y(j) = real(j,dp)

    do j = 1, n
       axis_parallel_hyper_ellips_2d(j) = sum( y(1:m) * x(1:m,j) ** 2 )
    end do

  end function axis_parallel_hyper_ellips_2d

  ! ------------------------------------------------------------------
  !
  !  The rotated hyper-ellipsoid function
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Discussion:
  !
  !    This function is also known as the weighted sphere model.
  !
  !    The function is continuous, convex, and unimodal.
  !
  !     There is a typographical error in Molga and Smutnicki, so that the
  !     formula for this function is given incorrectly.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 December 2011
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Marcin Molga, Czeslaw Smutnicki,
  !    Test functions for optimization needs.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function rotated_hyper_ellipsoid_2d(x)

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: rotated_hyper_ellipsoid_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: i
    integer(i4) :: j
    real(dp) :: x_sum

    m = size(x,1)
    n = size(x,2)
    do j = 1, n

       rotated_hyper_ellipsoid_2d(j) = 0.0_dp
       x_sum = 0.0_dp

       do i = 1, m
          x_sum = x_sum + x(i,j)
          rotated_hyper_ellipsoid_2d(j) = rotated_hyper_ellipsoid_2d(j) + x_sum**2
       end do

    end do

  end function rotated_hyper_ellipsoid_2d

  ! ------------------------------------------------------------------
  !
  !  Rosenbrock''s valley
  !  Solution: x(1:n) = 1.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 January 2012
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Howard Rosenbrock,
  !    An Automatic Method for Finding the Greatest or Least Value of a Function,
  !    Computer Journal,
  !    Volume 3, 1960, pages 175-184.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function rosenbrock_2d(x)

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: rosenbrock_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: j

    m = size(x,1)
    n = size(x,2)
    do j = 1, n
       rosenbrock_2d(j) = sum ( ( 1.0_dp - x(1:m,j) )**2 ) &
            + sum ( ( x(2:m,j) - x(1:m-1,j) )**2 )
    end do

  end function rosenbrock_2d

  ! ------------------------------------------------------------------
  !
  !  Rastrigin''s function
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 December 2011
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Marcin Molga, Czeslaw Smutnicki,
  !    Test functions for optimization needs.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function rastrigin_2d(x)

    use mo_constants, only: pi_dp

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: rastrigin_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: i
    integer(i4) :: j

    m = size(x,1)
    n = size(x,2)
    do j = 1, n

       rastrigin_2d(j) = real( 10 * m, dp )

       do i = 1, m
          rastrigin_2d(j) = rastrigin_2d(j) + x(i,j) ** 2 - 10.0_dp * cos( 2.0_dp * pi_dp * x(i,j) )
       end do

    end do

  end function rastrigin_2d

  ! ------------------------------------------------------------------
  !
  !  Schwefel''s function.
  !  Solution: x(1:n) = 420.9687_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 December 2011
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Hans-Paul Schwefel,
  !    Numerical optimization of computer models,
  !    Wiley, 1981,
  !    ISBN13: 978-0471099888,
  !    LC: QA402.5.S3813.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function schwefel_2d(x)

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: schwefel_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: j

    m = size(x,1)
    n = size(x,2)
    do j = 1, n
       schwefel_2d(j) = -sum( x(1:m,j) * sin( sqrt( abs( x(1:m,j) ) ) ) )
    end do

  end function schwefel_2d

  ! ------------------------------------------------------------------
  !
  !  Griewank''s function
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 December 2011
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Marcin Molga, Czeslaw Smutnicki,
  !    Test functions for optimization needs.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function griewank_2d(x)

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: griewank_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: j
    real(dp) :: y(size(x,1))

    m = size(x,1)
    n = size(x,2)
    forall(j=1:m) y(j) = real(j,dp)
    y(1:m) = sqrt( y(1:m) )

    do j = 1, n
       griewank_2d(j) = sum( x(1:m,j) ** 2 ) / 4000.0_dp &
            - product( cos( x(1:m,j) / y(1:m) ) ) + 1.0_dp
    end do

  end function griewank_2d

  ! ------------------------------------------------------------------
  !
  !  The power sum function.
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 December 2011
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Marcin Molga, Czeslaw Smutnicki,
  !    Test functions for optimization needs.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function power_sum_2d(x)

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: power_sum_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: j
    real(dp) :: y(size(x,1))

    m = size(x,1)
    n = size(x,2)
    forall(j=1:m) y(j) = real(j,dp)
    y(1:m) = y(1:m) + 1.0_dp

    do j = 1, n
       power_sum_2d(j) = sum( abs( x(1:m,j) ) ** y(1:m) )
    end do

  end function power_sum_2d

  ! ------------------------------------------------------------------
  !
  !  Ackley''s function
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 December 2011
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Marcin Molga, Czeslaw Smutnicki,
  !    Test functions for optimization needs.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function ackley_2d(x)

    use mo_constants, only: pi_dp

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: ackley_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: j
    real(dp), parameter :: a = 20.0_dp
    real(dp), parameter :: b = 0.2_dp
    real(dp), parameter :: c = 0.2_dp

    m = size(x,1)
    n = size(x,2)
    do j = 1, n
       ackley_2d(j) = -a * exp( -b * sqrt( sum( x(1:m,j)**2 ) &
            / real( m, dp ) ) ) &
            - exp( sum( cos( c * pi_dp * x(1:m,j) ) ) / real( m, dp ) ) &
            + a + exp( 1.0_dp )
    end do

  end function ackley_2d

  ! ------------------------------------------------------------------
  !
  !  Michalewicz''s function
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 December 2011
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Marcin Molga, Czeslaw Smutnicki,
  !    Test functions for optimization needs.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function michalewicz_2d(x)

    use mo_constants, only: pi_dp

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: michalewicz_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: j
    integer(i4), parameter :: p = 10
    real(dp) :: y(size(x,1))

    m = size(x,1)
    n = size(x,2)
    forall(j=1:m) y(j) = real(j,dp)

    do j = 1, n
       michalewicz_2d(j) = -sum( &
            sin( x(1:m,j) ) * ( sin( x(1:m,j)**2 * y(1:m) / pi_dp ) ) ** ( 2 * p ) &
            )
    end do

  end function michalewicz_2d

  ! ------------------------------------------------------------------
  !
  !  The drop wave function
  !  Solution: x(1:n) = 0.0_dp
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 January 2012
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Marcin Molga, Czeslaw Smutnicki,
  !    Test functions for optimization needs.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function drop_wave_2d(x)

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: drop_wave_2d

    integer(i4) :: m
    integer(i4) :: n
    integer(i4) :: j
    real(dp) :: rsq

    m = size(x,1)
    n = size(x,2)
    do j = 1, n

       rsq = sum( x(1:m,j)**2 )

       drop_wave_2d(j) = -( 1.0_dp + cos( 12.0_dp * sqrt( rsq ) ) ) &
            / ( 0.5_dp * rsq + 2.0_dp )

    end do

  end function drop_wave_2d

  ! ------------------------------------------------------------------
  !
  !  The deceptive function
  !  Solution: forall(i=1:n) x(i) = real(i,dp)/real(n+1,dp)
  !
  !  Discussion:
  !
  !    In dimension I, the function is a piecewise linear function with
  !    local minima at 0 and 1.0, and a global minimum at ALPHA(I) = I/(M+1).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 December 2011
  !
  !  Author:
  !
  !    John Burkardt
  !    Modified Jul 2012 Matthias Cuntz - function, dp, etc.
  !
  !  Reference:
  !
  !    Marcin Molga, Czeslaw Smutnicki,
  !    Test functions for optimization needs.
  !
  !  Parameters:
  !
  !    Input, real(dp), dimension(:,:) :: x, the arguments size (m,n), m spatial dimension, n number of arguments.
  !
  !    Output, real(dp), dimension(size(x,2)) :: f, the function evaluated at the arguments.
  !

  function deceptive_2d(x)

    use mo_utils, only: le

    implicit none

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,2)) :: deceptive_2d

    integer(i4) :: m
    integer(i4) :: n
    real(dp) :: g
    integer(i4) :: i
    integer(i4) :: j
    real(dp) :: alpha(size(x,1))
    real(dp), parameter :: beta = 2.0_dp
    !
    !  I'm just choosing ALPHA in [0,1] arbitrarily.
    !
    m = size(x,1)
    n = size(x,2)
    do i = 1, m
       alpha(i) = real( i, dp ) / real( m + 1, dp )
    end do

    do j = 1, n

       deceptive_2d(j) = 0.0_dp

       do i = 1, m

          if ( le(x(i,j),0.0_dp) ) then
             g = x(i,j)
          else if ( le(x(i,j),0.8_dp * alpha(i)) ) then
             g = 0.8_dp - x(i,j) / alpha(i)
          else if ( le(x(i,j),alpha(i)) ) then
             g = 5.0_dp * x(i,j) / alpha(i) - 4.0_dp
          else if ( le(x(i,j),( 1.0_dp + 4.0_dp * alpha(i) ) / 5.0_dp) ) then
             g = 1.0_dp + 5.0_dp * ( x(i,j) - alpha(i) ) / ( alpha(i) - 1.0_dp )
          else if ( le(x(i,j),1.0_dp) ) then
             g = 0.8_dp + ( x(i,j) - 1.0_dp ) / ( 1.0_dp - alpha(i) )
          else
             g = x(i,j) - 1.0_dp
          end if

          deceptive_2d(j) = deceptive_2d(j) + g

       end do

       deceptive_2d(j) = deceptive_2d(j) / real( m, dp )
       deceptive_2d(j) = -( deceptive_2d(j) ** beta )

    end do

  end function deceptive_2d

END MODULE mo_opt_functions
