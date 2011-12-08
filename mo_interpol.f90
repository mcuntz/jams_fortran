MODULE mo_interpol

  ! This module provides a linear interpolation routine for irregular grids

  ! Literature
  !   Inspired by interpolation routines of Numerical Recipes
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996
  !   and IDL routine interpol.pro Copyright (C) 1982-2004, Research Systems, Inc.

  ! Written March 2011, Matthias Cuntz

  USE mo_kind, ONLY: i4, sp, dp

  Implicit NONE

  PRIVATE

  PUBLIC :: interpol        ! Linearly interpolate vectors with an irregular grid

  INTERFACE interpol
     MODULE PROCEDURE interpol_sp, interpol_dp
  END INTERFACE interpol

  ! Private function, from numerical recipes
  INTERFACE locate
     MODULE PROCEDURE locate_0d_sp, locate_1d_sp, locate_0d_dp, locate_1d_dp
  END INTERFACE locate

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         interpol

  !     PURPOSE
  !         Linearly interpolate vectors with an irregular grid.
  !            yout = interpol(yin, xin, xout)
  !            be z = xin(i) then
  !            yout(i) = yin(i) + (z - floor(z)) * (yin(i+1) - yin(i))

  !     CALLING SEQUENCE
  !         yout = interpol(yin, xin, xout)
  
  !     INDENT(IN)
  !         real(sp/dp) :: yin(:)      1D-array with input values to interpolate
  !         real(sp/dp) :: xin(:)      1D-array with position of yin values
  !         real(sp/dp) :: xout(:)     1D-array with wanted positions of interpolated values

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp) :: yout        yin interpolated at xout

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Index arrays are not working with all compilers. Compile with -DABSOFT if problems occur.

  !     EXAMPLE
  !         ! Fill masked values of vector w with interpolated values
  !         w = interpol(pack(w,wmaske), pack(date,wmaske), date)
  !         -> see also example in test directory

  !     LITERATURE
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996
  !       IDL routine interpol.pro Copyright (C) 1982-2004, Research Systems, Inc.

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2011

  FUNCTION interpol_dp(v, x, u)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:), INTENT(IN) :: v
    REAL(dp),    DIMENSION(:), INTENT(IN) :: x
    REAL(dp),    DIMENSION(:), INTENT(IN) :: u
    REAL(dp),    DIMENSION(size(u))       :: interpol_dp

    INTEGER(i4), DIMENSION(size(u)) :: s
#ifdef ABSOFT
    REAL(dp)                        :: ums
#else
    REAL(dp),    DIMENSION(size(u)) :: ums
#endif
    INTEGER(i4) :: n, m
#ifdef ABSOFT
    INTEGER(i4) :: i
#endif

    m = size(v)
    if (size(x) /= m) stop 'Error interpol_dp: size(v) /= size(x)'
    n = size(u)

    s = min(max(locate(x, u), 1), m-1)  ! Subscript intervals

#ifdef ABSOFT
    do i=1, n
       ums = u(i) - x(s(i))
       if (abs(ums) < epsilon(1.0_dp)) ums = 0.0_dp
       interpol_dp(i) = v(s(i)) + ums*(v(s(i)+1)-v(s(i)))/(x(s(i)+1)-x(s(i)))
    end do
#else
    ums = (u-x(s)) ! distance from point before
    where (abs(ums) < epsilon(1.0_dp)) ums = 0.0_dp ! for numerical stability
    interpol_dp = v(s) + ums*(v(s+1)-v(s))/(x(s+1)-x(s))
#endif

  END FUNCTION interpol_dp


  FUNCTION interpol_sp(v, x, u)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:), INTENT(IN) :: v
    REAL(sp),    DIMENSION(:), INTENT(IN) :: x
    REAL(sp),    DIMENSION(:), INTENT(IN) :: u
    REAL(sp),    DIMENSION(size(u))       :: interpol_sp

    INTEGER(i4), DIMENSION(size(u)) :: s
#ifdef ABSOFT
    REAL(sp)                        :: ums
#else
    REAL(sp),    DIMENSION(size(u)) :: ums
#endif
    INTEGER(i4) :: n, m
#ifdef ABSOFT
    INTEGER(i4) :: i
#endif

    m = size(v)
    if (size(x) /= m) stop 'Error interpol_sp: size(v) /= size(x)'
    n = size(u)

    s = min(max(locate(x, u), 1), m-1)  ! Subscript intervals

#ifdef ABSOFT
    do i=1, n
       ums = u(i) - x(s(i))
       if (abs(ums) < epsilon(1.0_sp)) ums = 0.0_sp
       interpol_sp(i) = v(s(i)) + ums*(v(s(i)+1)-v(s(i)))/(x(s(i)+1)-x(s(i)))
    end do
#else
    ums = (u-x(s)) ! distance from point before
    where (abs(ums) < epsilon(1.0_sp)) ums = 0.0_sp ! for numerical stability
    interpol_sp = v(s) + ums*(v(s+1)-v(s))/(x(s+1)-x(s))
#endif

  END FUNCTION interpol_sp

  ! ------------------------------------------------------------------

  ! Given an array xx(1:N), and given a value x, returns a value j such that x is between 
  !  xx(j) and xx(j+1). xx must be monotonic, either increasing or decreasing. j=0 or 
  ! j=N is returned to indicate that x is out of range.
  FUNCTION locate_0d_dp(xx,x)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: xx
    REAL(dp),               INTENT(IN) :: x
    INTEGER(i4)                        :: locate_0d_dp

    INTEGER(i4) :: n, jl, jm, ju
    LOGICAL     :: ascnd

    n = size(xx)
    ascnd = (xx(n) >= xx(1))
    jl = 0
    ju = n+1
    do
       if (ju-jl <= 1) exit
       jm = (ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       end if
    end do
    if (x == xx(1)) then
       locate_0d_dp = 1
    else if (x == xx(n)) then
       locate_0d_dp = n-1
    else
       locate_0d_dp = jl
    end if

  END FUNCTION locate_0d_dp


  FUNCTION locate_1d_dp(xx,x)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: xx
    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), DIMENSION(size(x))    :: locate_1d_dp

    INTEGER(i4) :: n, jl, jm, ju
    LOGICAL     :: ascnd
    INTEGER(i4) :: nx, i

    n = size(xx)
    ascnd = (xx(n) >= xx(1))

    nx = size(x)
    do i=1, nx
       jl = 0
       ju = n+1
       do
          if (ju-jl <= 1) exit
          jm = (ju+jl)/2
          if (ascnd .eqv. (x(i) >= xx(jm))) then
             jl = jm
          else
             ju = jm
          end if
       end do
       if (x(i) == xx(1)) then
          locate_1d_dp(i) = 1
       else if (x(i) == xx(n)) then
          locate_1d_dp(i) = n-1
       else
          locate_1d_dp(i) = jl
       end if
    end do

  END FUNCTION locate_1d_dp


  FUNCTION locate_0d_sp(xx,x)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: xx
    REAL(sp),               INTENT(IN) :: x
    INTEGER(i4)                        :: locate_0d_sp

    INTEGER(i4) :: n, jl, jm, ju
    LOGICAL     :: ascnd

    n = size(xx)
    ascnd = (xx(n) >= xx(1))
    jl = 0
    ju = n+1
    do
       if (ju-jl <= 1) exit
       jm = (ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       end if
    end do
    if (x == xx(1)) then
       locate_0d_sp = 1
    else if (x == xx(n)) then
       locate_0d_sp = n-1
    else
       locate_0d_sp = jl
    end if

  END FUNCTION locate_0d_sp


  FUNCTION locate_1d_sp(xx,x)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: xx
    REAL(sp), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), DIMENSION(size(x))    :: locate_1d_sp

    INTEGER(i4) :: n, jl, jm, ju
    LOGICAL     :: ascnd
    INTEGER(i4) :: nx, i

    n = size(xx)
    ascnd = (xx(n) >= xx(1))

    nx = size(x)
    do i=1, nx
       jl = 0
       ju = n+1
       do
          if (ju-jl <= 1) exit
          jm = (ju+jl)/2
          if (ascnd .eqv. (x(i) >= xx(jm))) then
             jl = jm
          else
             ju = jm
          end if
       end do
       if (x(i) == xx(1)) then
          locate_1d_sp(i) = 1
       else if (x(i) == xx(n)) then
          locate_1d_sp(i) = n-1
       else
          locate_1d_sp(i) = jl
       end if
    end do

  END FUNCTION locate_1d_sp

  ! ------------------------------------------------------------------

END MODULE mo_interpol
